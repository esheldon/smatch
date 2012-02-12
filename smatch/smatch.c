// todo
//   store cos(rad) instead of rad in point?
//   make a npy_intp stack, and just re-use the memory
//   factor the domatch program
#include <Python.h>
#include <numpy/arrayobject.h> 

#include "math.h"
#include "defs.h"
#include "stack.h"
#include "match.h"
#include "tree.h"
#include "healpix.h"

struct point {
    double x;
    double y;
    double z;
    double rad; // radians
};

struct points {
    size_t size;
    struct point* data;
};

struct PyCatalogObject {
    PyObject_HEAD
    int64 maxmatch;
    struct points* pts;
    struct healpix* hpix;
    struct tree_node* tree;
};


/*
 * We assume these are in degrees, double precision, and are numpy arrays of
 * same length
 *
 * radii get converted to radians
 */

struct points* points_init(PyObject* raObj, PyObject* decObj, PyObject* radObj) {
    struct points* pts = NULL;
    struct point* pt = NULL;

    npy_intp n=0, i=0, nrad=0;
    n = PyArray_SIZE(raObj);
    if (n <= 0) {
        PyErr_SetString(PyExc_ValueError, "Entered lon/lat must have size > 0\n");
        return NULL;
    }
    nrad = PyArray_SIZE(radObj);
    if (nrad != n && nrad != 1) {
        PyErr_Format(PyExc_ValueError, 
                     "radii must be size 1 or same length as ra,ded (%ld).  Got %ld\n",n,nrad);
        return NULL;
    }


    pts = calloc(1,sizeof(struct points));
    if (pts == NULL) {
        PyErr_Format(PyExc_MemoryError, "Could not allocate struct points\n");
        return NULL;
    }
    pts->size = n;
    pts->data = calloc(n,sizeof(struct point));
    if (pts->data == NULL) {
        PyErr_Format(PyExc_MemoryError, "Could not allocate points\n");
        free(pts);
        return NULL;
    }

    pt=&pts->data[0];
    double *ra=NULL, *dec=NULL, *rad=NULL;
    for (i=0; i<n; i++) {
        ra=PyArray_GETPTR1(raObj, i);
        dec=PyArray_GETPTR1(decObj, i);

        hpix_eq2xyz(*ra, *dec, &pt->x, &pt->y, &pt->z);

        if (nrad > 1) {
            rad = PyArray_GETPTR1(radObj, i);
            *rad = (*rad)*D2R;
        } else if (i == 0) {
            rad = PyArray_GETPTR1(radObj, 0);
            *rad = (*rad)*D2R;
        }

        pt->rad = *rad;
        //fprintf(stderr,"%lf %lf %lf %lf %lf %e\n", *ra, *dec, pt->x, pt->y, pt->z, pt->rad);
        pt++;
    }



    return pts;
}



// create a tree based on the healpix id
struct tree_node* create_tree(struct healpix* hpix, struct points* pts) {
    struct tree_node* tree=NULL;
    struct i64stack* listpix=NULL;
    struct point* pt=NULL;
    int64* data=NULL;
    int64 half_npix=0;
    
    listpix = i64stack_new(0);

    // this will produce a more balanced tree across the whole sky
    half_npix = hpix->npix/2;

    size_t count=0;
    pt=pts->data;
    while (pt < pts->data + pts->size) {
        hpix_disc_intersect(hpix, pt->x, pt->y, pt->z, pt->rad, listpix);

        data=listpix->data;
        while (data < listpix->data + listpix->size) {
            tree_insert(&tree, (*data)-half_npix, count);
            data++;
        }

        pt++;
        count++;
    }
    listpix=i64stack_delete(listpix);

    return tree;
}


static int
PyCatalogObject_init(struct PyCatalogObject* self, PyObject *args, PyObject *kwds)
{
    PY_LONG_LONG nside=0;
    PY_LONG_LONG maxmatch=0;
    PyObject* raObj=NULL;
    PyObject* decObj=NULL;
    PyObject* radObj=NULL;
    int err=0;

    if (!PyArg_ParseTuple(args, (char*)"LLOOO", &nside, &maxmatch, &raObj, &decObj, &radObj)) {
        return -1;
    }

    self->maxmatch=maxmatch;
    self->tree=NULL;
    self->hpix=NULL;
    self->pts=NULL;

    //fprintf(stderr,"Initializing healpix struct, nside: %lld\n", nside);
    self->hpix = hpix_new((int64)nside);
    if (self->hpix==NULL) {
        err=1;
        goto _catalog_init_cleanup;
    }

    //fprintf(stderr,"Initializing x,y,z points\n");
    self->pts = points_init(raObj, decObj, radObj);
    if (self->pts==NULL) {
        err=1;
        goto _catalog_init_cleanup;
    }
    //fprintf(stderr,"Initializing tree\n");
    self->tree = create_tree(self->hpix, self->pts);
    if (self->tree==NULL) {
        err=1;
        goto _catalog_init_cleanup;
    }

_catalog_init_cleanup:
    if (err != 0) {
        free(self->pts);
        self->hpix = hpix_delete(self->hpix);
        self->tree = tree_delete(self->tree);
        return -1;
    }
    return 0;
}


static void
PyCatalogObject_dealloc(struct PyCatalogObject* self)
{
    //self->tree = tree_delete(self->tree);
    if (self->pts != NULL) {
        free(self->pts->data);
        free(self->pts);
    }
    self->hpix = hpix_delete(self->hpix);
    self->tree = tree_delete(self->tree);
    self->ob_type->tp_free((PyObject*)self);
}


static PyObject *
PyCatalogObject_repr(struct PyCatalogObject* self) {
    char repr[256];
    sprintf(repr, "Catalog\n    hpix nside: %ld", self->hpix->nside);
    return PyString_FromString(repr);
}



static PyObject *
PyCatalogObject_nside(struct PyCatalogObject* self) {

    PyObject* nsideObj=NULL;
    nsideObj = PyLong_FromLongLong( (long long)self->hpix->nside);
    return nsideObj;
}

void SetOwnData(PyObject* array) {
    PyArrayObject* tmp=NULL;

    tmp = (PyArrayObject* ) array;
    tmp->flags |= NPY_OWNDATA;
}

PyObject* pack_results(struct i64stack* cat_ind,
                       struct i64stack* input_ind,
                       struct f64stack* rad,
                       int get_dist) {

    npy_intp dims[1];
    PyObject* cat_indObj=NULL;
    PyObject* input_indObj=NULL;
    PyObject* radObj=NULL;

    PyObject* resTuple=NULL;

    int ntup=2;
    if (get_dist) {
        ntup=3;
    }
    // make sure sizes are exact
    if (cat_ind->size != cat_ind->allocated_size) {
        i64stack_realloc(cat_ind, cat_ind->size);
        i64stack_realloc(input_ind, input_ind->size);
        if (get_dist) {
            f64stack_realloc(rad, rad->size);
        }
    }
    if (cat_ind->size == 0) {
        // no results found, delete the stacks and output
        // empty arrays
        dims[0] = 0;
        cat_indObj = PyArray_ZEROS(1, dims, NPY_INT64, 0);
        input_indObj = PyArray_ZEROS(1, dims, NPY_INT64, 0);
        cat_ind = i64stack_delete(cat_ind);
        input_ind = i64stack_delete(input_ind);
        if (get_dist) {
            radObj = PyArray_ZEROS(1, dims, NPY_FLOAT64, 0);
            rad=f64stack_delete(rad);
        }
    } else {

        // We use SimpleNewFromData and set the flags so it is owned
        // by the array.  We then do not free the stack *data* sections
        dims[0] = cat_ind->size;
        cat_indObj = PyArray_SimpleNewFromData(1, dims, NPY_INT64, cat_ind->data);
        SetOwnData(cat_indObj);
        input_indObj = PyArray_SimpleNewFromData(1, dims, NPY_INT64, input_ind->data);
        SetOwnData(input_indObj);

        // this only frees the structure, not the data at which it is pointing
        free(cat_ind); cat_ind=NULL;
        free(input_ind); input_ind=NULL;
        if (get_dist) {
            radObj = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, rad->data);
            SetOwnData(radObj);
            free(rad); rad=NULL;
        }

    }
    resTuple=PyTuple_New(ntup);
    PyTuple_SetItem(resTuple, 0, cat_indObj);
    PyTuple_SetItem(resTuple, 1, input_indObj);
    if (get_dist) {
        PyTuple_SetItem(resTuple, 2, radObj);
    }

    return resTuple;
}
void domatch1(struct PyCatalogObject* self, 
              double ra, double dec, size_t input_ind,
              struct matchstack* matches) {

    int64 hpixid=0;
    struct tree_node* node=NULL;
    int64 half_npix=0;


    matchstack_resize(matches,0);

    half_npix = self->hpix->npix/2;
    hpixid = hpix_eq2pix(self->hpix, ra, dec) - half_npix;
    node = tree_find(self->tree, hpixid);

    if (node != NULL) {
        struct point* pt=NULL;
        size_t i=0, cat_ind=0;
        double x=0,y=0,z=0;
        double cos_radius=0;
        hpix_eq2xyz(ra,dec,&x,&y,&z);

        for (i=0; i<node->indices->size; i++) {
            // index into other list
            cat_ind = node->indices->data[i];

            pt = &self->pts->data[cat_ind];

            cos_radius = cos(pt->rad);
            double cos_angle = pt->x*x + pt->y*y + pt->z*z;

            if (cos_angle > cos_radius) {
                matchstack_push(matches, cat_ind, input_ind, cos_angle);
            }
        }
    }
}

PyObject* domatch(struct PyCatalogObject* self, PyObject* raObj, PyObject* decObj, int get_dist) {
    size_t i=0, j=0, n=0, oldsize=0, newsize=0;
    double *raptr=NULL, *decptr=NULL, trad=0;
    struct match* match=NULL;

    //npy_intp dims[1];

    PyObject* resTuple=NULL;
    //int ntup=2;

    struct matchstack* matches = matchstack_new();

    // these are outputs.  Using int64 because that
    // is common to stacks and numpy arrays
    struct i64stack* cat_ind   = i64stack_new(0);
    struct i64stack* input_ind = i64stack_new(0);
    struct f64stack* rad = NULL;

    if (get_dist) {
        // don't allocate unless we need it
        rad = f64stack_new(0);
    }

    n = PyArray_SIZE(raObj);
    for (i=0; i<n ; i++) {

        raptr=PyArray_GETPTR1(raObj, i);
        decptr=PyArray_GETPTR1(decObj, i);

        domatch1(self, *raptr, *decptr, i, matches);

        if (matches->size == 0) {
            continue;
        }

        if (self->maxmatch > 0) {
            // max match count was given
            // If we have too many matches, sort biggest first and take
            // the closest maxmatch matches
            if (self->maxmatch < matches->size) {
                matchstack_sort(matches);
                matchstack_resize(matches, self->maxmatch);
            }
        }

        oldsize=cat_ind->size;
        newsize = oldsize + matches->size;
        i64stack_resize(cat_ind, newsize);
        i64stack_resize(input_ind, newsize);
        if (get_dist) {
            f64stack_resize(rad, newsize);
        }

        for (j=0; j<matches->size; j++) {
            match=&matches->data[j];

            cat_ind->data[oldsize+j] = match->cat_ind;
            input_ind->data[oldsize+j] = match->input_ind;
            if (get_dist) {
                trad = match->cosdist;
                if (trad >= 1) {
                    trad=0;
                } else {
                    trad = acos(trad)*R2D;
                }
                rad->data[oldsize+j] = trad;
            }
        }
    }
    matches=matchstack_delete(matches);

    resTuple = pack_results(cat_ind, input_ind, rad, get_dist);
    return resTuple;
}

PyObject* PyCatalogObject_match(struct PyCatalogObject* self, PyObject *args)
{
    PyObject* raObj=NULL;
    PyObject* decObj=NULL;
    int get_dist=0;
    PyObject* resObj=NULL;

    if (!PyArg_ParseTuple(args, (char*)"OOi", &raObj, &decObj, &get_dist)) {
        Py_XINCREF(Py_None);
        return Py_None;
    }

    resObj = domatch(self, raObj, decObj, get_dist);
    return resObj;
}





static PyMethodDef PyCatalogObject_methods[] = {
    {"nside",              (PyCFunction)PyCatalogObject_nside,          METH_VARARGS,  "nside\n\nReturn the nside for healpix."},
    {"match",              (PyCFunction)PyCatalogObject_match,          METH_VARARGS,  "match\n\nMatch the catalog to the input ra,dec arrays."},
    {NULL}  /* Sentinel */
};



static PyTypeObject PyCatalogType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "_smatch.Catalog",             /*tp_name*/
    sizeof(struct PyCatalogObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)PyCatalogObject_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    //0,                         /*tp_repr*/
    (reprfunc)PyCatalogObject_repr,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Catalog Class",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    PyCatalogObject_methods,             /* tp_methods */
    0,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    //0,     /* tp_init */
    (initproc)PyCatalogObject_init,      /* tp_init */
    0,                         /* tp_alloc */
    PyType_GenericNew,                 /* tp_new */
};







static PyMethodDef smatchertype_methods[] = {
    {NULL}  /* Sentinel */
};


#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
init_smatch(void) 
{
    PyObject* m;

    PyCatalogType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&PyCatalogType) < 0)
        return;

    m = Py_InitModule3("_smatch", smatchertype_methods, "Define smatcher type and methods.");

    Py_INCREF(&PyCatalogType);
    PyModule_AddObject(m, "Catalog", (PyObject *)&PyCatalogType);

    import_array();
}
