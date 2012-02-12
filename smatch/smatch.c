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
        fprintf(stderr,"%lf %lf %lf %lf %lf %e\n", *ra, *dec, pt->x, pt->y, pt->z, pt->rad);
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

    fprintf(stderr,"Initializing healpix struct, nside: %lld\n", nside);
    self->hpix = hpix_new((int64)nside);
    if (self->hpix==NULL) {
        err=1;
        goto _catalog_init_cleanup;
    }

    fprintf(stderr,"Initializing x,y,z points\n");
    self->pts = points_init(raObj, decObj, radObj);
    if (self->pts==NULL) {
        err=1;
        goto _catalog_init_cleanup;
    }
    fprintf(stderr,"Initializing tree\n");
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

PyObject* domatch(struct PyCatalogObject* self, PyObject* raObj, PyObject* decObj) {
    size_t i=0, n=0;
    double *raptr=NULL, *decptr=NULL;

    PyObject* cat_indObj=NULL;
    PyObject* input_indObj=NULL;
    npy_intp dims[1];
    int dtype=NPY_INTP;

    PyObject* resTuple=NULL;

    struct matchstack* matches = matchstack_new();

    // these are outputs
    struct szstack* cat_ind   = szstack_new(0);
    struct szstack* input_ind = szstack_new(0);


    n = PyArray_SIZE(raObj);
    for (i=0; i<n ; i++) {
        raptr=PyArray_GETPTR1(raObj, i);
        decptr=PyArray_GETPTR1(decObj, i);
        domatch1(self, *raptr, *decptr, i, matches);

        if (self->maxmatch > 0) {
            // an exact max allowed was given
        } else {
            // all matches are kept
        }
    }
    matches=matchstack_delete(matches);

    if (matches->size == 0) {
        dims[0] = 0;
        cat_indObj = PyArray_ZEROS(1, dims, dtype, 0);
        input_indObj = PyArray_ZEROS(1, dims, dtype, 0);
    } else {
        npy_intp* captr=NULL;
        npy_intp* iaptr=NULL;
        size_t* cptr=NULL;
        size_t* iptr=NULL;
        dims[0] = cat_ind->size;
        cat_indObj = PyArray_EMPTY(1, dims, dtype, 0);
        input_indObj = PyArray_EMPTY(1, dims, dtype, 0);

        captr=PyArray_DATA(cat_indObj);
        iaptr=PyArray_DATA(input_indObj);
        cptr=cat_ind->data;
        iptr=input_ind->data;
        for (i=0; i<dims[0]; i++) {

            *captr = (npy_intp) (*cptr);
            *iaptr = (npy_intp) (*iptr);

            captr++;
            iaptr++;
            cptr++;
            iptr++;
        }
    }

    resTuple=PyTuple_New(2);
    PyTuple_SetItem(resTuple, 0, cat_indObj);
    PyTuple_SetItem(resTuple, 1, input_indObj);
    return resTuple;
}

PyObject* PyCatalogObject_match(struct PyCatalogObject* self, PyObject *args)
{
    PyObject* raObj=NULL;
    PyObject* decObj=NULL;
    PyObject* resObj=NULL;

    if (!PyArg_ParseTuple(args, (char*)"OO", &raObj, &decObj)) {
        Py_XINCREF(Py_None);
        return Py_None;
    }

    resObj = domatch(self, raObj, decObj);
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
