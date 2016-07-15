#include <Python.h>
#include <numpy/arrayobject.h> 

#include "math.h"
#include "defs.h"
#include "vector.h"
#include "tree.h"
#include "healpix.h"

struct PySMatchCat {
    PyObject_HEAD
    int64 maxmatch;

    point_vector *pts;
    struct healpix* hpix;
    struct tree_node* tree;

    match_vector *matches;

};

// sort functions for the matches
static int match_compare(const void *a, const void *b) {
    // we want to sort largest first, so will
    // reverse the normal trend
    double temp = 
        ((Match*)b)->cosdist
         -
        ((Match*)a)->cosdist;
    if (temp > 0)
        return 1;
    else if (temp < 0)
        return -1;
    else
        return 0;
}

static void match_vector_sort(match_vector* self) {
    qsort(self->data, self->size, sizeof(Match), match_compare);
}

static void push_matches(match_vector* self, const match_vector* matches)
{
    size_t i=0;
    const Match* match=NULL;

    for (i=0; i<matches->size; i++) {
        match=&matches->data[i];

        vector_push(self, *match);
    }
    
}



/*
 * We assume these are in degrees, double precision, and are numpy arrays of
 * same length
 *
 * radii get converted to radians
 */

static point_vector* points_init(PyObject* raObj, PyObject* decObj, PyObject* radObj) {
    point_vector *pts=NULL;
    Point* pt = NULL;
    double *raptr=NULL, *decptr=NULL, *radptr=NULL;

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


    pts = point_vector_new();

    vector_resize(pts, n); 

    if (nrad == 1) {
        radptr = PyArray_GETPTR1(radObj, 0);
        pt->radius = (*radptr)*D2R;
    }
    for (i=0; i<n; i++) {

        pt=&pts->data[i];

        raptr=PyArray_GETPTR1(raObj, i);
        decptr=PyArray_GETPTR1(decObj, i);

        hpix_eq2xyz(*raptr, *decptr, &pt->x, &pt->y, &pt->z);

        if (nrad > 1) {
            radptr = PyArray_GETPTR1(radObj, i);
            pt->radius = (*radptr)*D2R;
        }

    }

    return pts;
}

// create a tree based on the healpix id
static struct tree_node* create_tree(struct healpix* hpix, point_vector* pts) {
    struct tree_node* tree=NULL;
    lvector* listpix=NULL;
    Point* pt=NULL;
    int64* data=NULL;
    int64 half_npix=0;
    
    listpix = lvector_new();

    // this will produce a more balanced tree across the whole sky
    half_npix = hpix->npix/2;

    size_t count=0;
    pt=pts->data;
    while (pt < pts->data + pts->size) {
        hpix_disc_intersect(hpix, pt->x, pt->y, pt->z, pt->radius, listpix);

        data=listpix->data;
        while (data < listpix->data + listpix->size) {
            tree_insert(&tree, (*data)-half_npix, count);
            data++;
        }

        pt++;
        count++;
    }
    vector_free(listpix);

    return tree;
}


static int
PySMatchCat_init(struct PySMatchCat* self, PyObject *args, PyObject *kwds)
{
    PY_LONG_LONG nside=0;
    PyObject* raObj=NULL;
    PyObject* decObj=NULL;
    PyObject* radObj=NULL;
    int err=0;

    if (!PyArg_ParseTuple(args, (char*)"LOOO", &nside, &raObj, &decObj, &radObj)) {
        return -1;
    }

    self->tree=NULL;
    self->hpix=NULL;
    self->pts=NULL;

    self->matches = match_vector_new();

    self->hpix = hpix_new((int64)nside);
    if (self->hpix==NULL) {
        err=1;
        goto _catalog_init_cleanup;
    }

    self->pts = points_init(raObj, decObj, radObj);
    if (self->pts==NULL) {
        err=1;
        goto _catalog_init_cleanup;
    }

    self->tree = create_tree(self->hpix, self->pts);
    if (self->tree==NULL) {
        err=1;
        goto _catalog_init_cleanup;
    }

_catalog_init_cleanup:
    if (err != 0) {
        vector_free(self->pts);
        self->hpix = hpix_delete(self->hpix);
        self->tree = tree_delete(self->tree);
        return -1;
    }
    return 0;
}


static void
PySMatchCat_dealloc(struct PySMatchCat* self)
{
    //self->tree = tree_delete(self->tree);
    vector_free(self->pts);
    self->hpix = hpix_delete(self->hpix);
    self->tree = tree_delete(self->tree);
    self->ob_type->tp_free((PyObject*)self);
}


static PyObject *
PySMatchCat_repr(struct PySMatchCat* self) {
    char repr[256];
    sprintf(repr, "Catalog\n    hpix nside: %ld", self->hpix->nside);
    return PyString_FromString(repr);
}



static PyObject *
PySMatchCat_hpix_nside(struct PySMatchCat* self) {
    return Py_BuildValue("l", self->hpix->nside);
}
static PyObject *
PySMatchCat_hpix_area(struct PySMatchCat* self) {
    return Py_BuildValue("d", hpix_area(self->hpix->nside));
}


static PyObject *
PySMatchCat_nmatches(struct PySMatchCat* self) {
    return Py_BuildValue("l", self->matches->size);
}


static void domatch1(const struct PySMatchCat* self, 
                     double ra,
                     double dec,
                     size_t input_ind,
                     match_vector* matches) {

    int64 hpixid=0;
    struct tree_node* node=NULL;
    int64 half_npix=0;
    Match match={0};

    vector_resize(matches,0);

    half_npix = self->hpix->npix/2;
    hpixid = hpix_eq2pix(self->hpix, ra, dec) - half_npix;
    node = tree_find(self->tree, hpixid);

    if (node != NULL) {
        Point* pt=NULL;
        size_t i=0;
        int64_t cat_ind=0;
        double x=0,y=0,z=0;
        double cos_radius=0;
        hpix_eq2xyz(ra,dec,&x,&y,&z);

        for (i=0; i<node->indices->size; i++) {
            // index into other list
            cat_ind = node->indices->data[i];

            pt = &self->pts->data[cat_ind];

            cos_radius = cos(pt->radius);
            double cos_angle = pt->x*x + pt->y*y + pt->z*z;

            if (cos_angle > cos_radius) {
                match.cat_ind=cat_ind;
                match.input_ind=(int64_t)input_ind;
                match.cosdist=cos_angle;
                vector_push(matches, match);
            }
        }
    }

    if (self->maxmatch > 0 && self->maxmatch < matches->size) {
        // max match count was given
        // If we have too many matches, sort closest first and take
        // the closest maxmatch matches
        match_vector_sort(matches);
        vector_resize(matches, self->maxmatch);
    }

}

static int domatch(struct PySMatchCat* self, PyObject* raObj, PyObject* decObj) {
    size_t i=0, n=0;
    double *raptr=NULL, *decptr=NULL;
    match_vector* new_matches = match_vector_new();

    // always reset the match structure
    vector_clear(self->matches);

    n = PyArray_SIZE(raObj);
    for (i=0; i<n ; i++) {

        raptr=PyArray_GETPTR1(raObj, i);
        decptr=PyArray_GETPTR1(decObj, i);

        domatch1(self, *raptr, *decptr, i, new_matches);

        if (new_matches->size == 0) {
            continue;
        }

        push_matches(self->matches, new_matches);

    }
    vector_free(new_matches);

    return 1;
}

static PyObject* PySMatchCat_match(struct PySMatchCat* self, PyObject *args)
{
    PyObject* raObj=NULL;
    PyObject* decObj=NULL;
    PY_LONG_LONG maxmatch=0;

    if (!PyArg_ParseTuple(args, (char*)"LOO", &maxmatch,&raObj, &decObj)) {
        return NULL;
    }

    self->maxmatch=maxmatch;

    domatch(self, raObj, decObj);
    Py_RETURN_NONE;
}

static PyObject* PySMatchCat_copy_matches(struct PySMatchCat* self, PyObject *args)
{
    PyObject* matchesObj=NULL;
    Match* matches=NULL;

    if (!PyArg_ParseTuple(args, (char*)"O", &matchesObj)) {
        return NULL;
    }

    matches = PyArray_DATA(matchesObj);

    memmove(matches,
            self->matches->data,
            self->matches->size*sizeof(Match)
    );

    Py_RETURN_NONE;
}





static PyMethodDef PySMatchCat_methods[] = {
    {"get_nmatches",           (PyCFunction)PySMatchCat_nmatches,          METH_VARARGS,  "Get the number of matches."},
    {"get_hpix_nside",              (PyCFunction)PySMatchCat_hpix_nside,          METH_VARARGS,  "Get the nside for healpix."},
    {"get_hpix_area",              (PyCFunction)PySMatchCat_hpix_area,          METH_VARARGS,  "Get the nside for healpix."},
    {"match",              (PyCFunction)PySMatchCat_match,          METH_VARARGS,  
        "Match the catalog to the input ra,dec arrays."},
    {"_copy_matches",              (PyCFunction)PySMatchCat_copy_matches,          METH_VARARGS,  "Copy the matches into the input array."},
    {NULL}  /* Sentinel */
};



static PyTypeObject PyCatalogType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "_smatch.Catalog",             /*tp_name*/
    sizeof(struct PySMatchCat), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)PySMatchCat_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    //0,                         /*tp_repr*/
    (reprfunc)PySMatchCat_repr,                         /*tp_repr*/
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
    PySMatchCat_methods,             /* tp_methods */
    0,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    //0,     /* tp_init */
    (initproc)PySMatchCat_init,      /* tp_init */
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
