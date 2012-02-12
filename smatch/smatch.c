#include <Python.h>
#include <numpy/arrayobject.h> 

#include "math.h"
#include "defs.h"
#include "stack.h"
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



struct tree_node* create_tree(struct healpix* hpix, struct points* pts) {
    /*
    size_t i=0;

    struct point* pt=&pts->data[0];
    for (i=0; i< pts->size; i++) {
        hpix_disc_intersect(hpix, pt->x, pt->y, pt->z, pt->rad, listpix);

        pt++;
    }
    */
    return NULL;

}


static int
PyCatalogObject_init(struct PyCatalogObject* self, PyObject *args, PyObject *kwds)
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
    /*
    self->tree = create_tree(self->hpix, self->pts);
    if (self->tree==NULL) {
        err=1;
        goto _catalog_init_cleanup;
    }
    */

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

static PyMethodDef PyCatalogObject_methods[] = {
    {"nside",              (PyCFunction)PyCatalogObject_nside,          METH_VARARGS,  "nside\n\nReturn the nside for healpix."},
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
