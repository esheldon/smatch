#include <Python.h>
#include <numpy/arrayobject.h> 

#include "math.h"
#include "defs.h"
#include "vector.h"
#include "tree.h"
#include "healpix.h"

struct PySMatchCat {
    PyObject_HEAD
    int64_t maxmatch;

    point_vector *pts;
    struct healpix* hpix;
    struct tree_node* tree;

    // we keep this separately, for the case of writing
    // matches to a file
    int64_t nmatches;

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
        PyErr_SetString(PyExc_ValueError, "Entered ra/dec must have size > 0\n");
        return NULL;
    }
    nrad = PyArray_SIZE(radObj);
    if (nrad != n && nrad != 1) {
        PyErr_Format(PyExc_ValueError, 
                     "radii must be size 1 or same length as ra,dec (%ld).  Got %ld\n",n,nrad);
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
    int64_t hpix_id=0;
    int64_t half_npix=0;
    size_t index=0, ihpix=0;
    
    listpix = lvector_new();

    // this will produce a more balanced tree across the whole sky
    half_npix = hpix->npix/2;

    pt=pts->data;
    for (index=0; index < pts->size; index++) {

        pt = &pts->data[index];
        hpix_disc_intersect(hpix, pt->x, pt->y, pt->z, pt->radius, listpix);

        for (ihpix=0; ihpix < listpix->size; ihpix++) {
            hpix_id = listpix->data[ihpix];
            tree_insert(&tree, hpix_id-half_npix, index);
        }
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

    self->hpix = hpix_new((int64_t)nside);
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

    vector_free(self->pts);
    self->hpix = hpix_delete(self->hpix);
    self->tree = tree_delete(self->tree);

#if PY_MAJOR_VERSION >= 3
    Py_TYPE(self)->tp_free((PyObject*)self);
#else
    self->ob_type->tp_free((PyObject*)self);
#endif

}


static PyObject *
PySMatchCat_repr(struct PySMatchCat* self) {
    char repr[256];
    sprintf(repr, "Catalog\n    hpix nside: %ld", self->hpix->nside);
#if PY_MAJOR_VERSION >= 3
    // bytes
    return Py_BuildValue("y",repr);
#else
    return Py_BuildValue("s",repr);
#endif
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
    return Py_BuildValue("l", self->nmatches);
}

static PyObject* PySMatchCat_set_nmatches(struct PySMatchCat* self, PyObject *args)
{
    PY_LONG_LONG nmatches=0;

    if (!PyArg_ParseTuple(args, (char*)"L", &nmatches)) {
        return NULL;
    }

    self->nmatches=nmatches;
    Py_RETURN_NONE;
}


static void domatch1(const struct PySMatchCat* self, 
                     double ra,
                     double dec,
                     size_t input_ind,
                     match_vector* matches) {

    int64_t hpixid=0;
    struct tree_node* node=NULL;
    int64_t half_npix=0;

    size_t i=0;
    int64_t cat_ind=0;
    double x=0,y=0,z=0;
    double cos_radius=0, cos_angle=0;

    Match match={0};

    vector_resize(matches,0);

    half_npix = self->hpix->npix/2;
    hpixid = hpix_eq2pix(self->hpix, ra, dec) - half_npix;
    node = tree_find(self->tree, hpixid);

    if (node != NULL) {
        Point* pt=NULL;
        hpix_eq2xyz(ra,dec,&x,&y,&z);

        for (i=0; i<node->indices->size; i++) {
            // index into other list
            cat_ind = node->indices->data[i];

            pt = &self->pts->data[cat_ind];

            cos_radius = cos(pt->radius);
            cos_angle = pt->x*x + pt->y*y + pt->z*z;

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
    self->nmatches=0;

    n = PyArray_SIZE(raObj);
    for (i=0; i<n ; i++) {

        raptr=PyArray_GETPTR1(raObj, i);
        decptr=PyArray_GETPTR1(decObj, i);

        domatch1(self, *raptr, *decptr, i, new_matches);

        if (new_matches->size == 0) {
            continue;
        }

        push_matches(self->matches, new_matches);
        self->nmatches += new_matches->size;

    }
    vector_free(new_matches);

    return 1;
}

/*

   do the matching and fill in the self->matches array of structures

*/
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

static int domatch2file(struct PySMatchCat* self,
                        PyObject* raObj,
                        PyObject* decObj,
                        const char* filename) {
    size_t i=0, j=0, n=0;
    double *raptr=NULL, *decptr=NULL;
    FILE* fobj=NULL;
    const Match *m=NULL;
    match_vector* matches = match_vector_new();

    // always reset the match structure, even though
    // we are writing to a file
    vector_clear(self->matches);
    self->nmatches=0;

    fobj=fopen(filename, "w");
    if (fobj == NULL) {
        return 0;
    }

    n = PyArray_SIZE(raObj);
    for (i=0; i<n ; i++) {

        raptr=PyArray_GETPTR1(raObj, i);
        decptr=PyArray_GETPTR1(decObj, i);

        domatch1(self, *raptr, *decptr, i, matches);

        for (j=0; j<matches->size; j++) {
            m=&matches->data[j];
            fprintf(fobj, "%ld %ld %.16g\n", m->cat_ind, m->input_ind, m->cosdist);
        }
        self->nmatches += matches->size;

    }
    vector_free(matches);

    fclose(fobj);

    return 1;
}



/*

   do the matching and write to the indicated file

*/
static PyObject* PySMatchCat_match2file(struct PySMatchCat* self, PyObject *args)
{
    PyObject* raObj=NULL;
    PyObject* decObj=NULL;
    PY_LONG_LONG maxmatch=0;
    const char *filename=NULL;

    if (!PyArg_ParseTuple(args, (char*)"LOOs", &maxmatch, &raObj, &decObj, &filename)) {
        return NULL;
    }

    self->maxmatch=maxmatch;

    domatch2file(self, raObj, decObj, filename);
    Py_RETURN_NONE;
}

int64_t count_lines(FILE* fobj)
{

    char ch;
	int64_t nlines=0;

	while(!feof(fobj))
	{
		ch = fgetc(fobj);
		if(ch == '\n')
		{
			nlines++;
		}
	}

	return nlines;
}

static PyObject *
PySMatchCat_count_lines(PyObject* self, PyObject* args)
{
    const char *filename=NULL;
    FILE* fobj=NULL;
    int64_t nlines=0;

    if (!PyArg_ParseTuple(args, (char*)"s", &filename)) {
        return NULL;
    }

    fobj=fopen(filename, "r");
    if (fobj == NULL) {
        PyErr_Format(PyExc_IOError, "Could not open file: '%s'\n", filename);
        return NULL;
    }

    nlines = count_lines(fobj);

    fclose(fobj);

    return Py_BuildValue("l", nlines);

}

/*

   read from a match file

   i1 i2 cosdist

   i1->catalog index
   i2->secondary or input index for the match2file routine
   cosdist-> cos(angular distance)

*/

static PyObject* PySMatchCat_load_matches(PyObject* self, PyObject *args)
{
    const char *filename=NULL;
    FILE* fobj=NULL;
    int nread=0, ncols=3;
    npy_intp nmatches=0, i=0;
    Match *match=NULL;

    PyObject *matchesObj=NULL;

    if (!PyArg_ParseTuple(args, (char*)"sO", &filename, &matchesObj)) {
        return NULL;
    }

    fobj=fopen(filename, "r");
    if (fobj == NULL) {
        PyErr_Format(PyExc_IOError, "Could not open file: '%s'\n", filename);
        return NULL;
    }

    nmatches = PyArray_SIZE(matchesObj);

    for (i=0; i<nmatches; i++) {
        match = (Match* ) PyArray_GETPTR1(matchesObj, i);
        nread=fscanf(fobj,
                     "%ld %ld %lf\n", 
                     &match->cat_ind,
                     &match->input_ind,
                     &match->cosdist);
        if (nread != ncols) {
            PyErr_Format(PyExc_IOError,
                         "Error: only read %d at line %ld of file: '%s'\n", nread,
                         i+1,
                         filename);
            goto _load_matches_bail;
        }
    }

_load_matches_bail:
    fclose(fobj);
    if (nread != ncols) {
        return NULL;
    }

    Py_RETURN_NONE;
}



/*
   Copy into the given structured array.

   No error checking is performed here, set up the data in python

*/
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
    {"_set_nmatches",           (PyCFunction)PySMatchCat_set_nmatches,          METH_VARARGS,  "Set the number of matches, useful when reading from a file."},
    {"get_hpix_nside",              (PyCFunction)PySMatchCat_hpix_nside,          METH_VARARGS,  "Get the nside for healpix."},
    {"get_hpix_area",              (PyCFunction)PySMatchCat_hpix_area,          METH_VARARGS,  "Get the nside for healpix."},
    {"match",              (PyCFunction)PySMatchCat_match,          METH_VARARGS,  "Match the catalog to the input ra,dec arrays."},
    {"match2file",              (PyCFunction)PySMatchCat_match2file,          METH_VARARGS,  "Match the catalog to the input ra,dec arrays and write results to a file."},
    {"_copy_matches",              (PyCFunction)PySMatchCat_copy_matches,          METH_VARARGS,  "Copy the matches into the input array."},
    {NULL}  /* Sentinel */
};



static PyTypeObject PyCatalogType = {
#if PY_MAJOR_VERSION >= 3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
#endif
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







static PyMethodDef smatch_module_methods[] = {
    {"_count_lines",      (PyCFunction)PySMatchCat_count_lines, METH_VARARGS,  "count the lines in the specified file."},
    {"_load_matches",              (PyCFunction)PySMatchCat_load_matches,          METH_VARARGS,  "Load matches from the specifed filename."},
};


#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_smatch",      /* m_name */
        "Defines the catalog class and some methods",  /* m_doc */
        -1,                  /* m_size */
        smatch_module_methods, /* m_methods */
        NULL,                /* m_reload */
        NULL,                /* m_traverse */
        NULL,                /* m_clear */
        NULL,                /* m_free */
    };
#endif



#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
PyInit__smatch(void) 
#else
init_smatch(void) 
#endif
{
    PyObject* m;

    PyCatalogType.tp_new = PyType_GenericNew;


#if PY_MAJOR_VERSION >= 3
    if (PyType_Ready(&PyCatalogType) < 0) {
        return NULL;
    }
    m = PyModule_Create(&moduledef);
    if (m==NULL) {
        return NULL;
    }

#else

    if (PyType_Ready(&PyCatalogType) < 0)
        return;

    m = Py_InitModule3("_smatch", smatch_module_methods, "Define module methods.");
    if (m==NULL) {
        return;
    }
#endif

    Py_INCREF(&PyCatalogType);
    PyModule_AddObject(m, "Catalog", (PyObject *)&PyCatalogType);

    import_array();

#if PY_MAJOR_VERSION >= 3
    return m;
#endif
}
