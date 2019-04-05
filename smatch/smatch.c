/*
   C code to perform matches on the sphere using healpix
*/

#include <Python.h>
#include <numpy/arrayobject.h> 

#include "math.h"
#include "defs.h"
#include "vector.h"
#include "tree.h"
#include "healpix.h"
#include "catpoint.h"
#include "cat.h"

struct PySMatchCat {
    PyObject_HEAD
    int64_t maxmatch;
    int matching_self;

    //Catalog *cat;
    struct healpix* hpix;

    // we keep this separately, for the case of writing
    // matches to a file
    int64_t nmatches;

};

typedef struct {
    npy_intp size;            // number of elements that are visible to the user
    npy_intp capacity;        // number of allocated elements in data vector
    PyObject *data;
} np_match_vector;

static inline void np_match_vector_realloc(np_match_vector* self, npy_intp newcap)
{
	npy_intp dm[1];
	PyArray_Dims dims;

	dm[0] = newcap;
	dims.ptr = dm;
	dims.len = 1;
    // returns NULL on failure, otherwise Py_None
	PyArray_Resize((PyArrayObject *)self->data, &dims, 0, NPY_CORDER);

    self->capacity = newcap;
    if (self->size > self->capacity) {
        self->size = self->capacity;
    }
}

static inline void np_match_vector_push(np_match_vector *self, Match* match)
{
    Match *ptr=NULL;
    npy_intp newcap=0;
    double tmp=0;

    if (self->size == self->capacity) {
        tmp = ceil(self->capacity*1.5);
        newcap = (npy_intp)tmp;
        if (newcap < 2) {
            newcap = 2;
        }
        np_match_vector_realloc(self, newcap);
    }

    self->size++;
    ptr = (Match *) PyArray_GETPTR1(self->data, self->size-1);
    *ptr = *match;
}

static inline void np_match_vector_push_many(np_match_vector *self,
                                             match_vector * matches)
{
    size_t i=0;

    for (i=0; i<vector_size(matches); i++) {
        np_match_vector_push(self, &matches->data[i]);
    }
}

//
// build a heap in an existing match vector
//

static inline void match_build_heap(match_vector* self)
{
    size_t i=0, c=0, parent;
    Match* data = self->data;
    Match tmp;

    if (vector_size(self) <= 1) {
        // Its already a heap
        return;
    }

    for (i=1; i<vector_size(self)-1; i++) {
        c=i;
		do {
			parent  = (c - 1) / 2;
			if  (data[parent].cosdist < data[c].cosdist)
			{
				tmp =  data[parent];
				data[parent] = data[c];
				data[c]  = tmp;
			}
			c =  parent;
		} while (c !=  0);
    }
}

//
//   make sure the heap is heapified after putting a new element
//   at the root
//
//   it is assumed the data after the root already form a heap
//

static void match_heap_fix(match_vector* self)
{
    size_t n = vector_size(self)-1;

    const Match* v = &self->data[0];
    Match* data = self->data;

    size_t jhi = 0;
    size_t jlo = 1;

    while (jlo <= n) {
        if (jlo < n && data[jlo].cosdist > data[jlo+1].cosdist) {
            // The right node is smaller
            jlo += 1;
        }
        if (v->cosdist <= data[jlo].cosdist) {
            // It forms a heap already
            break;
        }

        data[jhi] = data[jlo]; // promotes the smaller of the branches
        jhi = jlo;             // move down the heap
        jlo = 2*jlo + 1;       // calculates position of left branch
    }

    data[jhi] = *v; // places v, vind at correct position in heap

}

//
//    possibly insert value, displacing larger values
//    it is asssumed the data are already a heap other than
//    the first item
//

static inline void match_heap_insert(match_vector* self, const Match* match)
{
    if (match->cosdist > self->data[0].cosdist) {
        self->data[0] = *match;
        if (vector_size(self) > 1) {
            match_heap_fix(self);
        }
    }
}

//
// uses more memory than we need, need to make a simpler point
// struct
point_vector* make_points(PyObject* raObj, PyObject* decObj, int* status)
{
    size_t i=0, n=0;
    point_vector* points=NULL;
    Point pt={0};
    double *raptr=NULL, *decptr=NULL;

    n = PyArray_SIZE(raObj);
    points = point_vector_new();
    vector_resize(points, n);

    for (i=0; i<n; i++) {

        raptr=PyArray_GETPTR1(raObj, i);
        decptr=PyArray_GETPTR1(decObj, i);

        *status=hpix_eq2xyz(*raptr, *decptr, &pt.x, &pt.y, &pt.z);
        if (! (*status) ) {
            goto _make_points_bail;
        }

        vector_set(points, i, pt); 

    }

    *status=1;

_make_points_bail:

    if (! (*status) ) {
        vector_free(points);
    }
    return points;
}

//
//   fill the entry with xyz and disc_pixels, as well
//   as resetting the matches
//

static int fill_catalog_entry(CatalogEntry* entry,
                              struct healpix* hpix,
                              double ra,
                              double dec,
                              double radius)
{
    int status=0;
    CatPoint* cpt=&entry->point;

    status=hpix_eq2xyz(ra, dec, &cpt->x, &cpt->y, &cpt->z);
    if (!status) {
        // we expect the error to be set already
        goto _fill_catalog_entry_bail;
    }

    cpt->radius = radius*D2R;
    cpt->cos_radius = cos( cpt->radius );

    hpix_disc_intersect(hpix, cpt->x, cpt->y, cpt->z, cpt->radius, entry->disc_pixels);
    vector_resize(entry->matches, 0);

    status=1;
_fill_catalog_entry_bail:

    return status;
}

//
// create the tree.  Put all positions into a binary
// tree based on healpix id; indices for the objects are held
// in the tree
//

static struct tree_node* create_hpix_tree(struct healpix* hpix,
                                          PyObject* raObj,
                                          PyObject* decObj,
                                          int *status)
{

    struct tree_node* tree=NULL;
    int64_t hpixid=0;
    int64_t half_npix=0;
    size_t i=0, n=0;
    double *raptr=NULL, *decptr=NULL;
    
    // this will produce a more balanced tree across the whole sky
    half_npix = hpix->npix/2;

    n = PyArray_SIZE(raObj);
    for (i=0; i<n ; i++) {

        raptr=PyArray_GETPTR1(raObj, i);
        decptr=PyArray_GETPTR1(decObj, i);

        hpixid = hpix_eq2pix(hpix, *raptr, *decptr, status);
        if ( !(*status) ) {
            PyErr_SetString(PyExc_ValueError, "Could not get hpix id, band ra,dec\n");
            goto _create_hpix_tree_bail;
        }

        tree_insert(&tree, hpixid-half_npix, i);

    }

_create_hpix_tree_bail:

    if ( !(*status) ) {
        tree = tree_delete(tree);
    }

    return tree;
}


//
// initialize the python catalog object
//

static int
PySMatchCat_init(struct PySMatchCat* self, PyObject *args, PyObject *kwds)
{
    PY_LONG_LONG nside=0;
    int err=0;

    if (!PyArg_ParseTuple(args, (char*)"L", &nside)) {
        return -1;
    }

    self->hpix = hpix_new((int64_t)nside);
    if (self->hpix==NULL) {
        err=1;
        goto _catalog_init_cleanup;
    }

_catalog_init_cleanup:
    if (err != 0) {
        self->hpix = hpix_delete(self->hpix);
        return -1;
    }
    return 0;
}


//
// deallocate the python object
//

static void
PySMatchCat_dealloc(struct PySMatchCat* self)
{

    self->hpix = hpix_delete(self->hpix);

#if PY_MAJOR_VERSION >= 3
    Py_TYPE(self)->tp_free((PyObject*)self);
#else
    self->ob_type->tp_free((PyObject*)self);
#endif

}

//
// a repr for the python object
//

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

//
// getters
//

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

/*

   Match the input catalog entry to the second set of points, the
   hpix ids of which have been put into a tree. 
   
   If no restriction is set on maximum number of matches, then matches are
   simply appended to the match vector in the catalog entry.

   If the number of matches is restricted (maxmatch > 0) then matches are
   appended up to the max allowed, then the match vector is converted to a heap
   and only matches closer than the farthest current match are added.

*/

static void domatch1(struct PySMatchCat* self, 
                     CatalogEntry* entry,
                     size_t cat_ind,
                     struct tree_node* tree,
                     point_vector* points)
{

    CatPoint *cpt=NULL;
    Point *pt=NULL;

    int64_t hpixid=0;
    struct tree_node* node=NULL;
    int64_t half_npix=0;

    size_t i=0, j=0, input_ind=0;
    double cos_angle=0;

    int64_t maxmatch = self->maxmatch;

    Match match={0};
    match_vector* matches=NULL;

    half_npix = self->hpix->npix/2;

    matches = entry->matches;
    cpt = &entry->point;

    // loop over pixels that intersected a disc around
    // this object

    for (i=0; i < vector_size(entry->disc_pixels); i++) {

        // get the tree node corresponding to this pixel
        hpixid = vector_get(entry->disc_pixels, i);
        node = tree_find(tree, hpixid-half_npix);

        if (node != NULL) {
            for (j=0; j < vector_size(node->indices); j++) {

                input_ind = vector_get(node->indices, j);

                if (self->matching_self && input_ind==cat_ind) {
                    continue;
                }

                pt = &points->data[input_ind];

                cos_angle = pt->x*cpt->x + pt->y*cpt->y + pt->z*cpt->z;

                if (cos_angle > cpt->cos_radius) {
                    match.cat_ind=cat_ind;
                    match.input_ind=(int64_t)input_ind;
                    match.cosdist=cos_angle;


                    if (maxmatch <= 0 || (int64_t)vector_size(matches) < maxmatch) {
                        // we increment for new matches
                        self->nmatches += 1;

                        // just keep adding entries
                        vector_push(matches, match);

                        // if we are now at capacity, heapify it unless maxmatch
                        // is size one, in which case it is already a heap
                        if (maxmatch > 1 && (int64_t)vector_size(matches)==maxmatch) {
                            match_build_heap(matches);
                        }
                    } else {
                        // note the number of matches is *not* incremented

                        // add only if closer than the farthest match
                        match_heap_insert(matches, &match);
                    }

                } // within distance

            } // loop over indices in node
        } // id found in tree
    } // loop over disc pixel ids

}

static int domatch1_2file_all(struct PySMatchCat* self, 
                              CatalogEntry* entry,
                              size_t cat_ind,
                              struct tree_node* tree, // second cat hpix tree
                              point_vector* points, // second cat points
                              FILE* fobj)
{
    int status=0, nret=0;

    CatPoint *cpt=NULL;
    Point *pt=NULL;

    int64_t hpixid=0;
    struct tree_node* node=NULL;
    int64_t half_npix=0;

    size_t i=0, j=0, input_ind=0;
    double cos_angle=0;

    half_npix = self->hpix->npix/2;

    cpt = &entry->point;

    // loop over pixels that intersected a disc around
    // this object

    for (i=0; i < vector_size(entry->disc_pixels); i++) {

        // get the tree node corresponding to this pixel
        hpixid = vector_get(entry->disc_pixels, i);
        node = tree_find(tree, hpixid-half_npix);

        if (node != NULL) {
            for (j=0; j < vector_size(node->indices); j++) {

                input_ind = vector_get(node->indices, j);
                if (self->matching_self && input_ind==cat_ind) {
                    continue;
                }

                pt = &points->data[input_ind];

                cos_angle = pt->x*cpt->x + pt->y*cpt->y + pt->z*cpt->z;

                if (cos_angle > cpt->cos_radius) {

                    self->nmatches += 1;
                    nret = fprintf(fobj, "%ld %ld %.16g\n", cat_ind, input_ind, cos_angle);
                    if (nret == 0) {
                        status = 0;
                        goto _domatch1_2file_all_bail;
                    }

                } // within distance

            } // loop over indices in node
        } // id found in tree
    } // loop over disc pixel ids

    status=1;

_domatch1_2file_all_bail:
    return status;
}


//
// do the match for each entered point
// all matches are saved in memory
// run match_prep() before calling this function
//

static int domatch(struct PySMatchCat* self,
                   PyObject* craObj,
                   PyObject* cdecObj,
                   PyObject* cradiusObj,
                   PyObject* raObj,
                   PyObject* decObj,
                   PyObject* matchesObj) {
    int status=0;
    size_t i=0, nrad=0, ncat=0;

    struct tree_node* tree=NULL;

    point_vector* points=NULL;
    double cra=0, cdec=0, crad=0;

    CatalogEntry *entry=NULL;

    np_match_vector nv = {0};

    nv.data = matchesObj;
    nv.capacity = PyArray_SIZE(matchesObj);
    nv.size = 0;

    ncat=(size_t)PyArray_SIZE(craObj);
    nrad=(size_t)PyArray_SIZE(cradiusObj);

    tree = create_hpix_tree(self->hpix, raObj, decObj, &status);
    if (!status) {
        goto _domatch_bail;
    }

    points = make_points(raObj, decObj, &status);
    if (!status) {
        goto _domatch_bail;
    }

    entry = cat_entry_new();

    self->nmatches=0;

    for (i=0; i < ncat ; i++) {
        cra = *(double *)PyArray_GETPTR1(craObj,i);
        cdec = *(double *)PyArray_GETPTR1(cdecObj,i);
        if (nrad==1) {
            crad = *(double *)PyArray_GETPTR1(cradiusObj,0);
        } else {
            crad = *(double *)PyArray_GETPTR1(cradiusObj,i);
        }

        fill_catalog_entry(entry, self->hpix, cra, cdec, crad);

        domatch1(self, entry, i, tree, points);

        np_match_vector_push_many(&nv, entry->matches);

    }

    // make sure final array has exactly the desired size
    if (nv.capacity > nv.size) {
        np_match_vector_realloc(&nv, nv.size);
    }

_domatch_bail:

    tree = tree_delete(tree);
    vector_free(points);
    cat_entry_free(entry);

    return status;
}

//
// match the "catalog" to another data set
//

static PyObject* PySMatchCat_match(struct PySMatchCat* self, PyObject *args)
{
    int status=0;
    PyObject* craObj=NULL;
    PyObject* cdecObj=NULL;
    PyObject* cradiusObj=NULL;
    PyObject* raObj=NULL;
    PyObject* decObj=NULL;

    PyObject* matchesObj=NULL;

    if (!PyArg_ParseTuple(args, (char*)"LiOOOOOO",
                          &self->maxmatch,
                          &self->matching_self,
                          &craObj,
                          &cdecObj,
                          &cradiusObj,
                          &raObj,
                          &decObj,
                          &matchesObj)) {
        return NULL;
    }

    status=domatch(self,
                   craObj,
                   cdecObj,
                   cradiusObj,
                   raObj,
                   decObj,
                   matchesObj);

    if (!status) {
        return NULL;
    } else {
        Py_RETURN_NONE;
    }
}


//
// write from a match vector to a file
//

static int write_matches(match_vector* matches, FILE *fobj)
{
    size_t i=0;
    Match* match=NULL;
    int status=0, nret=0;

    for (i=0; i<vector_size(matches); i++) {

        match = &vector_get(matches, i);
        nret = fprintf(fobj, "%ld %ld %.16g\n",
                       match->cat_ind, match->input_ind, match->cosdist);
        if (nret == 0) {
            status=0;
            goto _write_matches_bail;
        }

    }

    status=1;

_write_matches_bail:
    return status;
}

//
// do matching while writing to a file
//

static int domatch2file(struct PySMatchCat* self,
                        PyObject* craObj,
                        PyObject* cdecObj,
                        PyObject* cradiusObj,
                        PyObject* raObj,
                        PyObject* decObj,
                        const char* filename) {
    int status=0;
    size_t i=0, nrad=0, ncat=0;
    struct tree_node* tree=NULL;
    point_vector* points=NULL;

    double cra=0, cdec=0, crad=0;
    CatalogEntry *entry=NULL;

    FILE* fobj=NULL;

    fobj=fopen(filename, "w");
    if (fobj == NULL) {
        PyErr_Format(PyExc_IOError, "Could not open file for writing: '%s'", filename);
        goto _domatch2file_bail;
    }

    ncat=(size_t)PyArray_SIZE(craObj);
    nrad=(size_t)PyArray_SIZE(cradiusObj);

    tree = create_hpix_tree(self->hpix, raObj, decObj, &status);
    if (!status) {
        goto _domatch2file_bail;
    }

    points = make_points(raObj, decObj, &status);
    if (!status) {
        goto _domatch2file_bail;
    }

    entry = cat_entry_new();

    self->nmatches=0;
    status=1;

    for (i=0; i< ncat ; i++) {

        cra = *(double *)PyArray_GETPTR1(craObj,i);
        cdec = *(double *)PyArray_GETPTR1(cdecObj,i);
        if (nrad==1) {
            crad = *(double *)PyArray_GETPTR1(cradiusObj,0);
        } else {
            crad = *(double *)PyArray_GETPTR1(cradiusObj,i);
        }

        fill_catalog_entry(entry, self->hpix, cra, cdec, crad);

        if (self->maxmatch > 0) {

            domatch1(self, entry, i, tree, points);

            if (vector_size(entry->matches) > 0) {
                status = write_matches(entry->matches, fobj);
            }

        } else {
            status = domatch1_2file_all(self, entry, i, tree, points, fobj);
        }


        if (!status) {
            goto _domatch2file_bail;
        }

    }

_domatch2file_bail:

    if (fobj) {
        fclose(fobj);
    }

    tree = tree_delete(tree);
    vector_free(points);
    cat_entry_free(entry);

    return status;
}


/*

   do the matching and write to the indicated file

*/
static PyObject* PySMatchCat_match2file(struct PySMatchCat* self, PyObject *args)
{
    int status=0;
    PyObject* craObj=NULL;
    PyObject* cdecObj=NULL;
    PyObject* cradiusObj=NULL;
    PyObject* raObj=NULL;
    PyObject* decObj=NULL;
    const char *filename=NULL;

    if (!PyArg_ParseTuple(args, (char*)"LiOOOOOs",
                          &self->maxmatch,
                          &self->matching_self,
                          &craObj,
                          &cdecObj,
                          &cradiusObj,
                          &raObj,
                          &decObj,
                          &filename)) {
        return NULL;
    }

    status=domatch2file(self,
                        craObj,
                        cdecObj,
                        cradiusObj,
                        raObj,
                        decObj,
                        filename);

    if (!status) {
        return NULL;
    } else {
        Py_RETURN_NONE;
    }

}

//
// count lines in a file.  Used to read matches from a file
//

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
    nmatches = PyArray_SIZE(matchesObj);

    if (nmatches <= 0) {
        // nothing to do
        Py_RETURN_NONE;
    }

    fobj=fopen(filename, "r");
    if (fobj == NULL) {
        PyErr_Format(PyExc_IOError, "Could not open file: '%s'\n", filename);
        return NULL;
    }

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


static PyMethodDef PySMatchCat_methods[] = {
    {"get_nmatches",           (PyCFunction)PySMatchCat_nmatches,          METH_VARARGS,  "Get the number of matches."},
    {"get_hpix_nside",              (PyCFunction)PySMatchCat_hpix_nside,          METH_VARARGS,  "Get the nside for healpix."},
    {"get_hpix_area",              (PyCFunction)PySMatchCat_hpix_area,          METH_VARARGS,  "Get the nside for healpix."},
    {"match",              (PyCFunction)PySMatchCat_match,          METH_VARARGS,  "Match the catalog to the input ra,dec arrays."},
    {"match2file",              (PyCFunction)PySMatchCat_match2file,          METH_VARARGS,  "Match the catalog to the input ra,dec arrays and write results to a file."},
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
