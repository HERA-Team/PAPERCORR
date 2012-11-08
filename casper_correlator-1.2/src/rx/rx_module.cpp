#include "include/rx_module.h"

/*____                ____            _        _   
 / ___|___  _ __ _ __|  _ \ __ _  ___| | _____| |_ 
| |   / _ \| '__| '__| |_) / _` |/ __| |/ / _ \ __|
| |__| (_) | |  | |  |  __/ (_| | (__|   <  __/ |_ 
 \____\___/|_|  |_|  |_|   \__,_|\___|_|\_\___|\__|*/

// Deallocate memory when Python object is deleted
static void CorrPktObject_dealloc(CorrPktObject* self) {
    //PyObject_GC_UnTrack(self);
    self->ob_type->tp_free((PyObject*)self);
}

// Allocate memory for Python object 
static PyObject *CorrPktObject_new(PyTypeObject *type,
        PyObject *args, PyObject *kwds) {
    CorrPktObject *self;
    self = (CorrPktObject *) type->tp_alloc(type, 0);
    return (PyObject *) self;
}

// Initialize object (__init__)
static int CorrPktObject_init(CorrPktObject *self) {
    init_packet(&(self->pkt));
    return 0;
}

// Return HEADER_SIZE
static PyObject * CorrPktObject_header_size(CorrPktObject *self) {
    return PyInt_FromLong((long) HEADER_SIZE(self->pkt));
}

// Thin wrapper over data_size()
static PyObject * CorrPktObject_size(CorrPktObject *self) {
    return PyInt_FromLong((long) HEADER_SIZE(self->pkt) + self->pkt.pktinfo.len);
}

static PyObject * CorrPktObject_get_data(CorrPktObject *self) {
    return PyString_FromStringAndSize(self->pkt.data, self->pkt.pktinfo.len);
}

static PyObject * CorrPktObject_set_data(CorrPktObject *self, PyObject *args) {
    char *data;
    int size, i;
    PyObject *data_str;
    if (!PyArg_ParseTuple(args, "O", &data_str)) return NULL;
    CHK_STRING(data_str);
    size = PyString_Size(data_str);
    if (size > MAX_PAYLOAD_SIZE) {
        PyErr_Format(PyExc_ValueError, 
            "length of string (%d) exceeds maximum payload size of %d", 
            size, MAX_PAYLOAD_SIZE);
        return NULL;
    }
    data = PyString_AsString(data_str);
    for (i=0; i < size; i++) {
        self->pkt.data[i] = data[i];
    }
    self->pkt.pktinfo.len = size;
    Py_INCREF(Py_None);
    return Py_None;
}

/*
// Get packet information as a dict (so we can check it)
static PyObject * CorrPktObject_to_dict(CorrPktObject *self) {
    return Py_BuildValue("{s:i,s:{s:i,s:i,s:i,s:i,s:i,s:i},s:{s:i,s:i},s:i,s:i,s:i,s:i,s:s#}",
        "n_options", self->pkt.n_options,
        "insttype",
            "inst_id", self->pkt.insttype.inst_id,
            "is_complex", self->pkt.insttype.is_complex,
            "msb_first", self->pkt.insttype.msb_first,
            "data_type", self->pkt.insttype.data_type,
            "bits_per_val", self->pkt.insttype.bits_per_val,
            "vals_per_item", self->pkt.insttype.vals_per_item,
        "instsrc",
            "instance_id", self->pkt.instsrc.instance_id,
            "engine_id", self->pkt.instsrc.engine_id,
        "timestamp", self->pkt.timestamp,
        "datalen", self->pkt.datalen,
        "dataoff", self->pkt.dataoff,
        "preverr", self->pkt.preverr,
        "data", self->pkt.data, self->pkt.datalen);
}
*/

// Thin wrapper over unpack_header()
static PyObject * CorrPktObject_unpack_header(CorrPktObject *self,
        PyObject *args) {
    char *data;
    PyObject *data_str;
    if (!PyArg_ParseTuple(args, "O", &data_str)) return NULL;
    CHK_STRING(data_str);
    data = PyString_AsString(data_str);
    try {
        unpack_header(&(self->pkt), data);
    } catch (PacketError &e) {
        PyErr_Format(PyExc_ValueError, "%s", e.get_message());
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// Thin wrapper over unpack_data()
static PyObject * CorrPktObject_unpack_data(CorrPktObject *self,
        PyObject *args) {
    char *data;
    PyObject *data_str;
    if (!PyArg_ParseTuple(args, "O", &data_str)) return NULL;
    CHK_STRING(data_str);
    data = PyString_AsString(data_str);
    try {
        unpack_data(&(self->pkt), data);
    } catch (PacketError &e) {
        PyErr_Format(PyExc_ValueError, "%s", e.get_message());
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// Thin wrapper over unpack()
static PyObject * CorrPktObject_unpack(CorrPktObject *self,
        PyObject *args) {
    char *data;
    PyObject *data_str;
    if (!PyArg_ParseTuple(args, "O", &data_str)) return NULL;
    CHK_STRING(data_str);
    data = PyString_AsString(data_str);
    try {
        unpack(&(self->pkt), data);
    } catch (PacketError &e) {
        PyErr_Format(PyExc_ValueError, "%s", e.get_message());
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// Thin wrapper over pack_header()
static PyObject * CorrPktObject_pack_header(CorrPktObject *self) {
    char data[MAX_HEADER_SIZE];
    pack_header(self->pkt, data);
    return PyString_FromStringAndSize(data, HEADER_SIZE(self->pkt));
}

// Thin wrapper over pack_data()
static PyObject * CorrPktObject_pack_data(CorrPktObject *self) {
    char data[MAX_PAYLOAD_SIZE];
    pack_data(self->pkt, data);
    return PyString_FromStringAndSize(data, (long) self->pkt.pktinfo.len);
}

// Thin wrapper over pack()
static PyObject * CorrPktObject_pack(CorrPktObject *self) {
    char data[MAX_HEADER_SIZE+MAX_PAYLOAD_SIZE];
    pack(self->pkt, data);
    return PyString_FromStringAndSize(data,
        self->pkt.pktinfo.len+HEADER_SIZE(self->pkt));
}

// Bind methods to object
static PyMethodDef CorrPktObject_methods[] = {
    {"header_size", (PyCFunction)CorrPktObject_header_size, METH_NOARGS,
        "header_size()\nReturn the size (in bytes) of the header."},
    {"size", (PyCFunction)CorrPktObject_size, METH_NOARGS,
        "size()\nReturn the size (in bytes) of the packet."},
    {"get_data", (PyCFunction)CorrPktObject_get_data, METH_NOARGS,
        "get_data()\nReturn packet data as a binary string. This is different from pack_data() in that it does not change endianess."},
    {"set_data", (PyCFunction)CorrPktObject_set_data, METH_VARARGS,
        "set_data()\nSet packet data to a binary string. This is different from unpack_data() in that it does not change endianess."},
    {"unpack_header", (PyCFunction)CorrPktObject_unpack_header, METH_VARARGS,
        "unpack_header(data)\nSet header from binary string. Raise ValueError if data doesn't match packet format."},
    {"unpack_data", (PyCFunction)CorrPktObject_unpack_data, METH_VARARGS,
        "unpack_data(data)\nSet data from binary string."},
    {"unpack", (PyCFunction)CorrPktObject_unpack, METH_VARARGS,
        "unpack(data)\nSet packet from binary string. Raise ValueError if data doesn't match packet format."},
    {"pack_header", (PyCFunction)CorrPktObject_pack_header, METH_NOARGS,
        "pack_header()\nReturn header as a binary string."},
    {"pack_data", (PyCFunction)CorrPktObject_pack_data, METH_NOARGS,
        "pack_data()\nReturn data as a binary string."},
    {"pack", (PyCFunction)CorrPktObject_pack, METH_NOARGS,
        "pack()\nReturn packet as a binary string."},
    {NULL}  // Sentinel
};

static PyMemberDef CorrPktObject_members[] = {
    {"n_options", T_USHORT, offsetof(CorrPktObject, pkt) + 
        offsetof(CorrPacket, n_options), 0, "n_options"},
    {"instrument_id", T_USHORT, offsetof(CorrPktObject, pkt) + 
        offsetof(CorrPacket, instids) + offsetof(InstIds, instrument_id), 0, 
        "instrument_id"},
    {"instance_id", T_USHORT, offsetof(CorrPktObject, pkt) + 
        offsetof(CorrPacket, instids) + offsetof(InstIds, instance_id), 0, 
        "instance_id"},
    {"packet_len", T_UINT, offsetof(CorrPktObject, pkt) + 
	offsetof(CorrPacket, pktinfo) + offsetof(PktInfo, len), 0,
	"packet_len"},
    {"packet_count", T_UINT, offsetof(CorrPktObject, pkt) +
        offsetof(CorrPacket, pktinfo) + offsetof(PktInfo, count), 0,
        "packet_count"},
    {"engine_id", T_USHORT, offsetof(CorrPktObject, pkt) + 
        offsetof(CorrPacket, instids) + offsetof(InstIds, engine_id), 0, 
        "engine_id"},
    {"timestamp", T_ULONG, offsetof(CorrPktObject, pkt) + 
        offsetof(CorrPacket, timestamp), 0, "timestamp"},
    {"heap_len", T_ULONG, offsetof(CorrPktObject, pkt) + 
        offsetof(CorrPacket, heap_len), READONLY, "heap_len"},
    {"heap_off", T_ULONG, offsetof(CorrPktObject, pkt) + 
        offsetof(CorrPacket, heap_off), 0, "heap_off"},
    {"currerr", T_ULONG, offsetof(CorrPktObject, pkt) + 
        offsetof(CorrPacket, currerr), 0, "currerr"},
    {NULL}  /* Sentinel */
};

PyTypeObject CorrPktType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "CorrPacket", /*tp_name*/
    sizeof(CorrPktObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)CorrPktObject_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "This class provides the basic interfaces correlator packets.  CorrPacket()",       /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    CorrPktObject_methods,     /* tp_methods */
    CorrPktObject_members,     /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)CorrPktObject_init,      /* tp_init */
    0,                         /* tp_alloc */
    CorrPktObject_new,       /* tp_new */
};


/*____      _ _       _       ____         __  __           
 / ___|___ | | | __ _| |_ ___| __ ) _   _ / _|/ _| ___ _ __ 
| |   / _ \| | |/ _` | __/ _ \  _ \| | | | |_| |_ / _ \ '__|
| |__| (_) | | | (_| | ||  __/ |_) | |_| |  _|  _|  __/ |   
 \____\___/|_|_|\__,_|\__\___|____/ \__,_|_| |_|  \___|_|   */

// Deallocate memory when Python object is deleted
static void ColBufObject_dealloc(ColBufObject* self) {
    //PyObject_GC_UnTrack(self);
    free_collate_buffer(self->cb);
    if (self->pycallback != NULL) Py_DECREF(self->pycallback);
    self->ob_type->tp_free((PyObject*)self);
}

// Allocate memory for Python object 
static PyObject *ColBufObject_new(PyTypeObject *type,
        PyObject *args, PyObject *kwds) {
    ColBufObject *self;
    self = (ColBufObject *) type->tp_alloc(type, 0);
    return (PyObject *) self;
}

// Initialize object (__init__)
static int ColBufObject_init(ColBufObject *self,
        PyObject *args, PyObject *kwds) {
    int nant=16, nants_per_feng=4, nchan=2048, xeng_chan_mode=0, npol=4, nwin=1, sdisp=0, acc_len=1;
    char *sdisp_destination_ip = "127.0.0.1";
    static char *kwlist[] = {
      "nant",
      "nants_per_feng",
      "nchan",
      "xeng_chan_mode",
      "npol",
      "nwin",
      "sdisp",
      "sdisp_destination_ip",
      "acc_len",
      NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds,"|iiiiiiisi", kwlist, \
            &nant,
            &nants_per_feng,
            &nchan,
            &xeng_chan_mode,
            &npol,
            &nwin,
            &sdisp,
            &sdisp_destination_ip,
            &acc_len))
        return -1;
    try {
        init_collate_buffer(&(self->cb), nant, nants_per_feng, nchan, xeng_chan_mode, npol, nwin, sdisp, sdisp_destination_ip, acc_len);
    } catch (PacketError &e) {
        PyErr_Format(PyExc_MemoryError,
        "Couldn't allocate CollateBuffer for nant=%d,nchan=%d,npol=%d,nwin=%d",
            nant, nchan, npol, nwin);
        return -1;
    }
    self->pycallback = NULL;
    return 0;
}

// Get xengine baseline order from CollateBuffer
static PyObject * ColBufObject_xeng_bl_order(ColBufObject *self) {
    int i;
    PyObject *rv;
    rv = PyTuple_New(self->cb.nbl);
    for (i=0; i < self->cb.nbl; i++) {
        PyTuple_SET_ITEM(rv, i, Py_BuildValue("(ii)",
            self->cb.xeng_ai_order[i], self->cb.xeng_aj_order[i]));
    }
    return rv;
}

// A thin wrapper over collate_packet
static PyObject * ColBufObject_collate_packet(ColBufObject *self,
        PyObject *args) {
    int rv;
    CorrPktObject *pkt;
    if (!PyArg_ParseTuple(args, "O!", &CorrPktType, &pkt))
        return NULL;
    try {
        rv = collate_packet(&(self->cb), pkt->pkt);
    } catch (PacketError &e) {
        PyErr_Format(PyExc_ValueError, "%s", e.get_message());
        return NULL;
    }
    // If rv == 1, then we have a Python callback that raised an exception
    // so the error message should already be set.
    if (rv) return NULL;
    // Otherwise, we're good to go
    Py_INCREF(Py_None);
    return Py_None;
}

// Routine for inspecting a CollateBuffer's buffer
static PyObject * ColBufObject_get_buffer(ColBufObject *self) {
    PyArrayObject *dbuf, *fbuf;
    npy_intp npy_dims[1];
    npy_dims[0] = self->cb.buflen;
    dbuf = (PyArrayObject *) PyArray_SimpleNewFromData(1, npy_dims, 
            NPY_INT, self->cb.buf);
    npy_dims[0] = self->cb.buflen/2;
    fbuf = (PyArrayObject *) PyArray_SimpleNewFromData(1, npy_dims, 
            NPY_INT, self->cb.flagbuf);
    CHK_NULL(dbuf);
    CHK_NULL(fbuf);
    return Py_BuildValue("(OO)", PyArray_Return(dbuf), PyArray_Return(fbuf));
}

int wrap_cb_pycallback(int i, int j, int pol, int64_t t,
//int wrap_cb_pycallback(int i, int j, int pol, double t,
        float *data, int *flags, int nchan, void *userdata) {
    ColBufObject *cbo;
    PyArrayObject *ndata, *nflags;
    npy_intp npy_dims[1];
    PyObject *arglist, *rv;
    npy_dims[0] = nchan;
    PyGILState_STATE gstate;
    // Acquire Python Global Interpeter Lock
    //printf("wrap_cb_pycallback\n");
    gstate = PyGILState_Ensure();
    cbo = (ColBufObject *) userdata;  // Recast userdata as reference to a cb
    ndata = (PyArrayObject *) PyArray_SimpleNewFromData(1, npy_dims, 
            NPY_CFLOAT, data);
    nflags = (PyArrayObject *) PyArray_SimpleNewFromData(1, npy_dims, 
            NPY_INT, flags);
    // data/flags now embedded in ndata/nflags, but will NOT be freed when
    // these objects are deleted.  Underlying data may/will change and/or be
    // freed after the callback completes.  Callbacks wishing to hold onto data
    // must duplicate it themselves (but that is not a typical use case).
    arglist = Py_BuildValue("(iiilOO)", i, j, pol, t, 
            PyArray_Return(ndata), PyArray_Return(nflags));
    Py_DECREF(ndata); Py_DECREF(nflags);
    rv = PyEval_CallObject(cbo->pycallback, arglist);
    Py_DECREF(arglist);
    if (rv == NULL) {
        PyGILState_Release(gstate);
        return 1;
    }
    Py_DECREF(rv);
    // Release Python Global Interpeter Lock
    PyGILState_Release(gstate);
    return 0;
}

// Routine for setting a python callback for CollateBuffer data output
static PyObject * ColBufObject_set_callback(ColBufObject *self, PyObject *args){
    PyObject *cbk;
    if (!PyArg_ParseTuple(args, "O", &cbk))
        return NULL;
    if (!PyCallable_Check(cbk)) {
        PyErr_SetString(PyExc_TypeError, "parameter must be callable");
        return NULL;
    }
    Py_INCREF(cbk);
    if (self->pycallback != NULL) Py_DECREF(self->pycallback);
    (self->cb).userdata = (void *)self;
    self->pycallback = cbk;
    set_cb_callback(&self->cb, &wrap_cb_pycallback);
    Py_INCREF(Py_None);
    return Py_None;
}
    
// Routine for removing a python callback for CollateBuffer data output
static PyObject * ColBufObject_unset_callback(ColBufObject *self) {
    set_cb_callback(&self->cb, &default_callback);
    if (self->pycallback != NULL) Py_DECREF(self->pycallback);
    self->pycallback = NULL;
    (self->cb).userdata = NULL;
    Py_INCREF(Py_None);
    return Py_None;
}
    
// Bind methods to object
static PyMethodDef ColBufObject_methods[] = {
    {"xeng_bl_order", (PyCFunction)ColBufObject_xeng_bl_order, METH_NOARGS,
        "xeng_bl_order()\nReturn baseline order output by X engines."},
    {"collate_packet", (PyCFunction)ColBufObject_collate_packet, METH_VARARGS,
        "collate_packet(pkt)\nFile X engine packet into buffer."},
    {"get_buffer", (PyCFunction)ColBufObject_get_buffer, METH_NOARGS,
        "get_buffer()\nReturn entire packet buffer as a numpy array."},
    {"set_callback", (PyCFunction)ColBufObject_set_callback, METH_VARARGS,
        "set_callback()\nSet a callback for handling output data.  Callback should accept the arguments (i,j,pol,time,data,flags), where i,j are antenna indices, pol is a polarization index, time is a correlator integration count (NOT a real-world time), data is a complex numpy array representing the spectrum from a single integration, and flags is a boolean numpy array of data flags with the same shape as data."},
    {"unset_callback", (PyCFunction)ColBufObject_unset_callback, METH_NOARGS,
        "unset_callback()\nClear a callback (reseting to default) for handling output data."},
    {NULL}  // Sentinel
};

PyTypeObject ColBufType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "CollateBuffer", /*tp_name*/
    sizeof(ColBufObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)ColBufObject_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "This class collates correlator packets to produce per-baseline spectra.  CollateBuffer(n_ants=16,nchan=2048,npol=4,nwin=1)",       /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    ColBufObject_methods,             /* tp_methods */
    0,                     /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)ColBufObject_init,      /* tp_init */
    0,                         /* tp_alloc */
    ColBufObject_new,       /* tp_new */
};


/*___         __  __           ____             _        _   
| __ ) _   _ / _|/ _| ___ _ __/ ___|  ___   ___| | _____| |_ 
|  _ \| | | | |_| |_ / _ \ '__\___ \ / _ \ / __| |/ / _ \ __|
| |_) | |_| |  _|  _|  __/ |   ___) | (_) | (__|   <  __/ |_ 
|____/ \__,_|_| |_|  \___|_|  |____/ \___/ \___|_|\_\___|\__| */

// Deallocate memory when Python object is deleted
static void BsockObject_dealloc(BsockObject* self) {
    free_buffer_socket(&self->bs);
    if (self->pycallback) Py_DECREF(self->pycallback);
    self->ob_type->tp_free((PyObject*)self);
}

// Allocate memory for Python object 
static PyObject *BsockObject_new(PyTypeObject *type,
        PyObject *args, PyObject *kwds) {
    BsockObject *self;
    self = (BsockObject *) type->tp_alloc(type, 0);
    return (PyObject *) self;
}

// Initialize object (__init__)
static int BsockObject_init(BsockObject *self, PyObject *args, PyObject *kwds) {
    int item_count=128, payload_len=8192;
    static char *kwlist[] = {"item_count","payload_len", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds,"|ii", kwlist, \
            &item_count, &payload_len))
        return -1;
    init_buffer_socket(&self->bs, item_count, payload_len);
    self->pycallback = NULL;
    return 0;
}

static PyObject * BsockObject_start(BsockObject *self, PyObject *args) {
    int port;
    if (!PyArg_ParseTuple(args, "i", &port)) return NULL;
    PyEval_InitThreads();
    start(&(self->bs), port);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject * BsockObject_stop(BsockObject *self) {
    //PyThreadState *_save;
    // Release Python Global Interpreter Lock so that a python callback end
    Py_BEGIN_ALLOW_THREADS
    stop(&(self->bs));
    // Reacquire Python Global Interpreter Lock
    Py_END_ALLOW_THREADS
    Py_INCREF(Py_None);
    return Py_None;
}

int wrap_bs_pycallback(char *data, size_t size, void *userdata) {
    BsockObject *bso;
    PyObject *arglist, *rv;
    PyGILState_STATE gstate;
    // Acquire Python Global Interpeter Lock
    gstate = PyGILState_Ensure();
    bso = (BsockObject *) userdata;  // Recast userdata as reference to a bs
    arglist = Py_BuildValue("(s#)", data, (int) size);
    rv = PyEval_CallObject(bso->pycallback, arglist);
    Py_DECREF(arglist);
    if (rv == NULL) {
        PyGILState_Release(gstate);
        return 1;
    }
    Py_DECREF(rv);
    // Release Python Global Interpeter Lock
    PyGILState_Release(gstate);
    return 0;
}

int bs_collatebuffer_callback(char *data, size_t size, void *userdata) {
    CollateBuffer *cb;
    CorrPacket pkt;
    cb = (CollateBuffer *) userdata;
    //printf("in bs_collatebuffer_callback: size=%d\n", size);
    init_packet(&pkt);
    try {
        unpack(&pkt, data);
    } catch (PacketError &e) {
        fprintf(stderr, "%s\n", e.get_message());
        return 1;
    }
    if (pkt.sync_time > 0) {
     cb->sync_time = pkt.sync_time;
     printf("Got sync_time packet\n");
     return 0;
    } 
    //printf("packet accepted: pkt->heap_off=%ld\n", pkt.heap_off);
    return collate_packet(cb, pkt);
}

// Routine for setting a python callback for BufferSocket data output
static PyObject * BsockObject_set_callback(BsockObject *self, PyObject *args){
    PyObject *cbk;
    ColBufObject *cbo;
    if (PyArg_ParseTuple(args, "O!", &ColBufType, &cbo)) {
        if (self->pycallback != NULL) Py_DECREF(self->pycallback);
        Py_INCREF(cbo);
        self->pycallback = (PyObject *)cbo;
        self->bs.userdata = (void *)&(cbo->cb);
        set_callback(&self->bs, &bs_collatebuffer_callback);
    } else {
        PyErr_Clear();
        if (!PyArg_ParseTuple(args, "O", &cbk))
            return NULL;

        if (!PyCallable_Check(cbk)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be callable");
            return NULL;
        }
        Py_INCREF(cbk);
        if (self->pycallback != NULL) Py_DECREF(self->pycallback);
        (self->bs).userdata = (void *)self;
        self->pycallback = cbk;
        set_callback(&self->bs, &wrap_bs_pycallback);
    }
    Py_INCREF(Py_None);
    return Py_None;
}
    
// Routine for removing a python callback for data output
static PyObject * BsockObject_unset_callback(BsockObject *self) {
    set_callback(&self->bs, &default_callback);
    if (self->pycallback != NULL) Py_DECREF(self->pycallback);
    self->pycallback = NULL;
    (self->bs).userdata = NULL;
    Py_INCREF(Py_None);
    return Py_None;
}
    

// Bind methods to object
static PyMethodDef BsockObject_methods[] = {
    {"start", (PyCFunction)BsockObject_start, METH_VARARGS,
        "start(port)\nBegin listening for UDP packets on the specified port."},
    {"stop", (PyCFunction)BsockObject_stop, METH_NOARGS,
        "stop()\nHalt listening for UDP packets."},
    {"set_callback", (PyCFunction)BsockObject_set_callback, METH_VARARGS,
        "set_callback(cbk)\nSet a callback function for output data from a BufferSocket.  If cbk is a CollateBuffer, a special handler is used that feeds data into the CollateBuffer without entering back into Python (for speed).  Otherwise, cbk should be a function that accepts a single argument: a binary string containing packet data."},
    {"unset_callback", (PyCFunction)BsockObject_unset_callback, METH_NOARGS,
        "unset_callback()\nReset the callback to the default."},
    {NULL}  // Sentinel
};

PyTypeObject BsockType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "BufferSocket", /*tp_name*/
    sizeof(BsockObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)BsockObject_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "A ring-buffered, multi-threaded socket interface for holding UDP packets.  BufferSocket(item_count=128, payload_len=8192)",       /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    BsockObject_methods,     /* tp_methods */
    0,                       /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)BsockObject_init,      /* tp_init */
    0,                         /* tp_alloc */
    BsockObject_new,       /* tp_new */
};

// Module methods
static PyMethodDef rx_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

// Module init
PyMODINIT_FUNC initrx(void) {
    PyObject* m;
    CorrPktType.tp_new = PyType_GenericNew;
    ColBufType.tp_new = PyType_GenericNew;
    BsockType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&CorrPktType) < 0) return;
    if (PyType_Ready(&ColBufType) < 0) return;
    if (PyType_Ready(&BsockType) < 0) return;
    m = Py_InitModule3("rx", rx_methods,
    "A module for handling low-level (high performance) packet rx from a CASPER correlator.");
    import_array();
    Py_INCREF(&BsockType);
    PyModule_AddObject(m, "BufferSocket", (PyObject *)&BsockType);
    Py_INCREF(&ColBufType);
    PyModule_AddObject(m, "CollateBuffer", (PyObject *)&ColBufType);
    Py_INCREF(&CorrPktType);
    PyModule_AddObject(m, "CorrPacket", (PyObject *)&CorrPktType);
}

