from ctypes import cdll, cast
from ctypes import c_void_p, c_double, c_int, c_long, c_size_t, c_char_p
import os
import pickle


_path = os.path.dirname(os.path.abspath(__file__))
_lib = cdll.LoadLibrary(os.path.join(_path, 'libparham.so'))

_lib.compute_hamming_matrix.argtypes = (
        c_size_t, c_size_t,
        c_char_p, c_void_p, c_void_p,
        c_char_p, c_char_p)

def compute_hamming_matrix(seqss, pos_freq, pos_ref, out=None, net=None):
    n = len(seqss)
    m = len(seqss[0])
    seq_str = ''.join([str(s) for s in seqss])
    pos_freq_arr = [pos_freq[i] for i in range(len(pos_freq))]
    pos_freq_arr = (c_double * len(pos_freq_arr))(*pos_freq_arr)
    pos_ref_arr = [pos_ref[i] for i in range(len(pos_ref))]
    pos_ref_arr = (c_int * len(pos_ref_arr))(*pos_ref_arr)
    if out is not None:
        out = out.encode('ascii')
    if net is not None:
        net = net.encode('ascii')
    _lib.compute_hamming_matrix(n, m,
            seq_str.encode('ascii'), pos_freq_arr, pos_ref_arr,
            out, net)
