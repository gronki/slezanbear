from ctypes import CDLL, POINTER
from ctypes import c_double, c_float, c_int, c_char
from numpy.ctypeslib import ndpointer

def CFACTORY(lib, name, args = None, restype = None):
    from ctypes import CFUNCTYPE
    flags = { 'in':1, 'out':2, 'inout':3 }
    args = args if args else []
    argtypes = tuple( a_type for a_type,a_intent,a_name in args )
    arglist = tuple( (flags[a_intent], a_name) for a_type,a_intent,a_name in args )
    proto = CFUNCTYPE(restype, *argtypes)
    f = proto((name,lib), arglist)
    f.__doc__ = '{}({})'.format(name, ', '.join([
        '{} {}'.format(a_type.__name__,a_name) \
        for a_type,a_intent,a_name in args  ]))
    return f
