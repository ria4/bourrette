#!/usr/bin/env python
from gimpfu import *
from ctypes import *
from sys import argv, platform
import sys
#sys.stderr = open('C:/temp/python-fu-output.txt','a')
#sys.stdout=sys.stderr # So that they both go to the same file

def load_library (library_name):
    if platform == "linux" or platform == "linux2":
        library_name = library_name + '.so.0'
    elif platform == "win32":
        from ctypes.util import find_library
        library_name = find_library (library_name + "-0")
    else:
        raise BaseException ("TODO")
    return CDLL (library_name)
#import sys
#print (sys.version)
gimp = load_library ('libgimp-2.0')
#print (gimp)
#gimp = load_library ('libgimp-3.0')

class GimpParamData (Union):
    _fields_ = [ ('d_float'   , c_double),
                 ('d_string'  , c_char_p),
                 ('d_color'   , c_double * 4),
                 ('d_drawable', c_int),
                 ('d_status'  , c_int),
                 ('d_int'     , c_int)]

class GimpParam (Structure):
    _fields_ = [ ('type', c_int),
                 ('data', GimpParamData) ]

GimpRunProc = CFUNCTYPE (None,
                         c_char_p,
                         c_int,
                         POINTER (GimpParam),
                         POINTER (c_int),
                         POINTER (POINTER (GimpParam)))

void_FUNC_void = CFUNCTYPE (None)

class GimpPlugInInfo (Structure):
    _fields_ = [ ('init_proc', void_FUNC_void),
                 ('quit_proc', void_FUNC_void),
                 ('query_proc', void_FUNC_void),
                 ('run_proc', GimpRunProc) ]
    def __init__ (self,
                  query_proc,
                  run_proc,
                  init_proc = void_FUNC_void (),
                  quit_proc = void_FUNC_void ()):
        Structure.__init__ (self)
        self.init_proc  = init_proc
        self.quit_proc  = quit_proc
        self.query_proc = query_proc
        self.run_proc   = run_proc

class GimpParamDef (Structure):
    _fields_ = [ ('type', c_int),
                 ('name', c_char_p),
                 ('description', c_char_p) ]

GimpParam2 = GimpParam * 2
ret_values = GimpParam2 ()

def run (proc_name, n_params, params, n_return_vals, return_vals):
    class GeglBuffer(Structure):
        pass
    return_vals[0] = ret_values
    ret_values[0].type = 21 # GIMP_PDB_STATUS

    drawable = params[0].data.d_drawable
    if gimp.gimp_item_is_group (drawable) or gimp.gimp_item_is_vectors (drawable):
        n_return_vals[0] = 2

        ret_values[0].data.d_status = PDB_CALLING_ERROR # 1

        ret_values[1].type          = PDB_STRING # 4
        ret_values[1].data.d_string = b"filter " + proc_name + b" called on a layer-group/vectors item"
        return

    n_return_vals[0] = 1
    ret_values[0].data.d_status = PDB_SUCCESS # 3

    gegl = load_library ('libgegl-0.4')
    sucess = gegl.gegl_init (None, None)
#    gegl.gegl_list_operations.restype=POINTER(c_char_p)
#    fred=gegl.gegl_list_operations()
#    i=0
#    while fred[i]!= None:
#      print(fred[i])
#      i+=1
    gimp.gimp_drawable_get_shadow_buffer.restype = POINTER (GeglBuffer)
    gimp.gimp_drawable_get_buffer.restype        = POINTER (GeglBuffer)

### List of parameter names and values
    args = specific_args(gegl,proc_name,params)
    
    x = c_int (); y = c_int (); w = c_int (); h = c_int ()
    if gimp.gimp_drawable_mask_intersect (drawable, byref (x), byref (y), byref (w), byref (h)):
        source = gimp.gimp_drawable_get_buffer (drawable)
        target = gimp.gimp_drawable_get_shadow_buffer (drawable)
#        print("gegl_render_op: ", proc_name)
        sucess = gegl.gegl_render_op (source, target, proc_name, *args)
#        print('sucess=%d' % (sucess))
        gegl.gegl_buffer_flush (target)
        gimp.gimp_drawable_merge_shadow (drawable, PushUndo = True)
        gimp.gimp_drawable_update (drawable, x, y, w, h)
        gimp.gimp_displays_flush ()

###
def specific_args (gegl,gegl_op_name,params):
    run_args = [c_void_p ()]
    # List of arg names and value pairs, terminated by a null
    if (gegl_op_name == "gegl:gegl"): 
      run_args = [b"string", c_char_p (params[1].data.d_string),
                  c_void_p ()]   
#    print(gegl_op_name, " run_args: ",run_args)
#    sys.stdout.flush()
    return run_args  
      
def specific_params (gegl_op_name):
    procedure_details = {'name': gegl_op_name}
    procedure_details['blurb'] = b"Invokes " + procedure_details['name']
    procedure_details['help'] = b"Runs "  + procedure_details['name']
    procedure_details['author'] = b"drawoC suomynonA kjp"
    procedure_details['copyright'] = b"Public Domain"
    procedure_details['date'] = b"2020"
 
    # Constants from enums.pl
    if (gegl_op_name == "gegl:gegl"): 
      params = [ GimpParamDef (PDB_DRAWABLE, b'layer', b'Input layer'),
                 GimpParamDef (PDB_STRING, b'string', b'String'),
               ]
      procedure_details['menu_label'] = b"GEGL graph"
      procedure_details['image_types'] = b"*"
    #elif (gegl_op_name == "gegl:simplex-noise"): 
    #  params = [ GimpParamDef (PDB_DRAWABLE, b'layer', b'Input layer'),
    #             GimpParamDef (PDB_INT, b'seed', b'Seed'),
    #             GimpParamDef (PDB_FLOAT, b'scale', b'Scale'),
    #             GimpParamDef (PDB_INT, b'iterations', b'Iterations'),
    #           ]
    #  procedure_details['menu_label'] = b"GEGL Simplex Noise"
    #  procedure_details['image_types'] = b"*"
#    print(params,procedure_details)
    return params,procedure_details
    
def query ():
    gegl_op_list = [
                    "gegl:gegl",
    #                "gegl:simplex-noise",
                   ]
    for gegl_op_name in gegl_op_list:
      params,procedure_details = specific_params(gegl_op_name)
      generic_query(params, procedure_details)
###
    
def generic_query (params, procedure_details):   
    GimpParamDefN = GimpParamDef * len (params)
    gimp.gimp_install_procedure (procedure_details['name'],         # name
                                 procedure_details['blurb'],        # blurb
                                 procedure_details['help'],         # help
                                 procedure_details['author'],       # author
                                 procedure_details['copyright'],    # copyright
                                 procedure_details['date'],         # date
                                 procedure_details['menu_label'],   # menu_label
                                 procedure_details['image_types'],  # image_types
                                 1,                         # GimpPDBProcType type GIMP_PLUGIN == 1
                                 len (params),              # n_params
                                 0,                         # n_return_vals
                                 GimpParamDefN (*params),   # params
                                 None)


PLUG_IN_INFO = GimpPlugInInfo (void_FUNC_void (query), GimpRunProc (run))

Argc = len (argv)

charArrayN = c_char_p * Argc
Argv = charArrayN (*[ai.encode ('utf8') for ai in argv])

#print ("Here")
gimp.gimp_main (byref (PLUG_IN_INFO), Argc, Argv)
