#!/usr/bin/python3
import subprocess
import argparse
import core_functions as core

"""
This is the main function that calls the core_functions module and extracts the bend angle from 
"""

def onehelix(file_basename, secstr_output,helix_type):

    filename_hr = (file_basename + ".1hr")
    filename_res = (file_basename + ".1res")
    filename_format = (file_basename + ".1format")
    filename_angle = (file_basename + ".1angle")

    print('start onehelixres')
    list_of_helices, helix_secstr, helix_resno, helix_resname = core.single_helix_parser(secstr_output)
    if not list_of_helices:
        print("%s has no one helix"%file_basename)
        pass
    else:
        proline_res_d = None

        core.filewrite_nestedlist(filename_hr, list_of_helices)
        print('finish onehelixre with %d helix ' % (len(list_of_helices)))
        subprocess.run(['pdbgetresidues', filename_hr, sel_filename, filename_res])
        dict_of_ca_in_helix = core.ca_atom_organiser(filename_hr, filename_res, filename_format)

        core.proline_and_nonproline_middle_angle(filename_format, filename_angle, helix_type,proline_res_d, verbose)

def twohelix(file_basename, secstr_output,helix_type):
    filename_hr = (file_basename + ".2hr")
    filename_res = (file_basename + ".2res")
    filename_format = (file_basename + ".2format")
    filename_angle = (file_basename + ".2angle")

    print('start twohelixres')
    list_of_helices, helix_secstr, helix_resno, helix_resname, proline_res_d = core.double_helix_parser(secstr_output)
    print("proline_res_d %s "%proline_res_d)
    if not list_of_helices:
        print("%s has no two helix"%file_basename)
        pass
    else:
        core.filewrite_nestedlist(filename_hr, list_of_helices)
        print("finish twohelixres with %d helix" % (len(list_of_helices)))
        subprocess.run(['pdbgetresidues', filename_hr, sel_filename, filename_res])
        dict_of_ca_in_helix = core.ca_atom_organiser(filename_hr, filename_res, filename_format)

        core.proline_and_nonproline_middle_angle(filename_format, filename_angle, helix_type, proline_res_d, verbose)


if __name__=="__main__":

    pdb_filename = None
    sec_filename = None
    sel_filename = None
    file_basename = None
    helix_type = None
    verbose = False



    pdb_filename, helix_type = core.optional_linux_argument()

    file_basename = pdb_filename.split('.')[0]
    sec_filename = (file_basename + ".sec")
    sel_filename = (file_basename + ".sel")
    print("Pdb file %s" % file_basename)

    secstr_process = subprocess.run(['pdbsecstr', pdb_filename], check=True, stdout=subprocess.PIPE, universal_newlines=True)
    secstr_output = secstr_process.stdout

    atomsel_process = subprocess.run(['pdbatomselect', pdb_filename,sel_filename])


    if helix_type == 1:

        onehelix(file_basename, secstr_output,1)

    if helix_type == 2:
        twohelix(file_basename, secstr_output,2)

    if helix_type == 0:
        onehelix(file_basename,secstr_output,1)
        twohelix(file_basename,secstr_output,2)


