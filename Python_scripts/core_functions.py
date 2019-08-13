#!/usr/bin/python3
import os
import sys
import re
import numpy as np

import argparse
import subprocess

"""
Created 06.2019 by Oliver Roberts 
Changelog 
v 1.0 
v 1.1 23.07.2019


v1.1 
modified double/single_helix_res_parser 
changed input from (input_file_name) to an open file type 




functions in this module and their origin 

optional_linux_argument	    one/twohelixres.py
linux_arguments  		    one/twohelixres.py
single_helix_parser		    one/twohelixres.py
double_helix_parser		    one/twohelixres.py

ca_atom_organiser		    ca_atom_organiser
write_helix_dict		    ca_atom_organiser
first_residue_pdblines	    ca_atom_organiser
aa_chains_split			    ca_atom_organiser

filewrite_nestedlist	    ca_atom_organiser
string_nestedlist		    one/twohelixres.py
keychain_value_str		    one/twohelixres.py
keychain_value_list		    one/twohelixres.py
fileread			        one/twohelixres.py

middle_angle_linux_arguments            proline_middle_angle.py
residue_pairs_for_pdbline               proline_middle_angle.py
proline_segment_first_mid_last_finder   proline_middle_angle.py
lineinformation_extractor               proline_middle_angle.py
shell_interface                         proline_middle_angle.py
commandline_wrapper                     proline_middle_angle.py
create_pdbline_string                   proline_middle_angle.py
calculate_midpoint                      proline_middle_angle.py
calculate_angle                         proline_middle_angle.py
segment_middle_res                      proline_middle_angle.py
residue_extractor_from_ca_pdb           proline_middle_angle.py
non_proline_segment_first_mid_last_finder   proline_middle_angle.py
proline_and_nonproline_middle_angle     proline_middle_angle.py
"""

def optional_linux_argument():
    """
    Takes the first arg as a file name and the second optional argument -p to specify
    non proline (1) or proline (3) containing helices. defualt val is 1
    :return:
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'pdb_file',
        help=
        'The pdbfile that we will use for pdbsecstr '
    )
    parser.add_argument("-p", "--proline", type=int, choices=[0, 1, 2 ],
                        default=1,
                        help="1 for [H|h] helix or 2 for helix with proline HHHPHHH for example")

    args = parser.parse_args()


    if args.proline == 2:
        print("Checking for one segment helices ")

        return(args.pdb_file,args.proline)
    elif args.proline == 1:
        print("Checking for two segment helices ".format(args))

        return(args.pdb_file,args.proline)

    elif args.proline == 0:
        print("Both proline and non proline containing helices ".format(args))

        return(args.pdb_file,args.proline)
    else:
        print("Incorrect -p value (should be 0 1 or 2) selected")
        print(args)
        return(args.pdb_file,args.proline)


def linux_arguments():
    """
    uses argparse and returns the 1st argument as the imput file, and 2nd
    argument as the output file.
    Input: .sec file made using psbsecstr on a pdbfile
    output: .1hr or .2hr file dependend on the helix parser.
    returns: input_file,output_file

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'sec_file',
        help=
        'The pdbsecstr output of a pdb file that contains residues and secondary structure '
    )
    parser.add_argument(
        '1hr_file',
        help='output filename which will contain single alpha helicies residues'
    )

    # This creates a namespace object which allows you to treat files as if they are open
    args = parser.parse_args()

    sec_name = vars(args)['sec_file']
    onehr_name = vars(args)['1hr_file']
    return (sec_name, onehr_name)


def single_helix_parser(input_file, helicies_length=13):
    """

    :param input_file string:
    :param helicies_length:
    :return:one_helix_l,chains_sec_str_d,chains_res_no_d,chains_res_name_d

    examples:
            one_helix_l (helix in list format),
       [['A8', 'A9', 'A10', 'A11', 'A12', 'A13', 'A14'], ['A116', 'A117', 'A118', 'A119', 'A120']]
            chains_sec_str_d : dict
        {'A': '-----hHHHHHHHHHHHHHHHHHHHHHHHht--EEEEEeTTthHHHHHHHHHH'}
            chains_res_no_d
        {'A': ['A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10']}
            chains_res_name_d
        {'A': ['THR', 'GLY', 'ILE', 'THR', 'TYR', 'ASP', 'GLU', 'ASP', 'ARG', 'LYS', 'THR']}

    """

    res_no_l = []  # for residue names
    res_name_l = []  # for amino acid names
    sec_str_l = []  # for sec structure prediction

    one_helix_l = []  # contains one a list aminoacids (also a list)

    #text = fileread(input_file)
    text= input_file

    # Extracts the residue no, amino acid and secstr and signs to variables
    rx_seq = re.compile(r"^(\w+?)\s+?(\w+?)\s+?(\S)", re.MULTILINE)

    for match in rx_seq.finditer(text):
        res_no, res_name, sec_str = match.groups()

        res_no_l.append(res_no)
        res_name_l.append(res_name)
        sec_str_l += sec_str

    chains_sec_str_d = keychain_value_str(res_no_l, sec_str_l)

    chains_res_no_d = keychain_value_list(res_no_l, res_no_l)

    chains_res_name_d = keychain_value_list(res_no_l, res_name_l)

    # only adds if a proline is found in the gap
    # contains 2 groups, the 1st group being the whole helix and group 2 being the gap
    for x in chains_sec_str_d:
        # print(x)
        regex = "([H|h]{" + str(helicies_length) + ",})"
        p = re.compile(r"" + regex + "")

        # if one is found it prints out the residues numbers of that helix
        for match in p.finditer(chains_sec_str_d[x]):
            # adjusted to check for Proline around the gap 1 before and 1 after

            one_helix_l += [
                chains_res_no_d[x][(match.start(1)):(match.end(1))]
            ]
    return (one_helix_l, chains_sec_str_d, chains_res_no_d, chains_res_name_d)



def double_helix_parser(input_file, helicies_length=6, helix_gap=3, pro_eitherside=3):
    """ the purpose of this function is to create a list of residues names A11
    of proteins that are made of two helixes seperated by a gap
    pro_eitherside is how many side of the gap I should search

   :param input_file string:
    :param helicies_length:
    :return:two_helix_l,chains_sec_str_d,chains_res_no_d,chains_res_name_d

    examples:
            two_helix_l (helix in list format),
       [['A8', 'A9', 'A10', 'A11', 'A12', 'A13', 'A14'], ['A116', 'A117', 'A118', 'A119', 'A120']]
            chains_sec_str_d : dict
        {'A': '-----hHHHHHHHHHHHHHHHHHHHHHHHht--EEEEEeTTthHHHHHHHHHH'}
            chains_res_no_d
        {'A': ['A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10']}
            chains_res_name_d
        {'A': ['THR', 'GLY', 'ILE', 'THR', 'TYR', 'ASP', 'GLU', 'ASP', 'ARG', 'LYS', 'THR']}

    """
    res_no_l = []  # for residue names
    res_name_l = []  # for amino acid names
    sec_str_l = []  # for sec structure prediction


    two_helix_l = []  # contains one a list aminoacids (also a list)
    proline_secstr = {}
    proline_res_index = []

    # Extracts the residue no, amino acid and secstr and signs to variables
    rx_seq = re.compile(r"^(\w+?)\s+?(\w+?)\s+?(\S)", re.MULTILINE)
    #text = fileread(input_file)

    text = input_file

    # assign the matched groups in the text to the res_no_l, res_name_l and sec_str_str
    for match in rx_seq.finditer(text):

        res_no, res_name, sec_str = match.groups()

        res_no_l.append(res_no)
        res_name_l.append(res_name)
        sec_str_l += sec_str

    # creates dictionaries for each with the chain as the key
    chains_sec_str_d = keychain_value_str(res_no_l, sec_str_l)
    chains_res_no_d = keychain_value_list(res_no_l, res_no_l)
    chains_res_name_d = keychain_value_list(res_no_l, res_name_l)
    mod_chains_sec_str_d = chains_sec_str_d.copy()

    # which a Pro is found a in the res_name_d[chain] its secstr in sec_str_d is replaced with a P
    # We will then search for this P later on

    counter = 0
    for chain in chains_res_name_d:
        # print(chains_res_name_d[chain])
        counter = 0
        for residue in chains_res_name_d[chain]:
            # print(chains_res_name_d[chain][counter])

            # default is
            # #if residue == 'PRO':

            #If you want to recive any proline within secstr H (no gap)
            # if residue == 'PRO' and chains_sec_str_d[chain][counter] =='H':
            # Proline with secstr h (no gap)
            # if residue == 'PRO' and chains_sec_str_d[chain][counter] =='h':
            #or


            # and to capture only 1 gap
            #if residue == 'PRO' and (chains_sec_str_d[chain][counter] !='h' and chains_sec_str_d[chain][counter] != 'H'):

            if residue == 'PRO':

                #chains_sec_str_d[chain] = chains_sec_str_d[chain][:counter] + 'P' + chains_sec_str_d[chain][
                #                                                                    counter + 1:]

                mod_chains_sec_str_d[chain] = mod_chains_sec_str_d[chain][:counter] + 'P' + mod_chains_sec_str_d[chain][counter + 1:]
            counter += 1

            # only adds if a proline is found in the gap
    # contains 2 groups, the 1st group being the whole helix and group 2 being the gap

    for chain in mod_chains_sec_str_d:
        # The {1} either side of the P creates a 3 gap
        #  so it captures anything in a gap             regex = "([h|H]{5,}(?:.?){1}(P)(?:.?){1}[h|H]{5,})"
        # to capture ONLY three gap                     regex = "([h|H]{5,}[^Hh](P)[^Hh][h|H]{5,})"
        # to capture ONLY two gap                       regex = "([h|H]{5,}[^Hh](P)[h|H]{6,}|[h|H]{6,}(P)[^Hh][h|H]{5,})"
        # To capture ONLY 1 gap (only allowing proline) regex = "([h|H]{6,}(P)[h|H]{6,})"

        regex = "([h|H]{5,}[^Hh](P)[^Hh][h|H]{5,})"
        p = re.compile(r"" + regex + "")

        # if one is found it prints out the residues numbers of that helix
        for match in p.finditer(mod_chains_sec_str_d[chain]):
            # adjusted to check for Proline around the gap 1 before and 1 after
            two_helix_l += [chains_res_no_d[chain][(match.start(1)): (match.end(1))]]
            # if chain in proline_res_index:
            #     proline_res_index[chain] = proline_res_index[chain].append(match.start(2))
            # else:
            #     proline_res_index[chain] = match.start(2)

            temp_list = list((chains_res_no_d[chain][match.start(2)], chains_sec_str_d[chain][match.start(2)]))

            proline_res_index.append(temp_list)
            #proline_secstr[chain] =   (chains_sec_str_d[chain][match.start(2)],proline_res_index)

    return(two_helix_l,mod_chains_sec_str_d,chains_res_no_d,chains_res_name_d, proline_res_index)

def ca_atom_organiser(aa_chains, ca_residues, output_file):
    # opens aa_chains file and splits by blank line and line
    file_list=[]
    aa_chain_residues = fileread(aa_chains)

    first_aa_in_heli = aa_chains_split(aa_chain_residues)

    # opens file of alpha carbons and splits by line
    ca_chain_string = fileread(ca_residues)
    ca_chain_list = ca_chain_string.splitlines()

    # if there is a non value for first_residue_pdblines end the function
    # It is a Dictionary of lists, each key the first residue no of a helix, each value
    # ca atoms in that helix
    pdb_output = first_residue_pdblines(first_aa_in_heli, ca_chain_list)
    write_helix_dict(output_file,pdb_output)



    return (pdb_output)


def write_helix_dict(output_file, ca_atom_of_helix_dict):
    tempstring = ""
    output = open(output_file, 'w')

    # This gives the pdb file name at the start of the document
    # pdbfilename = aa_chains[:4]
    # tempstring += pdbfilename

    for key in sorted(ca_atom_of_helix_dict):
        # tempstring +="\n"
        # tempstring +="\n" + key +"\n"
        tempstring += key + "\n"
        for value in ca_atom_of_helix_dict[key]:
            tempstring += value
        tempstring += "\n"
    output.write(tempstring)
    # print(tempstring)

    # print(tempstring)
    output.close()


def first_residue_pdblines(aa_list, pdb_ca_list):
    """
     seachres the chain list for the starting residues, when found fills a dict
     with key the first residue in a helix and the value being the pdb file lines
     of associated residues
    """

    list_of_residues = []
    location = []
    temp_dict = {}
    temp_dict2 = {}

    # captures the residue numbers with chain for pdb format
    pattern = re.compile(r'^ATOM\s+?\d+?\s+?CA\s+?\w+?\s(\w+?\s*?\d+?)\s')
    counter = 0

    # This creates a list with the index locations of each item in aa_list
    # in ca_list
    for line in pdb_ca_list:

        if pattern.search(line):
            hit = pattern.search(line)

            list_of_residues += [hit.group(1)]

    # if the file has no fist residues in a helix an empty list is returned.
    # if a res no from aa_list is in list of residues (from pdb CA list)
    # note the pdb CA line index in location
    for x in aa_list:
        if x == '':
            return
        if x == ['']:
            return
        if x in list_of_residues:
            location.append(list_of_residues.index(x))
        else:
            return
    # This reverses the AA list (Turns [A 7, A 20, A 118] to [A 118, A20, A7]
    # It captures those CA residues up to its location and then deletes them.

    location.reverse()
    counter = 0

    for x in reversed(aa_list):
        location[counter]
        temp_dict[x] = pdb_ca_list[location[counter]:]
        del pdb_ca_list[location[counter]:]
        counter += 1

    # This creates a new dict to reverses the order of the 1st dict
    # which currently goes last amino acid to first

    for key in sorted(temp_dict):
        temp_dict2[key] = temp_dict[key]
    return (temp_dict2)


def aa_chains_split(chains):
    """
    Splits the lists of chains by empty newline and then takes the
    first amino acid from each chain for naming and organising

    """

    chains = chains.lstrip()
    chains = chains.rstrip()

    myarray = chains.split("\n\n")

    temparray = []
    temparray2 = []
    temparray3 = []

    for x in myarray:
        y = x.split("\n")
        temparray += [y]

    # temparray = list(filter(None, temparray))

    for x in temparray:
        temparray2 += [x[0]]

    for x in temparray2:

        # these spaces are needed to mimic the space in pdb files between chain
        # and res number
        if (len(x[1:])) == 1:
            y = x[0] + "   " + x[1:]
        if (len(x[1:])) == 2:
            y = x[0] + "  " + x[1:]
        if (len(x[1:])) == 3:
            y = x[0] + " " + x[1:]
        if (len(x[1:])) >= 4:
            y = x[0] + "" + x[1:]
        temparray3 += [y]

    temparray3 = sorted(temparray3)

    return (temparray3)

def filewrite_nestedlist(filename, outerlist):
    """
    input:
    filename:name of file to be written to
    outerlist: list of lists
    output: a file named filename, with a /m seperated list of lists
    This writes out the elements of a list of lists in order with a blank linke between each set of list elements and
    each list element on a newline
    """
    tempstr = ""

    for nestedlist in outerlist:
        for element in nestedlist:
            tempstr += (element + "\n")
        tempstr += ("\n")

    with open(filename, 'w') as out_file:
        out_file.write(tempstr)

def string_nestedlist(outerlist):
    """
    input:
    outerlist: list of lists
    output: returns a string with a /m seperated list of lists
    This writes out the elements of a list of lists in order with a blank linke between each set of list elements and
    each list element on a newline
    """
    tempstr = ""

    for nestedlist in outerlist:
        for element in nestedlist:
            tempstr += (element + "\n")
        tempstr += ("\n")

    return(tempstr)


def keychain_value_str(key_provider, dict_values):
    """
    creates a dictionary where the list key providesr first nonwhitespace is
     used as the key in this case it is the chain of the residue number. The
      key value is made into string associated with each chain.
    :param key_provider: takes the first character from key provider
    :param dict_values: either res_no,pdbsecstr,res_type
    :return:
    """
    counter = 0
    new_dict = {}
    for x in key_provider:
        #   this checks if a chain exists and then adds secstr to a dictoinary
        # of that chain letter for : secondary structutre, residue number and residue name
        if x[0] in new_dict:
            new_dict[x[0]] += dict_values[counter]
            counter += 1
        else:
            new_dict[x[0]] = dict_values[counter]
            counter += 1
    return new_dict


def keychain_value_list(key_provider, dict_values):
    counter = 0
    new_dict = {}
    for x in key_provider:
        if x[0] in new_dict:
            new_dict[x[0]] += [dict_values[counter]]
            counter += 1
        else:
            new_dict[x[0]] = [dict_values[counter]]
            counter += 1
    return new_dict

def fileread(filename):
    file = open(filename, 'r')
    output = file.read()
    file.close()
    return (output)

# Functions below were from proline_non_proline_middle_angle
def middle_angle_linux_arguments():
    """
    This determins the arguments for the program when in a linux enviroment
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('format_file', help='ca atoms of proteins seperated by chains')
    parser.add_argument('angle_file', help='output file which will have the pdb name and 1st residue with ABC angle')
    parser.add_argument('helix_type', type=int, choices=[1, 2], default='1',
                        help='specifies if it is a proline or non proline helix')
    parser.add_argument('--pdbline', action='store_true',
                        help='also creates pdbfile containing pdblines used in bend angle calculation')

    # This creates a namespace object which allows you to treat files as if they are open
    args = parser.parse_args()

    format_name = vars(args)['format_file']
    angle_name = vars(args)['angle_file']
    helix_type = vars(args)['helix_type']

    if args.pdbline:
        print('Pdbline creation specified')

        pdbline_option = True
        print(str(pdbline_option))
    else:
        pdbline_option = False

    return (format_name, angle_name, helix_type, pdbline_option)


def residue_pairs_for_pdbline(input_file, helix_type,proline_res_d):
    """
    input: string
    the file name

    output: list of objects


    this calls the input files and outputs a list of helicies.
    it calls first_pro_last_finder() and information_extracter()
    to determine the 3 residues and their attributes

    it calls either proline/non_proline_segment_first_mid_last_finder depending on the
    helix type
    """

    file_txt = fileread(input_file)

    helix_start_mid_end = []

    first_res_no = []
    last_res_no = []
    whole_helix_single_residues = []
    pdbline_single_residues = []
    counter = 0
    # match an object where the data is broken up by line breaks and chain res no linebreak
    # captures 3 groups per helix: the first resno, the CA info for helix, The last resno
    # The  usees lookahead to capture the CA inf while also capturing the last resno
    pattern = re.compile(r'(\w+\s*?\d+?)\n(?=(.*))(?:.*?)(?:(\w\s*?\d+?)\s+?(?:-*?\d+?\.\d+\s*?){1,}C\s+?)\n')
    match = pattern.findall(file_txt)

    # Each x is a helix. This helix contains the ca pdb atoms of that helix
    for helix_residues in match:

        temp_list = []
        first_res_no.append(helix_residues[0])
        last_res_no.append(helix_residues[2])
        pdb_ca = helix_residues[1]

        # has the first and last residue which will be used to
        # create the pdblines. if no proline is found that helix is skipped

        # Positions are the index locations of the helix segment used for calculating
        # bend angle made

        if helix_type == 1:
            positions = non_proline_segment_first_mid_last_finder(helix_residues[1])
            #print("non proline")

        elif helix_type == 2:
            positions = proline_segment_first_mid_last_finder(helix_residues[1],proline_res_d[counter])
            #print("proline")

        if positions == None:
            break

        atom_inf = lineinformation_extractor(positions, pdb_ca)
        for helix_residues in atom_inf:
            temp_list += [helix_residues]

        helix_start_mid_end.append(temp_list)

        whole_helix_single_residues.append(residue_extractor_from_ca_pdb(pdb_ca))

        pdbline_single_residues.append(whole_helix_single_residues[0][positions[0]:positions[2] + 1])

        counter += 1
    # print("The index of the first residues is " +str(positions[0]) +" and "+ str(positions[2]+1))

    return (helix_start_mid_end, first_res_no, last_res_no, whole_helix_single_residues, pdbline_single_residues)


def proline_segment_first_mid_last_finder(protein_pdb, proline_information):
    """
    input: a string

    output: a tupple
    containing 3 integers

    this finds the position that will be used for pdbline for creating
    lines of best fit. The first residue of the first pdbline, the proline
    at the middle of the second pdbline and the last residue of the third
    pdb line.

    """
    res_p = re.compile(r'ATOM\s+?\d+?\s+?\w+?\s+?(\w+?)\s(\w\s*?\w+?)\s')
    res = res_p.findall(protein_pdb)

    counter = 0
    gap_window = int(0.5 * len(res)) - 6
    pro_res = None

    while counter <= (gap_window):

        midpoint = int((len(res) - 1) / 2)
        current_resno = res[midpoint - counter][1].replace(" ", "")
        proline_resiude_no = proline_information[0]
        if current_resno == proline_resiude_no:



            pro_res = (midpoint - counter)
            counter += 1
            break

        elif res[midpoint + counter][0] == "PRO":

            pro_res = (midpoint + counter)
            counter += 1
            break

        else:

            counter += 1

    if pro_res == None:

        return (None)

    return (pro_res - 6, pro_res, pro_res + 6)


def lineinformation_extractor(selected_ca, pdb_txt):
    """
    input: Tupple
    containing the positions of the mid -6 /PRO/mid +6 residue in the
    aminoacid sequence (for calculating the bend angle between them)

    output:a list of 3 lists,
    each containing the information for first/mid/last res

    This extracts the the aminoacid, resno and the xyz cords for creating the atom
    objects. The
    """
    first = selected_ca[0]
    middle = selected_ca[1]
    last = selected_ca[2]

    first_six = ""
    middle_five = ""
    last_six = ""

    cord_p = re.compile(
        r'ATOM\s+?\d+?\s+?CA\s+?(\w+?)\s(\w\s*?\d+?)\s+?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)')

    cords = cord_p.findall(pdb_txt)

    temp_str = (str(cords[first][1]))
    temp_str = temp_str.replace(" ", "")
    first_six += temp_str
    first_six += (" ")
    temp_str = (str(cords[first + 5][1]))
    temp_str = temp_str.replace(" ", "")
    first_six += temp_str

    temp_str = (str(cords[middle - 2][1]))
    temp_str = temp_str.replace(" ", "")
    middle_five += temp_str
    middle_five += (" ")
    temp_str = (str(cords[middle + 2][1]))
    temp_str = temp_str.replace(" ", "")
    middle_five += temp_str

    temp_str = (str(cords[last - 5][1]))
    temp_str = temp_str.replace(" ", "")
    last_six += temp_str
    last_six += (" ")
    temp_str = (str(cords[last][1]))
    temp_str = temp_str.replace(" ", "")
    last_six += temp_str
    return [first_six, middle_five, last_six]


def shell_interface(residue_pos, input_file, pdbline_option):
    """
    The purpose of this function is to wrap some residue chain and residue numbers with a wrapper that
    calls the pdbline in the command line

    input: list of 3 strings 		e.g ['E115 E120', 'E132 E136', 'E148 E153']

    each string is the start and end of the resdidues involved in pdbline. The first ist the starting pdbline,
    the second is the middle pdbline and the third is the end pdb line

    output: tupple of 3 strings
    e.g. ('pdbline A52 A57 1m1j.pdb 1m1j_line_s.pdb', 'pdbline A64 A68 1m1j.pdb 1m1j_line_m.pdb', 'pdbline A74 A79 1m1j.pdb 1m1j_line_e.pdb')

    each of these are a string which will be called on the command line

    """
    print("##################################")

    filename = input_file

    start_helix = commandline_wrapper(residue_pos[0], filename)
    start_helix_str = create_pdbline_string(start_helix)
    start_pdbline_midpoint = calculate_midpoint(start_helix_str)

    mid_helix = commandline_wrapper(residue_pos[1], filename)
    mid_helix_str = create_pdbline_string(mid_helix)
    middle_pdbline_midpoint = calculate_midpoint(mid_helix_str)

    end_helix = commandline_wrapper(residue_pos[2], filename)
    end_helix_str = create_pdbline_string(end_helix)
    end_pdbline_midpoint = calculate_midpoint(end_helix_str)

    """

        # This is added as part of the optional commandline argument --pdbline, if so it writes to a new file for every
        # is based on the horrible global variable pdbline_option which is true if argument "--pdbline" is given

        # input variables;
         residue_pos (list of lists coordiantes) [['A14 A19', 'A18 A22', 'A21 A26']
         , ['A117 A122', 'A121 A125', 'A124 A129']
         , ['A156 A161', 'A160 A164', 'A163 A168'],
          ['A199 A204', 'A203 A207', 'A206 A211']]

          filename : string "1ct5"
          start_helix_str contains the atom information from pdbline for that residue pair
        """
    if pdbline_option == True:
        print('roger roger')
        temp_res_list = []
        for pair in residue_pos:
            temp_res_list += pair.split()
        temp_filename = (str(filename) + '_' + str(temp_res_list[0] + '_' + str(temp_res_list[-1])) + '_pdbline')

        # Accesses original .pdb filename
        pdbline_file_list = []
        filename += '.pdb'
        # pdbline_file_list.append(fileread(filename))

        # this is required to remove b'string' from start/mid/end_helix_str which has been converted over when it was in byte form
        start_helix_str = start_helix_str.strip('b')
        start_helix_str = start_helix_str.strip("''")
        # start_helix_str = end_helix_str.split('\n')

        mid_helix_str = mid_helix_str.strip('b')
        mid_helix_str = mid_helix_str.strip("''")
        # mid_helix_str = end_helix_str.split('\n')

        end_helix_str = end_helix_str.strip('b')
        end_helix_str = end_helix_str.strip("''")
        # end_helix_str = end_helix_str.split('\n')
        # print(end_helix_str)

        # start_helix_str =start_helix_str.decode('utf-8')

        # this is all an attempts to combine into one list and to replace '\\n' with '\n',
        # i could have just used find and replace
        pdbline_file_list.append(start_helix_str)
        pdbline_file_list.append(mid_helix_str)
        pdbline_file_list.append(end_helix_str)
        # print(pdbline_file_list)
        pdbline_file_str = ('').join(pdbline_file_list)
        pdbline_file_list = pdbline_file_str.split('\\n')

        # writes the awakward list format to a file nicly
        with open(temp_filename, 'w') as text_file:
            text_file.writelines("%s\n" % line for line in pdbline_file_list)

    # text_file.write(str(pdbline_file_str))

    # print(temp_filename)
    return (start_pdbline_midpoint, middle_pdbline_midpoint, end_pdbline_midpoint)


def commandline_wrapper(res_pair, input_name):
    temp_str = ""

    temp_str += "pdbline"
    temp_str += " "
    temp_str += str(res_pair)
    temp_str += " "
    temp_str += input_name
    temp_str += ".pdb"

    print(temp_str)
    return (temp_str)


def create_pdbline_string(string):
    retval = subprocess.check_output(string, shell=True)
    retval = str(retval)  # Convert from byte string
    return (retval)


def calculate_midpoint(pdbline_str):
    """ This calculate the bend angle by measureing pdblines created at the tips and
    the middle
    """

    pdbline = re.compile(
        r'ATOM\s+?(\d+?)\s+?X\s+?\w+?\s\w\s*?\d+?\s+?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)')
    pdbline_xyz = pdbline.findall(pdbline_str)

    # print(" .line file findall :",pdbline_xyz)
    # print("##################")
    mid = int((len(pdbline_xyz) - 1) / 2)

    # print(pdbline_xyz[mid])
    # print(mid)

    x = pdbline_xyz[mid][1]

    y = pdbline_xyz[mid][2]

    z = pdbline_xyz[mid][3]
    # print('x,y,z :',x,y,z)
    return (x, y, z)


def calculate_angle(first_xyz, second_xyz, third_xyz):
    first_xyz = list(first_xyz)
    second_xyz = list(second_xyz)
    third_xyz = list(third_xyz)

    float_first = [float(i) for i in first_xyz]
    float_second = [float(i) for i in second_xyz]
    float_third = [float(i) for i in third_xyz]

    a = np.array([float_first[0], float_first[1], float_first[2]])
    b = np.array([float_second[0], float_second[1], float_second[2]])
    c = np.array([float_third[0], float_third[1], float_third[2]])

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    angle = (np.degrees(angle))
    return (angle)


def segment_middle_res(pdb_line_res_pair):
    """
    finds the proline res by taking the pdbline res pair e.g. A174 A178,
    takes the first and adds two to get A176
    """

    res_pair = pdb_line_res_pair
    res_pair = res_pair.split()

    first_res = res_pair[0]

    first_res_no = first_res[1:]
    first_res_chain = first_res[:1]

    proline_res_no = int(first_res_no) + 2
    proline_res = (first_res_chain + str(proline_res_no))

    return (proline_res)

    def residue_extractor_from_ca_pdb(helix_pdb_info):
        helix_pdb_info


def residue_extractor_from_ca_pdb(ca_string):
    single_letter_res = []
    residue_finder = re.compile(r'ATOM\s+?\d+?\s+?\w+?\s+?(\w+?)\s')
    residues = residue_finder.findall(ca_string)

    three_letter_res_d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                          'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                          'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                          'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    for residue in residues:
        single_letter_res.append(three_letter_res_d[residue])
    single_letter_res = "".join(single_letter_res)
    return (single_letter_res)


def non_proline_segment_first_mid_last_finder(protein_pdb):
    """
    input: a string

    output: a tupple
    containing 3 integers

    this finds the position that will be used for pdbline for creating
    lines of best fit. The first residue of the first pdbline, the proline
    at the middle of the second pdbline and the last residue of the third
    pdb line.

    """
    res_p = re.compile(r'ATOM\s+?\d+?\s+?\w+?\s+?(\w+?)\s')
    res = res_p.findall(protein_pdb)

    mid = int((len(res) - 1) / 2)

    return (mid - 6, mid, mid + 6)

def proline_and_nonproline_middle_angle(format_file, angle_file, helix_type,proline_res_d, pdbline_option):


    """the purpose of this function is to create a list of residues names A11
    of proteins that are made of two helixes seperated by a gap
    pro_eitherside is how many side of the gap I should search
    1h3l A28 A28 A58 100.71272063185602
    """
    tempstring = ""
    pdbname = str(format_file)[:4]

    pdbline_res = residue_pairs_for_pdbline(format_file, helix_type,proline_res_d)
    print("Checking helix type  %s" % helix_type)

    helix_counter = 0
    for pdbline_segment in pdbline_res[0]:

        # I would insert the extracting_ca_inf_with_first_res.py scratch here
        #print("x or pdbline{0] is:")
        #print(pdbline_segment)
        center_pdbline_segment = segment_middle_res(pdbline_segment[1])

        one, two, three = shell_interface(pdbline_segment, pdbname, pdbline_option)
        angle = calculate_angle(one, two, three)

        # pdbline_res takes the residue for the first pdbline and splits by space, giving the first
        whole_helix_first_res = str(pdbline_res[1][helix_counter])
        whole_helix_first_res = whole_helix_first_res.replace(" ", "")
        whole_helix_last_res = str(pdbline_res[2][helix_counter])
        whole_helix_last_res = whole_helix_last_res.replace(" ", "")

        pdbline_segment_first_res = pdbline_segment[0].split()[0]
        pdbline_segment_last_res = pdbline_segment[2].split()[1]

        # print("pdbline_segment first res is " + str(pdbline_segment_first_res) + "with last res " + pdbline_segment_last_res)

        whole_helix_single_residues = pdbline_res[3][helix_counter]
        pdbline_segment_residues = pdbline_res[4][helix_counter]

        tempstring += format_file[:4] + " " + str(whole_helix_first_res) + " " + str(whole_helix_last_res) + " " \
                      + center_pdbline_segment + " " + whole_helix_single_residues + " " + " " + str(
            angle) + " " + pdbline_segment_first_res + " " + pdbline_segment_last_res + " " + pdbline_segment_residues +\
                      " "
        if helix_type == 2:
            tempstring += " " +  str(proline_res_d[helix_counter][1]) +"\n"
        else:
            tempstring +="\n"
        helix_counter += 1

    print(tempstring)
    file = open(angle_file, 'w')
    file.write(tempstring)
    file.close()

