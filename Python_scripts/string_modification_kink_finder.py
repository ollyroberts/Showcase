#!/usr/bin/python3

import subprocess
"""
Lets seperate this into blocks
i)  I wish to grab all rwos with the same pdb file
ii) Extract the res number pairs
iii)run the multiple lines
vi) parse output for correct files
"""
"""

i)

1cza N414 N435 N421 GSLYKTHPQYSRRFHKTLRRLV  110.72633798240781 N415 N427 SLYKTHPQYSRRF
1cza N862 N883 N869 GTLYKLHPHFSRIMHQTVKELS  115.0510975191965 N863 N875 SLYKTHPQYSRRF
1mpx A358 A377 A370 DTARQFRHDVLRPFFDQYLV  170.41578509514642 A364 A376 RHDVLRPFFDQYL
1mpx B358 B377 B370 DTARQFRHDVLRPFFDQYLV  163.26882219236043 B364 B376 RHDVLRPFFDQYL
1mpx C358 C377 C370 DTARQFRHDVLRPFFDQYLV  169.75858964683323 C364 C376 RHDVLRPFFDQYL
1mpx D358 D377 D370 DTARQFRHDVLRPFFDQYLV  169.70217884633593 D364 D376 RHDVLRPFFDQYL
2qsf X279 X297 X288 LLENISARYPQLREHIMAN  125.02727967075862 X282 X294 NISARYPQLREHI
2qss B49 B74 B57 TADAVMNNPKVKAHGKKVLDSFSNGM  126.29398428495547 B51 B63 DAVMNNPKVKAHG
2qss D49 D75 D57 TADAVMNNPKVKAHGKKVLDSFSNGMK  130.805723901276 D51 D63 DAVMNNPKVKAHG
3f7m A326 A345 A332 GTSMATPHIAGLAAYLFGLE  172.38848852332777 A326 A338 GTSMATPHIAGLA
3tui A5 A41 A33 MMWLLVRGVWETLAMTFVSGFFGFVIGLPVGVLLYVT  166.23373611484138 A27 A39 GFVIGLPVGVLLY
3tui A49 A76 A67 NAKLYRTVSAIVNIFRSIPFIILLVWMI  128.944465353281 A61 A73 LAMTFVSGFFGFV
3tui A88 A101 A95 GLQAAIVPLTVGAA  158.25144769124705 A89 A101 MWLLVRGVWETLA
3tui B5 B42 B33 MMWLLVRGVWETLAMTFVSGFFGFVIGLPVGVLLYVTR  172.2003074438012 B27 B39 GFVIGLPVGVLLY
3tui B49 B76 B67 NAKLYRTVSAIVNIFRSIPFIILLVWMI  130.79059197004446 B61 B73 LAMTFVSGFFGFV
3tui B88 B101 B95 GLQAAIVPLTVGAA  170.58423113754733 B89 B101 MWLLVRGVWETLA
3tui D275 D290 D281 GQSVDAPLLSETARRF  159.36102036930097 D275 D287 MMWLLVRGVWETL
3tui E5 E41 E33 MMWLLVRGVWETLAMTFVSGFFGFVIGLPVGVLLYVT  175.24690658954773 E27 E39 GFVIGLPVGVLLY
3tui E49 E76 E67 NAKLYRTVSAIVNIFRSIPFIILLVWMI  128.7884633010175 E61 E73 LAMTFVSGFFGFV
3tui E88 E101 E95 GLQAAIVPLTVGAA  162.16789187922018 E89 E101 MWLLVRGVWETLA
3tui F5 F41 F33 MMWLLVRGVWETLAMTFVSGFFGFVIGLPVGVLLYVT  170.50291121093792 F27 F39 GFVIGLPVGVLLY
3tui F49 F76 F67 NAKLYRTVSAIVNIFRSIPFIILLVWMI  127.27364651355367 F61 F73 LAMTFVSGFFGFV
3tui F88 F101 F95 GLQAAIVPLTVGAA  166.65984944083547 F89 F101 MWLLVRGVWETLA
3tui G275 G291 G281 GQSVDAPLLSETARRFN  151.09335744975058 G275 G287 MMWLLVRGVWETL
4jkz A116 A141 A130 HEDTAARLREGFSRPLIESVRDRLRD  168.14334031212854 A124 A136 REGFSRPLIESVR
4jl5 A31 A56 A44 TGDILREAVQKGTPLGKKAKEYMERG  98.61101124200661 A38 A50 AVQKGTPLGKKAK
4jl5 A156 A178 A171 EVIKKRLEVYREQTAPLIEYYKK  170.00500703766653 A165 A177 QKGTPLGKKAKEY
4jl5 B31 B54 B44 TGDILREAVQKGTPLGKKAKEYME  99.09266785107177 B38 B50 AVQKGTPLGKKAK
4jl5 B156 B178 B171 EVIKKRLEVYREQTAPLIEYYKK  173.9415848414965 B165 B177 QKGTPLGKKAKEY


ii) input:

3tui A5 A41 A33 MMWLLVRGVWETLAMTFVSGFFGFVIGLPVGVLLYVT  166.23373611484138 A27 A39 GFVIGLPVGVLLY
3tui A49 A76 A67 NAKLYRTVSAIVNIFRSIPFIILLVWMI  128.944465353281 A61 A73 LAMTFVSGFFGFV
3tui A88 A101 A95 GLQAAIVPLTVGAA  158.25144769124705 A89 A101 MWLLVRGVWETLA
3tui B5 B42 B33 MMWLLVRGVWETLAMTFVSGFFGFVIGLPVGVLLYVTR  172.2003074438012 B27 B39 GFVIGLPVGVLLY
3tui B49 B76 B67 NAKLYRTVSAIVNIFRSIPFIILLVWMI  130.79059197004446 B61 B73 LAMTFVSGFFGFV
3tui B88 B101 B95 GLQAAIVPLTVGAA  170.58423113754733 B89 B101 MWLLVRGVWETLA
3tui D275 D290 D281 GQSVDAPLLSETARRF  159.36102036930097 D275 D287 MMWLLVRGVWETL

    output:

    var2  ["3tui","A","27 39", "61 73", "89 101",]
    var3  ["3tui","B","27 39", "61 73", "89 101", "275 287"]

iii)
res_extractor
python2.7 ~/Downloads/KF_err_lin/Kink_Finder.py -f ~/Git/Project/refactored_code/kink_finder_tests/3tui.pdb -o ~/Git/Project/refactored_code/kink_finder_tests/3tui.pdb_output/ -l '27-39 61-73 89-101 27-39 61-73 89-101 275-287' -d"

iv)
in folder 3tui.pdb_output/angles.csv
~/Git/Project/refactored_code/kink_finder_tests/3tui.pdb_output/angles.csv
capture splitline[0]  splitline[7]
A27

3tuiA27,0,0,0,0,0,11.756,8.262,0,0,0,0,0,0
3tuiA61,0,0,0,0,0,101.500,102.289,0,0,0,0,0,0
3tuiA89,0,0,0,0,0,20.402,20.514,0,0,0,0,0,0
3tuiA27,0,0,0,0,0,11.756,8.262,0,0,0,0,0,0
3tuiA61,0,0,0,0,0,101.500,102.289,0,0,0,0,0,0
3tuiA89,0,0,0,0,0,20.402,20.514,0,0,0,0,0,0
3tuiA89,0,0,0,0,0,20.402,20.514,0,0,0,0,0,0
3tuiB27,0,0,0,0,0,14.052,9.205,0,0,0,0,0,0
3tuiB61,0,0,0,0,0,98.530,98.464,0,0,0,0,0,0
3tuiB89,0,0,0,0,0,18.067,17.529,0,0,0,0,0,0
3tuiB27,0,0,0,0,0,14.052,9.205,0,0,0,0,0,0
3tuiB61,0,0,0,0,0,98.530,98.464,0,0,0,0,0,0
3tuiB89,0,0,0,0,0,18.067,17.529,0,0,0,0,0,0
3tuiB89,0,0,0,0,0,18.067,17.529,0,0,0,0,0,0
3tuiC27,0,0,0,0,0,35.818,15.246,0,0,0,0,0,0
3tuiC61,0,0,0,0,0,114.796,84.070,0,0,0,0,0,0
3tuiC89,0,0,0,0,0,92.479,107.820,0,0,0,0,0,0
3tuiC27,0,0,0,0,0,35.818,15.246,0,0,0,0,0,0
3tuiC61,0,0,0,0,0,114.796,84.070,0,0,0,0,0,0
3tuiC89,0,0,0,0,0,92.479,107.820,0,0,0,0,0,0
3tuiC275,0,0,0,0,0,54.444,38.479,0,0,0,0,0,0
3tuiD27,0,0,0,0,0,44.166,15.404,0,0,0,0,0,0
3tuiD61,0,0,0,0,0,114.657,83.712,0,0,0,0,0,0
3tuiD89,0,0,0,0,0,93.025,103.236,0,0,0,0,0,0
3tuiD27,0,0,0,0,0,44.166,15.404,0,0,0,0,0,0
3tuiD61,0,0,0,0,0,114.657,83.712,0,0,0,0,0,0
3tuiD89,0,0,0,0,0,93.025,103.236,0,0,0,0,0,0
3tuiD275,0,0,0,0,0,53.538,36.553,0,0,0,0,0,0
3tuiE27,0,0,0,0,0,14.536,9.574,0,0,0,0,0,0
3tuiE61,0,0,0,0,0,98.076,97.801,0,0,0,0,0,0
3tuiE89,0,0,0,0,0,22.013,24.342,0,0,0,0,0,0
3tuiE27,0,0,0,0,0,14.536,9.574,0,0,0,0,0,0
3tuiE61,0,0,0,0,0,98.076,97.801,0,0,0,0,0,0
3tuiE89,0,0,0,0,0,22.013,24.342,0,0,0,0,0,0
3tuiE89,0,0,0,0,0,22.013,24.342,0,0,0,0,0,0
3tuiF27,0,0,0,0,0,14.467,8.486,0,0,0,0,0,0
3tuiF61,0,0,0,0,0,97.362,99.188,0,0,0,0,0,0
3tuiF89,0,0,0,0,0,21.659,21.206,0,0,0,0,0,0
3tuiF27,0,0,0,0,0,14.467,8.486,0,0,0,0,0,0
3tuiF61,0,0,0,0,0,97.362,99.188,0,0,0,0,0,0
3tuiF89,0,0,0,0,0,21.659,21.206,0,0,0,0,0,0
3tuiF89,0,0,0,0,0,21.659,21.206,0,0,0,0,0,0
3tuiG27,0,0,0,0,0,117.876,133.768,0,0,0,0,0,0
3tuiG61,0,0,0,0,0,153.126,102.596,0,0,0,0,0,0
3tuiG89,0,0,0,0,0,94.166,107.852,0,0,0,0,0,0
3tuiG27,0,0,0,0,0,117.876,133.768,0,0,0,0,0,0
3tuiG61,0,0,0,0,0,153.126,102.596,0,0,0,0,0,0
3tuiG89,0,0,0,0,0,94.166,107.852,0,0,0,0,0,0
3tuiG275,0,0,0,0,0,53.339,37.997,0,0,0,0,0,0
3tuiH27,0,0,0,0,0,125.865,23.567,0,0,0,0,0,0
3tuiD89,0,0,0,0,0,93.025,103.236,0,0,0,0,0,0
3tuiD275,0,0,0,0,0,53.538,36.553,0,0,0,0,0,0
3tuiE27,0,0,0,0,0,14.536,9.574,0,0,0,0,0,0
3tuiE61,0,0,0,0,0,98.076,97.801,0,0,0,0,0,0
3tuiE89,0,0,0,0,0,22.013,24.342,0,0,0,0,0,0
3tuiE27,0,0,0,0,0,14.536,9.574,0,0,0,0,0,0
3tuiE61,0,0,0,0,0,98.076,97.801,0,0,0,0,0,0
3tuiE89,0,0,0,0,0,22.013,24.342,0,0,0,0,0,0
3tuiE89,0,0,0,0,0,22.013,24.342,0,0,0,0,0,0
3tuiF27,0,0,0,0,0,14.467,8.486,0,0,0,0,0,0
3tuiF61,0,0,0,0,0,97.362,99.188,0,0,0,0,0,0
3tuiF89,0,0,0,0,0,21.659,21.206,0,0,0,0,0,0
3tuiF27,0,0,0,0,0,14.467,8.486,0,0,0,0,0,0
3tuiF61,0,0,0,0,0,97.362,99.188,0,0,0,0,0,0
3tuiF89,0,0,0,0,0,21.659,21.206,0,0,0,0,0,0
3tuiF89,0,0,0,0,0,21.659,21.206,0,0,0,0,0,0
3tuiG27,0,0,0,0,0,117.876,133.768,0,0,0,0,0,0
3tuiG61,0,0,0,0,0,153.126,102.596,0,0,0,0,0,0
3tuiG89,0,0,0,0,0,94.166,107.852,0,0,0,0,0,0
3tuiG27,0,0,0,0,0,117.876,133.768,0,0,0,0,0,0
3tuiG61,0,0,0,0,0,153.126,102.596,0,0,0,0,0,0
3tuiG89,0,0,0,0,0,94.166,107.852,0,0,0,0,0,0
3tuiG275,0,0,0,0,0,53.339,37.997,0,0,0,0,0,0
3tuiH27,0,0,0,0,0,125.865,23.567,0,0,0,0,0,0
3tuiH61,0,0,0,0,0,105.761,78.980,0,0,0,0,0,0
3tuiH89,0,0,0,0,0,92.214,107.109,0,0,0,0,0,0
3tuiH275,0,0,0,0,0,63.794,43.708,0,0,0,0,0,0


v)
create a new  file which takes the original and add the .angle


two_helix_res_input ="1cza N414 N435 N421 GSLYKTHPQYSRRFHKTLRRLV  110.72633798240781 N415 N427 SLYKTHPQYSRRF"


"python2.7 ~/Downloads/KF_err_lin/Kink_Finder.py -f ~/Git/Project/refactored_code/kink_finder_tests/2angle/1cza.pdb -o ~/Git/Project/refactored_code/kink_finder_tests/2angle/1cza.pdb_output/ -l '415-427' -d"
target_output =  "python2.7 ~/Downloads/KF_err_lin/Kink_Finder.py -f " \
                 "~/Git/Project/refactored_code/kink_finder_tests/2angle/1cza.pdb -o " \
                 "~/Git/Project/refactored_code/kink_finder_tests/2angle/1cza.pdb_output/ -l '415-427' -d"
#print("Target output")
#print(target_output)


split_output = two_helix_res_input.split()
pdb_file, chain, first_pair_res, second_pair_res =(split_output[0] + ".pdb",split_output[6][:1], split_output[6][1:],split_output[7][1:])
#print(pdb_file, chain, first_pair_res, second_pair_res)

pdb_file ="1cza.pdb"
residues ="'415-427'"
chain = None

folder_directory = "~/Git/Project/refactored_code/kink_finder_tests/2angle/"

kink_finder_str = "python2.7 ~/Downloads/KF_err_lin/Kink_Finder.py"
pdb_file_location_str = "-f "+ folder_directory + pdb_file
kink_output_location_str = "-o "+ folder_directory +"/1cza.pdb_output/"
residues_display_str ="-l "+ residues +" -d"


list_of_strings=[kink_finder_str, pdb_file_location_str, kink_output_location_str, residues_display_str]
final_string = " ".join(list_of_strings)
#print(final_string)
"""


def pdb_res_pair(angle_imput_string):
    """

    :param angle_imput_string:
    :return: dict with pdb as key, list of list containing values
    """
    with open(angle_imput_string,'r') as openfile:
        file_contents = openfile.read()

    all_lines = file_contents.splitlines()
    counter = 0
    pdb_dict={}

    for line in all_lines:
        line = line.strip("\n")

    for line in all_lines:

        split_line = line.split()
        pdb, chain, first_pair_res, second_pair_res = (split_line[0] , split_line[6][:1], split_line[6][1:], split_line[7][1:])
        line_information = [chain, first_pair_res, second_pair_res]

        if pdb in pdb_dict:
            pdb_dict[pdb].append(line_information)

        else:
            pdb_dict[pdb] = [line_information]

    return(pdb_dict)

def kink_string_maker(pdb_dict, folder_directory, output_directory):
    """
    pdb_file ="1cza.pdb"
    residues ="'415-427'"
    chain = None

    kink_finder_str = "python2.7 ~/Downloads/KF_err_lin/Kink_Finder.py"
    pdb_file_location_str = "-f ~/Git/Project/refactored_code/kink_finder_tests/2angle/" + pdb_file
    kink_output_location_str = "-o ~/Git/Project/refactored_code/kink_finder_tests/2angle/1cza.pdb_output/"
    residues_display_str ="-l "+ residues +" -d"


    list_of_strings=[kink_finder_str, pdb_file_location_str, kink_output_location_str, residues_display_str]
    final_string = " ".join(list_of_strings)
    """
    list_of_strings = []

    for pdb_name in pdb_dict:
        residue_pairs = []

        for pdb_res in pdb_dict[pdb_name]:
            new_pair = pdb_res[1]+ "-" + pdb_res[2]
            if new_pair in residue_pairs:
                pass
            else:
                residue_pairs.append(new_pair)
        residue_pairs_string = (" ").join(residue_pairs)

        kink_finder_str = "python2.7 ~/Downloads/KF_err_lin/Kink_Finder.py"
        pdb_file_location_str = "-f "+ folder_directory + pdb_name +".pdb"
        kink_output_location_str = "-o " + output_directory + pdb_name + "_output/"
        residues_display_str = "-l '" + residue_pairs_string + "' -d"

        temp_string = [kink_finder_str, pdb_file_location_str, kink_output_location_str, residues_display_str]
        commandline_string = " ".join(temp_string)
        list_of_strings.append(commandline_string)


    return(list_of_strings)


def bash_command_process(list_of_strings,output_directory):
    for bash_command in list_of_strings:

        subprocess.run([bash_command], shell=True)
        angles_csv = output_directory +"angles.scv"
        print(angles_csv)
    return()

def function4():
    # do blargh
    return()

#angle_filename = "/home/oliver/Git/Project/refactored_code/kink_finder_tests/small_sample.2angle"
#folder_directory = "~/Git/Project/refactored_code/kink_finder_tests/"
#output_directory = folder_directory +"kink_output/"

angle_filename = "/home/oliver/Git/Project/helix_counts/normal_helices/normal_helix_03.08.2019_1angle.txt"
folder_directory = "/home/oliver/Documents/one_dir_all_pdbs/"
output_directory = folder_directory +"normal_helix_03.08.2019_1angle_kink_angles/"


pdb_dict=(pdb_res_pair(angle_filename))
list_of_commandline_strings = kink_string_maker(pdb_dict, folder_directory, output_directory)
for x in list_of_commandline_strings:
    print(x)

bash_command_process(list_of_commandline_strings,output_directory)
