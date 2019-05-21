import json
import sys
import re
from Bio.Align.Applications import MuscleCommandline
import collections
from subprocess import *
import math
import os

#function calls ruby script in local environment and returns
# boolean (if there are gaps in alignment), the aligned sequence
# the hxb2 sequence, and the hxb2 position
#PASSES ALL TESTS!
def locate_sequence(sequence):
    slave = Popen(['ruby', 'sequence_hxb2.rb', sequence], stdin=PIPE, stdout=PIPE, stderr=STDOUT)
    stdout_ = slave.communicate()[0]
    out_data = stdout_.split()
    return out_data

#aligns sequence fasta file and generates new out file
#PASSES ALL TESTS!
#TODO: put path of where muscle is
def muscle_aligner_func(in_m_file, out_m_file):
    muscle_exe = r""
    muscle_cline = MuscleCommandline(muscle_exe, input=in_m_file, out=out_m_file)
    #print(muscle_cline)
    stdout, stderr = muscle_cline()

#converts codon to its amino acid pair
#input: Codon
#output: Amino acid
#test: codon_converter("ATA") 
#PASSES ALL TESTS!
def codon_converter(codon):
    
    dna_amino_keys = {('ATT', 'ATC', 'ATA'): 'I', ('CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'): 'L', ('GTT', 'GTC', 'GTA', 'GTG'): 'V', ('TTT', 'TTC'): 'F',
                      'ATG': 'M', ('TGT', 'TGC'): 'C', ('GCT', 'GCC', 'GCA', 'GCG'): 'A', ('GGT', 'GGC', 'GGA', 'GGG'): 'G', ('CCT', 'CCC', 'CCA', 'CCG'): 'P',
                      ('ACT', 'ACC', 'ACA', 'ACG'): 'T', ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'): 'S', ('TAT', 'TAC'): 'Y', 'TGG': 'W', ('CAA', 'CAG'): 'Q',
                      ('AAT', 'AAC'): 'N', ('CAT', 'CAC'):'H', ('GAA', 'GAG'): 'E', ('GAT', 'GAC'):'D', ('AAA', 'AAG'): 'K', ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'): 'R'}
    
    return next(v for k, v in dna_amino_keys.items() if codon in k)

#identifies if arrays have the same values in the same order
#input: arr1 and arr2
#output: true or false
#test: is_match(["N", "o","d", "e", "8"], ["N", "o", "d", "e", "83"]) 
#PASSES ALL TESTS!
def is_match(new_char, list_node_arr):
    for i in range(len(new_char)):
        if new_char[i] != list_node_arr[i]:
            #print "the arrays are NOT the same!"
            return 0
    #print "the arrays ARE the same!"
    #pause = raw_input()
    return 1

#this function advances the indicy to a position where there as a nucleotide
#this is not used 
def advance_nuc_marker(q_sequence, indx):
    #print "Inside advance nuc marker function !"
    while(q_sequence[indx] == '-'):
        #print "index = " + str(indx)
        indx = indx + 1
    return indx

#identifies the sequence from the unaligned asr file 
#input: a node and a string that holds the sequences
#output: string with the string
#test:PASSES ALL TESTS!
# asr_sequence_finder("Node33", ">Node3\nAAT\n>Node3\nATATEEGS\n")
def asr_sequence_finder(cur_node, entire_doc):
    #print "this is the entire_doc variable"
    #print entire_doc
    #pause  = raw_input()
    ed_arr = entire_doc.split()
    #print "the is the ed_arr variable..."
    #print ed_arr
    #pause = raw_input()
    if '>' + cur_node in ed_arr: #finds the node that was used a parameter for this function 
        indx = ed_arr.index('>' + cur_node) #get that index 
    #print indx
    #pause = raw_input()
    return str(ed_arr[indx + 1]) #returing the sequence

#this function gets the nucleotide(s) positions that were
#mutated 
#parameters: sergei's site and a ancestral sequence 
#test: PASSES ALL TESTS!
def find_nucleotide(sergei_indx, asr_sequence):
    #print "this is sergei's site: " + str(sergei_indx)
    #pause = raw_input()
    no_int_hxb2 = ""
    
    n = 3
    asr_codons = [asr_sequence[i:i+n] for i in range(0, len(asr_sequence), n)]
    #asr_codons = [asr_sequence[i:i+n] for i in range(0, len(asr_sequence), n)]
    #raw_input("This is the codon from the asr sequence from sergei's indx!" + str(asr_codons[sergei_indx]))
    
    new_codon = "XXX"    
    #pause = raw_input("this is the new codon! " + new_codon)
    asr_codons[sergei_indx] = new_codon
    asr_seq_ = ''.join(asr_codons)
    #print "This is the asr_seq_ with the XXX's: " + asr_seq_ 
    #j = len(asr_seq_)
    nuc_numbers = 0
    #if len(nuc_diff) == 1:
    for i in range(len(asr_seq_)):
        if asr_seq_[i] == 'X':
            nuc_numbers = i 
            break
    return nuc_numbers

#this function returns the hxb2 number associaited with a nucleotide position 
#parameters: nucleotide position, the starting point where the sequence aligns
# to in the envelope, and the hxb2 sequence
#tests: PASSES ALL TESTS
#hxb2_amino_converter(4, 6560, "ATA--ATA")
def hxb2_amino_converter(nuc_num, hxb2_start, hxb2_sequence):
    counter = 0
    #print "this is the nucleotide number" + str(nuc_num)
    #print('The is the len of hxb2_sequence: ' + str(len(hxb2_sequence)))
    #raw_input('This is the hxb2_sequence: ' + str(hxb2_sequence)) 
    #raw_input('This is the hxb2_start before the while loop: ' + str(hxb2_start))
    indx = 0 #once this equals the nuc number than the nucleotide in the hxb2_sequence has been found
    hxb2_num = int(hxb2_start) - 1 #first index will be the begining of where the sequence aligns to  
    while int(nuc_num) > indx: #loop until indexer is greater than zero

        if hxb2_sequence[counter] != '-': #only increment counter if the the char is not a hyphen
            hxb2_num = hxb2_num + 1
            indx = indx + 1
            #raw_input('Index incrememnted to: ' + str(indx) + 'this is the nucleotide: ' + str(hxb2_sequence[counter]) + 'this site hxb2 site: ' + str(hxb2_num))
        counter = counter + 1
        if counter == len(hxb2_sequence) - 1:
            return str(-10000)
        #print  "this is the counter" + str(counter)
    #raw_input('This is the hxb2_converter' + str(hxb2_num))
    hxb2_nuc_site = hxb2_num - 6224
    #raw_input('This is the hxb2_nuc_site/3.0:' + str(hxb2_nuc_site/3.0))
    hxb2_amino_site = int(math.ceil(hxb2_nuc_site/3.0))
    return str(hxb2_amino_site)

def find_changes(ref_codon, sub_codon):
    arr = []
    if ref_codon[0] != sub_codon[0]:
        arr.append(0)
    if ref_codon[1] != sub_codon[1]:
        arr.append(1)
    if ref_codon[2] != sub_codon[2]:
        arr.append(2)
    return arr


#TODO: function to handle if there are insertions in the hxb2 site that are at sergei's reported site
if __name__ == "__main__":
    dir = sys.argv[1].split('.')[0]
    #raw_input('this is dir: ' + str(dir) + 'this is dir\'s type' + str(type(dir)) )
    os.system('mkdir ' + dir)
    #locate_sequence(test_sequence)
    #pause = raw_input()
    #output file generation
    asr_fname = sys.argv[1].split('.')[0] + "_ancestral_seq.fas"
    outf_name = sys.argv[1].split('.')[0] + ".tree"
    muscle_outfile = sys.argv[1].split('.')[0] + "_muscle_aligned.fas"
    subs_all_file = sys.argv[1].split('.')[0] + "_substitutions.txt"
    hxb2_out_file = sys.argv[1].split('.')[0] + "_hxb2_sequences.fas"
    errors = sys.argv[1].split('.')[0] + "_errors_log"
    #get the json output data from the hyphy run
    with open(sys.argv[1]) as json_file:
        json_data = json.load(json_file)

    #gets the ancestral sequences from the json file and writes them to **_ancestral_seq.fas 
    unaligned_asr_file_lines = ""
    with open(asr_fname, 'w+') as asr_file:
        for node in json_data['ancestral_sequences']:
            asr_file.write('>' + node + '\n')
            asr_file.write(json_data['ancestral_sequences'][node] + '\n')
            unaligned_asr_file_lines = unaligned_asr_file_lines + str('>' + node + '\n' + json_data['ancestral_sequences'][node] + '\n')
    asr_file.close()
    
    #errors log file
    errors_log = open(errors, 'w+') 
    #this file is for the alignment lines that will be created through calls to the ruby program
    aligned_asr_file = open(muscle_outfile, 'w+') 
    #put the header on the substition file
    subs_all_f = open(subs_all_file, 'w+')
    subs_all_f.write("Node" + '\t' + "Amino Acid Reference" + '\t' + "Amino Acid Change" + '\t' +  "HXB2 Site Number" + '\t' + "% sequence match" + '\n')
    #this file is the for the hxb2 alignment sequences
    hxb2_aligned_file = open(hxb2_out_file, "w+")
    
    #dictionary holds all the substitutions for each node
    #keys = node, values = substitutions
    string_dict_nodes = {}

    #newick format tree data from json
    tree = str(json_data['tree'])
    #gets instances of nodes and adds to array
    nodes_list = re.findall(r'(Node\w+)', tree)
    node_count = len(nodes_list)
    #print nodes_list
    #pause = raw_input()

    #map of substituions 
    subs = json_data['substitution_map']

    node_txt = "Node"
  
    for site in subs:  #loop through map to find the substituions at nodes
        #print site
        #pause = raw_input()
        search_node = str(subs[site].values())
        #print subs[site]
        #pause = raw_input()
    
        #if there is a node in the substitution map at site X go
        #through the map locating the nodes, amino acid changes and HXB2 sites
        if re.search(node_txt, search_node):
            #print re.search(node_txt, search_node)
            #print type(re.search(node_txt, search_node))
            #pause = raw_input()
            #print("Node Found!")
            #print( str(subs[site].keys()[0]) + str(site)  + str(subs[site].values()[0].keys()[0]))
            #pause = raw_input()
            my_sub_nodes = re.findall(r'(Node\w+)', search_node)
            #loop through the number of nodes at this substitution site and get the amino acid substitution 
            #the reference amino acid and convert the site to hxb2 numbers
            for num_nodes in range(len(my_sub_nodes)):
                #"reference" amino acid
                ref = codon_converter(str(subs[site].keys()[0]))
                #amino acid substitution                                                                                                         
                sub = codon_converter(str(subs[site].values()[0].keys()[0]))
                print "ASR Substiution Site: " + str(site) + " Reference Codon Sequence: " + subs[site].keys()[0] +  " Substituion Codon Sequence: " + subs[site].values()[0].keys()[0]
                #pause = raw_input('---------------------------------------')
                #print "This is the node that im on!VVV"
                #print my_sub_nodes[num_nodes]
                #pause = raw_input()
                #find sequence from this node and pass it to the sequence locator for hxb2 sequence processing
                node_seq = asr_sequence_finder(my_sub_nodes[num_nodes],  unaligned_asr_file_lines)
                arr_hxb2_seq_range = locate_sequence(node_seq) 
                #find the nucleotide associated with this position and pass it to the hxb2 number finder
                #print arr_hxb2_seq_range[5]
                nuc_number = find_nucleotide(int(site), node_seq)
                hxb2_site = hxb2_amino_converter(nuc_number, arr_hxb2_seq_range[0], arr_hxb2_seq_range[5])
                if int(hxb2_site) == -10000:
                    errors_log.write(my_sub_nodes[num_nodes] +  ref + '\t\t\t' + sub + '\t\t\t' + "HXB2 site could not be found" + arr_hxb2_seq_range[2] + '\n')
                else:
                    subs_all_f.write(my_sub_nodes[num_nodes] + '\t' + ref + '\t\t\t' + sub + '\t\t\t' + hxb2_site + '\t\t\t' + arr_hxb2_seq_range[2] + '\n'\
)

                #raw_input('This is the final nucleotide site' + str(hxb2_site))
                #adds new node to global dictionary
                if my_sub_nodes[num_nodes] not in string_dict_nodes:
                    string_dict_nodes[my_sub_nodes[num_nodes]] = [ref] + [hxb2_site] + [sub] + ['_']
                #adds characters to existing values
                else:
                    string_dict_nodes[my_sub_nodes[num_nodes]].append(ref)
                    string_dict_nodes[my_sub_nodes[num_nodes]].append(hxb2_site)
                    string_dict_nodes[my_sub_nodes[num_nodes]].append(sub)
                    string_dict_nodes[my_sub_nodes[num_nodes]].append('_')
                #write alignments to the files
                aligned_asr_file.write('>' + my_sub_nodes[num_nodes] + '\n' + arr_hxb2_seq_range[4] + '\n')
                hxb2_aligned_file.write('>' + my_sub_nodes[num_nodes] + " " +  arr_hxb2_seq_range[0] + '-' +  arr_hxb2_seq_range[1] + '\n' + arr_hxb2_seq_range[5] + '\n')
                
#                string_dict_nodes[my_sub_nodes[num_nodes]] = ''.join(string_dict_nodes[my_sub_nodes[num_nodes]])
                #print  string_dict_nodes[my_sub_nodes[num_nodes]]
                #pause = raw_input()
                 #print(string_dict_nodes) 
    #update nodes in newick string to proper amino acid format
    aligned_asr_file.close()
    subs_all_f.close()
    i = 0
    nodes_list.sort()
    #order the nodes 
    ordered = collections.OrderedDict(sorted(string_dict_nodes.items()))

    #turn tree string into a char array
    char_tree = list(tree)
    counter  = len(char_tree) - 1

    flag = 0 #this is the flag for the is match function
    for cur_node in nodes_list:
        #print cur_node
        #print "^^^ this is the current node"
        #print ("flag at the begininng of the loop..." + str(flag))
    #pause = raw_input()
        if cur_node in string_dict_nodes: #if substitutions in the node were found through the map then create that string from the values 

            newick_node = ''.join(string_dict_nodes[cur_node])
        
        #reads string in reverse to match current node from list of nodes
        #with the string in newick tree format
        
        #print newick_node
        #pause = raw_input()
            while (counter > 4 and flag == 0):
            #print "im inside the while loop rn..."
            #print("this is the value of flag..." + str(flag))
                list_cur_node = list(cur_node)
            #print ('before the if statement that looks at node number')
            #pause  = raw_input()
                if char_tree[counter] == list_cur_node[len(list_cur_node) - 1]:
                
                #print ("this is char. from the tree: " + char_tree[counter] + " this is the char from the node list: " + list_cur_node[len(list_cur_node) - 1])
                #print char_tree[(counter - (len(list_cur_node) -1)): (counter + 1)]
                    #print ("match not found!")
                    #if the nodes are the same then switch out the strings
                    if is_match(char_tree[(counter - (len(list_cur_node) -1)): (counter + 1)], list_cur_node) and char_tree[counter + 1] == ':':
#                    print char_tree[counter+ 1]
#                    print("match has been found")
#                    pause = raw_input()
                        flag = 1
                        char_tree[(counter - (len(list_cur_node) -1))] = newick_node
                    #char_tree[(counter - (len(list_cur_node))): (counter + 1)]
                        char_tree[(counter + 2) - (len(list_cur_node)): (counter + 1)] = ''
#                    print char_tree[(counter - (len(list_cur_node) -1)): (counter + 1)]
#                    print ("flag is true..." + str(flag))
#                    pause = raw_input()
                #print char_tree[(counter - (len(list_cur_node) -1)): (counter + 1)] 
                #pause = raw_input()
                counter = counter - 1
                #print "this is the counter..." + str(counter)
                #print "this is the value of flag..." + str(flag)

            #pause = raw_input()
            #print "this is the END of the finding the node in the newick tree..."
            flag = 0
            counter = len(char_tree) - 1

    #write the new newick string to the file
    with open(outf_name, 'w+') as f:                                                                                                       
        for char_item in char_tree:
            f.write('%s' % char_item )
    json_file.close()
    f.close()

    os.system('mv ' + asr_fname + ' ' + dir)
    os.system('mv ' + outf_name + ' ' + dir)
    os.system('mv ' + muscle_outfile + ' ' + dir)
    os.system('mv ' + subs_all_file + ' ' + dir)
    os.system('mv ' + hxb2_out_file + ' ' + dir)
