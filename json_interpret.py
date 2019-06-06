import json
import sys
import re
#from Bio.Align.Applications import MuscleCommandline
import collections
from subprocess import *
import math
import os
#TODO: 
# modify script to display on tree the number 
# a substitution occurs x number of times  
#

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
    #print 'This is line 120' 
    #print hxb2_start
    #raw_input(type(hxb2_start))
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

#this function finds all the mutations that occur only once throughout the substituion tree
def single_mutant_dict_gen(full_list):
    
    sing_mut_dict = {}
    arr = []
    nodes = full_list.keys()
    for i in range(len(nodes)):
        vals = full_list[nodes[i]]
    #    print vals
        for  j in range(len(vals)):
            arr.append(vals[j])
    #print arr
    uniq = []
    temp = arr
    for i in range(len(arr)):
        #print(arr[i])
        #print i
        cur_mut = arr[i]
    #    raw_input(cur_mut)
        temp[i] = 'Flag'
    #    raw_input(temp)
        if cur_mut not in temp:
            uniq.append(cur_mut)
    #        print 'THere is a unique thing here'
        temp[i] = cur_mut
#    print temp
   # print uniq
    for i in range(len(uniq)):
#        print i
    #    raw_input(uniq[i])
        for j in range(len(nodes)):
#            print full_list[nodes[j]]
     #       raw_input()
     #       raw_input(uniq[i])
            if uniq[i] in full_list[nodes[j]]:
                if nodes[j] not in sing_mut_dict:
                    sing_mut_dict[nodes[j]] = [uniq[i]]
                else:
                    sing_mut_dict[nodes[j]].append(uniq[i])
    #print sing_mut_dict
    return sing_mut_dict

            


#this funtion finds the mutations that occur more than once throughout the substituion map multiple times
#TEST: PASSES ALL TEST
def co_occur_dict_gen(dict):
    multi_subs_dict = {}
    #print dict                                                                                                                                       
    #this is an array of keys                                                                                                                        
    ks = dict.keys()
    #print ks                                                                                                                                        
    #raw_input(dict.items())                                                                                                                         
    #loop over all the nodes in the dictionary                                                                                                      
    for i in range(len(ks)):
    #store the values of current key in an array                                                                                                     
        values = dict[ks[i]]
        #raw_input(values)                                                                                                                            
        #print("this is the array in the dictionary")
        #raw_input(dict[ks[i]])
        del dict[ks[i]]
        #raw_input("dict value erased here: " )
        #print dict
        #print values              
        for j in range(len(values)):
            #raw_input("This is the current array element im on: " +  values[j])
            for (key, items) in dict.items():
                if values[j] in items:
                    # neither key are in the co-occuring dict
                    if key not in multi_subs_dict and ks[i] not in multi_subs_dict:
                        multi_subs_dict[key] = [values[j]]
                        multi_subs_dict[ks[i]] = [values[j]]
                        # both keys are already in co-occuring dict                                                                                   
                    elif key in multi_subs_dict and ks[i] in multi_subs_dict:
                        #raw_input(multi_subs_dict[ks[i]])
                    #if the co occuring value is already in, then do nothing but
                #if it is not then add it to key/ks[i] array of values                                                                             
                        if values[j] not in multi_subs_dict[ks[i]]:
                            multi_subs_dict[ks[i]].append(values[j])
                        if values[j] not in multi_subs_dict[key]:
                            multi_subs_dict[key].append(values[j])
                #if one key is present and the other one isn't conditionals                                                                           
                    elif key not in multi_subs_dict and ks[i] in multi_subs_dict:
                        multi_subs_dict[key] = [values[j]]
                        if values[j] not in multi_subs_dict[ks[i]]: #chekck if the value                                                          
                            multi_subs_dict[ks[i]].append(values[j])
                    elif key in multi_subs_dict and ks[i] not in multi_subs_dict:
                        multi_subs_dict[ks[i]] = [values[j]]
                        if values[j] not in multi_subs_dict[key]: #check if the value is already stored inside the co-occuring key value array      
                            multi_subs_dict[key].append(values[j])
                else:
                    #raw_input("You are now in the else statement!")
                    pass
        #raw_input(multi_subs_dict)
        dict[ks[i]] = values
    #raw_input(multi_subs_dict)
    return multi_subs_dict

def find_changes(ref_codon, sub_codon):
    arr = []
    if ref_codon[0] != sub_codon[0]:
        arr.append(0)
    if ref_codon[1] != sub_codon[1]:
        arr.append(1)
    if ref_codon[2] != sub_codon[2]:
        arr.append(2)
    return arr

def non_synon_dict_gen(this_dict):
    keys = this_dict.keys()
    non_synon_nodes = {}
    for i in range(len(keys)):
        cur_key = keys[i]
        for j in range(len(this_dict[cur_key])):

            if this_dict[cur_key][j][0] != this_dict[cur_key][j][len(this_dict[cur_key][j]) - 2]:
                if cur_key not in non_synon_nodes:
                    non_synon_nodes[cur_key] = [this_dict[cur_key][j]]
#                    raw_input('added new shit to the new key ' + cur_key)
                else:
#                    raw_input(non_synon_nodes[cur_key])
#                    raw_input(non_synon_nodes)
                    non_synon_nodes[cur_key].append(this_dict[cur_key][j])
#                    raw_input('added new shit to the old key ' + cur_key)
#    raw_input(non_synon_nodes)
    return non_synon_nodes
#TODO: function to handle if there are insertions in the hxb2 site that are at sergei's reported site
if __name__ == "__main__":
    #os.system('module load muscle')
    
    dir = sys.argv[1].split('.')[0]
    #raw_input('this is dir: ' + str(dir) + 'this is dir\'s type' + str(type(dir)) )
    os.system('mkdir ' + dir)
    #locate_sequence(test_sequence)
    #pause = raw_input()
    #output file generation
    asr_fname = sys.argv[1].split('.')[0] + "_ancestral_seq.fas"
    outf_name = sys.argv[1].split('.')[0] + "multiple.tree"
    outf_single_name = sys.argv[1].split('.')[0] + "single.tree"
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
                #raw_input('This is the var nuc_number: ' + str(nuc_number))
                #print('This is the arr_hxb2_seq_range var: ')
                #raw_input(arr_hxb2_seq_range)
                hxb2_site = hxb2_amino_converter(nuc_number, arr_hxb2_seq_range[0], arr_hxb2_seq_range[5])
                if int(hxb2_site) == -10000:
                    errors_log.write(my_sub_nodes[num_nodes] +  ref + '\t\t\t' + sub + '\t\t\t' + "HXB2 site could not be found" + arr_hxb2_seq_range[2] + '\n')
                else:
                    subs_all_f.write(my_sub_nodes[num_nodes] + '\t' + ref + '\t\t\t' + sub + '\t\t\t' + hxb2_site + '\t\t\t' + arr_hxb2_seq_range[2] + '\n')
                    if my_sub_nodes[num_nodes] not in string_dict_nodes:
                        string_dict_nodes[my_sub_nodes[num_nodes]] = [ref + hxb2_site + sub + '_']
                    else:
                        string_dict_nodes[my_sub_nodes[num_nodes]].append(ref + hxb2_site + sub + '_')
                #raw_input('This is the final nucleotide site' + str(hxb2_site))
#                if hxb2_site != -10000:
                    #adds new node to global dictionary
#                    if my_sub_nodes[num_nodes] not in string_dict_nodes:
#                        string_dict_nodes[my_sub_nodes[num_nodes]] = [ref + hxb2_site + sub + ' ']
                    #adds characters to existing values
#                    else:
#                    print 'Inside else condition so that means that a value should be appended to an existing key!'
#                    print my_sub_nodes[num_nodes]

#                        string_dict_nodes[my_sub_nodes[num_nodes]].append(ref + hxb2_site + sub + ' ')
                    #raw_input (string_dict_nodes[my_sub_nodes[num_nodes]])
                    #string_dict_nodes[my_sub_nodes[num_nodes]].append(hxb2_site)
                    #string_dict_nodes[my_sub_nodes[num_nodes]].append(sub)
                    #string_dict_nodes[my_sub_nodes[num_nodes]].append('_')
                #write alignments to the files
                #print(string_dict_nodes)
                    aligned_asr_file.write('>' + my_sub_nodes[num_nodes] + '\n' + arr_hxb2_seq_range[4] + '\n')
                    hxb2_aligned_file.write('>' + my_sub_nodes[num_nodes] + " " +  arr_hxb2_seq_range[0] + '-' +  arr_hxb2_seq_range[1] + '\n' + arr_hxb2_seq_range[5] + '\n')
                
#                string_dict_nodes[my_sub_nodes[num_nodes]] = ''.join(string_dict_nodes[my_sub_nodes[num_nodes]])
                #print  string_dict_nodes[my_sub_nodes[num_nodes]]
                #pause = raw_input()
                 #print(string_dict_nodes) 
    #update nodes in newick string to proper amino acid format
    aligned_asr_file.close()
    subs_all_f.close()
    
############## UNCOMMENT THIS BLOCK IF SOMETHING GOES WRONG ##############
    #raw_input(nodes_list)
    #i = 0
    nodes_list.sort() #UNCOMMENT THIS LINE IF YOU NEED TO GO BACK
    #order the nodes 
    #ordered = collections.OrderedDict(sorted(string_dict_nodes.items()))
##########################################################################
    #turn tree string into a char array
    char_tree = list(tree)
    counter  = len(char_tree) - 1
####################### NEW BLOCK OF UPDATES #####################
    #get rid of synonymous mutations
    no_synon_mutant =  non_synon_dict_gen(string_dict_nodes)
    #raw_input(no_synon_mutant)
    singl_occur_dict = single_mutant_dict_gen(no_synon_mutant)
    multi_occur_dict = co_occur_dict_gen(no_synon_mutant) #call to 
    
 #   raw_input(singl_occur_dict)

    list_of_keys = multi_occur_dict.keys()
    list_of_sing_keys = singl_occur_dict.keys()
    #raw_input(multi_occur_dict)
    flag = 0 #this is the flag for the is match function
    for cur_node in nodes_list:
     
        #print cur_node
        #print "^^^ this is the current node"
        #print ("flag at the begininng of the loop..." + str(flag))
    #pause = raw_input()
        #if cur_node in string_dict_nodes: #if substitutions in the node were found through the map then create that string from the values 
        if cur_node in multi_occur_dict: 
            #newick_node = ''.join(string_dict_nodes[cur_node])
            newick_node = ''.join(multi_occur_dict[cur_node])
        #reads string in reverse to match current node from list of nodes
        #with the string in newick tree format
        
#############PICK UP HERE #######################
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
    #json_file.close()
    f.close()

    flag = 0
    counter = len(char_tree) - 1
    char_tree = list(tree)
    for cur_node in nodes_list:
        if cur_node in singl_occur_dict:
            newick_node = ''.join(singl_occur_dict[cur_node])
#############PICK UP HERE #######################                                                                                                                                                                                      
    
        #pause = raw_input()                                                                                                                                                                                                    
    
            while (counter > 4 and flag == 0):
            #print "im inside the while loop rn..."                                                                                                                                                                               
            #print("this is the value of flag..." + str(flag))                                                                                                                                                                        
                list_cur_node = list(cur_node)
            #print ('before the if statement that looks at node number')                                                                                                                                                         
                if char_tree[counter] == list_cur_node[len(list_cur_node) - 1]:

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
                #pause = raw_input()                                                                                                                                                                                
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
    with open(outf_single_name, 'w+') as o_single_f:                                                                                                       
        for char_item in char_tree:
            o_single_f.write('%s' % char_item )
    json_file.close()
    o_single_f.close()

    os.system('mv ' + asr_fname + ' ' + dir)
    os.system('mv ' + outf_name + ' ' + dir)
    os.system('mv ' + muscle_outfile + ' ' + dir)
    os.system('mv ' + subs_all_file + ' ' + dir)
    os.system('mv ' + hxb2_out_file + ' ' + dir)
    os.system('mv ' + outf_single_name + ' ' + dir)
