#this is a test to create a method to find co occuring values 
#in a dictionary to display them on the tree only

dict = {'Node1': ['G432B', 'H454C', 'T345T'], 'Node2': ['C321A', 'G432Y', 'H454C']}
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
    print("this is the array in the dictionary")
    raw_input(dict[ks[i]])
    del dict[ks[i]]
    raw_input("dict value erased here: " )
    print dict
    #print values
    for j in range(len(values)):
        raw_input("This is the current array element im on: " +  values[j])
        for (key, items) in dict.items():
            #raw_input("this is items: " )
            #print items
            #raw_input("this is the key: " )
            #print key
            if values[j] in items:
                # neither key are in the co-occuring dict
                raw_input(multi_subs_dict)
                if key not in multi_subs_dict and ks[i] not in multi_subs_dict:
                    multi_subs_dict[key] = [values[j]]
                    multi_subs_dict[ks[i]] = [values[j]]
                # both keys are already in co-occuring dict
                elif key in multi_subs_dict and ks[i] in multi_subs_dict:
                    raw_input(multi_subs_dict[ks[i]])
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
                    multi_subs[ks[i]] = [values[j]] 
                    if values[j] not in multi_subs_dict[key]: #check if the value is already stored inside the co-occuring key value array 
                        multi_subs_dict[key].append(values[j])
            else:
                raw_input("You are now in the else statement!")
            raw_input(multi_subs_dict)
                
                
    dict[ks[i]] = values
 #   raw_input("dict value added back to dictionary vvv:\n " )
 #   print (dict) 
    
print multi_subs_dict    
