import sys
 
def species_control(args):
    values_list = args.values_list
    
    if values_list == None:
        if args.command == None:
            print("please provide a list of species you want to write the tests for. ")
            sys.exit(1)
        else:
            separate_cmd = args.command.split(" ")
            try:
                index_species = separate_cmd.index("-s")
                values_list = separate_cmd[index_species + 1]
            except:
                print("please provide a list of species you want to write the tests for. You can use the argument --values_list for that")
                sys.exit(1) 
    else:
        # checked in check_individual_command in CreateOutputs
        pass
                
    return values_list