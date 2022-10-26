#!/usr/local/bin/python3.8

import filter_append as fa
import glob

def main():
    for name in glob.glob("../data/Chimera_sites/*.txt"): #input is the vcf file 
        inputfile = str(name)
        outputfile = str(name)[:-4]+"_remove_artifacts.txt"
        print(outputfile )
        
        with open(inputfile,'r') as f:
            data_in = f.read().rstrip().split('\n') # Read file and create list split by newline 
        
        # remove indel
        SNVs=filter_SNV(data_in)
        
        # remove homozygous SNVs from parents, remove pTumor DP<30
        both_carrier_count, mom_carrier_count, dad_carrier_count=filter_AFDP(SNVs)
        print(len(both_carrier_count), len(mom_carrier_count), len(dad_carrier_count))
       
        with open(outputfile,'w') as f:
            for lst in [both_carrier_count, mom_carrier_count, dad_carrier_count]:
                for line in lst:
                    f.write(str(line)+'\n')

def filter_SNV(data_in):
    SNVs=[]
    removed_num=0
    for line in data_in:
        field=line.split('\t')
        if len(str(field[3]))==1 and len(str(field[4]))==1 and "." not in field[2]:
            SNVs.append(line)
        else:
            removed_num +=1 
    print("number of indel: ", removed_num)  
    return SNVs 



def filter_AFDP(SNVs):
    both_carrier_count=[]
    mom_carrier_count=[]
    dad_carrier_count=[]
    remove_num=0
    for line in SNVs:
        field=line.split('\t')
        if "dad" in field[-1]:
            if float(field[21])>30 and float(field[16])>0.3 and float(field[16])<0.7:
                dad_carrier_count.append(line)
            else:
                remove_num +=1
        elif "mom" in field[-1]:
            if float(field[21])>30 and float(field[14])>0.25 and float(field[14])<0.75:
                mom_carrier_count.append(line)
            else:
                remove_num +=1
        elif "both" in field[-1]:
            if float(field[21])>30 and float(field[14])>0.25 and float(field[14])<0.75 and float(field[16])>0.25 and float(field[16])<0.75:
                both_carrier_count.append(line)
            else:
                remove_num +=1
        else:
            print("field[-1]")
    print("likely aretifacts: ", remove_num)
    return both_carrier_count, mom_carrier_count, dad_carrier_count
                


                    
main()