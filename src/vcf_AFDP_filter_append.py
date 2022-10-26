#!/usr/local/bin/python3.8

import filter_append as fa
import glob

def main():
    for name in glob.glob("../data/VQSR_filtered/snp.recalibrated.vcf"): #input is the vcf file 
        inputfile = str(name)
        outputfile = str(name)[:-17]+"_VQSR_AFDP_filtered_dad_carrier_invisible.txt"
        print(outputfile )
        
        with open(inputfile,'r') as f:
            data_in = f.read().rstrip().split('\n') # Read file and create list split by newline 

        filtered=fa.FilterAppend(data_in) # build a new instance
        header_line_number, header = filtered.find_header()
        print("header_line_number:", len(header),header_line_number)
        print("Total input SNPs:", len(data_in)-len(header))
        passed=filtered.filter_PASS(header_line_number)
        
        
        filtered1=fa.FilterAppend(passed) # build a new instance
        AF_filtered=filtered1.filter_AF()

        
        filtered2=fa.FilterAppend(AF_filtered) # build a new instance
        high_quality_snp=filtered2.append_AF()
        print("High quality SNPs:",len(high_quality_snp)) 
        
          
        filtered3=fa.FilterAppend(high_quality_snp) # build a new instance
        somatic, both_num, mom_carrier_num, dad_carrier_num, dad_invisible_num, both, mom_carrier, dad_carrier,dad_invisible=filtered3.find_chimera()
        with open(outputfile,'w') as f:
            for lst in [dad_carrier,dad_invisible]:
                for line in lst:
                    f.write(str(line)+'\n')
                    
main()