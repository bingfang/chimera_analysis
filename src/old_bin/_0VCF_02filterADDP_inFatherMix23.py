#!/usr/local/bin/python3.6

import glob

#save only "passed" and AF>0.2, DP>50 variants


def main():
    for name in glob.glob("./snp.recalibrated.vcf"): #input is the vcf file 
        inputfile = str(name)
        outputfile = str(name)[:-17]+"_VQSR_AFDP_filtered.txt"
        print(inputfile)
        
        with open(inputfile,'r') as f:
            data_in = f.read().rstrip().split('\n') # Read file and create list split by newline 
        print(len(data_in))
        
        header_line_number, header = find_header(data_in)
        print("header_line_number:", len(header))
        print("Total input SNPs:", len(data_in)-len(header))
        passed=filter_PASS(data_in,header_line_number)
        
        AF_filtered=filter_AF(passed)
        print("High quality SNPs:",len(AF_filtered))   
        
        mix_AF=find_AF(AF_filtered)
        
        somatic, low_DP, both_num, mom_carrier_num, dad_carrier_num, both, mom_carrier, dad_carrier=find_mix(mix_AF)
        with open(outputfile,'w') as f:
            for lst in [both, mom_carrier, dad_carrier]:
                for line in lst:
                    f.write(str(line)+'\n')

# append header lines to the header list.    
def find_header(data_in):
    header=[] 
    for i in range(len(data_in)):   
        if '#CHROM' not in data_in[i][0:6]:
            header.append(data_in[i])
        else:
            header.append(data_in[i])
            header_line_number = i
            break 
    return header_line_number, header

# append lines with quality "PASS" to the filtered list.   
def filter_PASS(data_in,header_line_number):
    passed=[]    
    for line in data_in[(header_line_number + 1):(len(data_in))]:
        field=line.split('\t')
        if field[6]=="PASS":
            passed.append(line)
    return passed
    
# append lines with AF>20%, DP>50 to the filtered list.    
def filter_AF(passed):
    AF_filtered=[]    
    for line in passed:
        Max_AF_allele = 0
        field=line.split('\t')
        for sample in field[9:len(field)]:
            reads=sample.split(':')
            AD=reads[1]
            DP=reads[2]
            if DP != "." and float(DP) > 50:
                AD_allele=AD.split(",")
                for j in range(1,len(AD_allele)): # first AD_allele is WT
                    AF_allele=float(AD_allele[j])/float(DP)
                    if float(AF_allele) > Max_AF_allele:
                        Max_AF_allele = float(AF_allele)
        if float(Max_AF_allele) > 0.20:
            AF_filtered.append(line)
     
    return AF_filtered

# append AF and DP to the end of line
def find_AF(AF_filtered):
    mix_AF=[]    
    for line in AF_filtered:
        AF_list=[]
        field=line.split('\t')
        for sample in field[9:len(field)]:
            reads=sample.split(':')
            AD=reads[1]
            DP=reads[2]
            AD_allele=AD.split(",")
            if len(AD_allele) == 2 and "." not in DP:
                if float(DP) != 0:
                    AF=float(AD_allele[1])/float(DP)
            elif len(AD_allele) > 2:
                AF ="multiAD"
            else:
                AF="LowDP"
            AF_list.append(str(AF) + '\t'+DP)
        line=line+'\t'+('\t').join(AF_list)
        mix_AF.append(line)      
    return mix_AF 
       
# select line with Pblood AF >0.2   
def find_mix(mix_AF):
    somatic = 0
    low_DP = 0
    both_num = 0
    mom_carrier_num = 0
    dad_carrier_num = 0
    both=[]
    mom_carrier=[] 
    dad_carrier=[]

    for line in mix_AF:
        field=line.split('\t')
        for elem in field[14:]
            try:
                float(element)
            except ValueError:
                low_DP += 1
        mother_AF=field[14]
        father_AF=field[16]
        p_blood_AF=field[18]
        p_tumor_AF=field[20]
        brother_AF=field[22]
        
        
        if p_tumor_AF > 0 and p_blood_AF==0 and mother_AF==0 and father_AF==0:
            somatic += 1
        elif p_tumor_AF > 0 and p_blood_AF==0 and mother_AF>0 and father_AF>0:
            both_num += 1
            both.append(line+'\tboth carriers)')
        elif p_tumor_AF > 0 and p_blood_AF==0 and mother_AF>0 and father_AF==0:
            mom_carrier_num += 1
            mom_carrier.append(line+'\tmom carriers)')
        elif p_tumor_AF > 0 and p_blood_AF==0 and mother_AF==0 and father_AF>0:
            dad_carrier_num += 1
            dad_carrier.append(line+'\tdad carriers)')    
            
    print("low DP ", low_DP)  
    print("somatic ", somatic)
    print("in p-tumor, in father, not in p-blood and mother:", dad_carrier_num)
    print("in p-tumor, in mother, not in p-blood and father:", mom_carrier_num)
    print("in p-tumor, in mother and father, not in p-blood :", both_num)
    return somatic, low_DP, both_num, mom_carrier_num, dad_carrier_num, both, mom_carrier, dad_carrier

main()
