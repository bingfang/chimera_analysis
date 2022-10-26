#!/usr/local/bin/python3.6

import glob

#save only "passed" and AF>0.2, DP>50 variants


def main():
    for name in glob.glob("./snp.recalibrated.vcf"): #input is the vcf file 
        inputfile = str(name)
        outputfile = str(name)[:-4]+"_inFather_mix2_3.txt"
        with open(inputfile,'r') as f:
            data_in = f.read().rstrip().split('\n') # Read file and create list split by newline 
        print(len(data_in))
        
        header_line_number, header = find_header(data_in)
        print("header_line_number:", len(header))
        print("Total input variants:", len(data_in)-len(header))
        passed=filter_PASS(data_in,header_line_number)
        print("Total input SNPs:", len(passed))

        
        AF_filtered=filter_AF(passed)
        print("High quality SNPs:",len(AF_filtered))   
        
        mix_AF=find_AF(AF_filtered)
        
        mix1,mix2,mix3,mix4,mix5, mix6, mix7,mix8,mix9, mix10=find_mix(mix_AF)
        with open(outputfile,'w') as f:
            for lst in [mix2,mix3]:
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
    print(len(header)) 
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

# append AF to the end
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
                AF=""
            AF_list.append(str(AF) + '\t'+DP)
        line=line+'\t'+('\t').join(AF_list)
        mix_AF.append(line)      
    return mix_AF 
       
# select line with Pblood AF >0.2   
def find_mix(mix_AF):
    mix1=[]
    mix2=[]
    mix3=[] 
    mix4=[]
    mix5=[] 
    mix6=[]
    mix7=[]
    mix8=[] 
    mix9=[]
    mix10=[]
    uncommon=0
    p_tumor_mix = 0
    for line in mix_AF:
        field=line.split('\t')
        mother=field[9]
        father=field[10]
        p_blood=field[11]
        p_tumor=field[12]
        brother=field[13]
        if "0/0" not in p_tumor and "0/0" in p_blood:
            p_tumor_mix += 1
            if "0/1" in p_tumor:
                mix1.append(line+'\t1')
            else:
                uncommon += 1 
        if "0/1" in p_tumor and "0/0" in p_blood and "0/0" in mother and "0/1" in father:
            mix2.append(line+'\t1\t2')
        if "0/1" in p_tumor and "0/0" in p_blood and "0/0" in father and "0/1" in mother:
            mix3.append(line+'\t1\t3')
        if "0/1" in p_tumor and "0/0" in p_blood and "0/1" in brother:
            if "0/1" in mother or "0/1" in father:
                mix4.append(line+'\t1\t4')
        if "0/1" in p_tumor and "0/0" in p_blood:
            if "0/1" in mother and "0/1" in father: 
                mix5.append(line+'\t5')
        if "0/1" in p_tumor and "0/0" in p_blood:
            if "0/1" in mother or "0/1"  in father :
                mix6.append(line+'\t6')                        
        if "0/0" in p_blood and "0/0" in mother and "0/0" in father and "0/1" in p_tumor and "0/0" in brother :
            mix7.append(line+'\t7')
            
        if "0/0" in p_blood and "0/1" in p_tumor and "0/1" in mother and "0/1" in father:
            mix8.append(line+'\t8') 
        if "0/0" in p_blood and "0/0" in p_tumor and "0/1" in mother and "0/1" in father:
            mix9.append(line+'\t9') 

        if "0/0" in p_blood and "0/1" in mother and "0/1" in father:
            mix10.append(line+'\t10')     
    print("in p-tumor, not in p-blood: total %d uncommon %d common %d %f2" % (p_tumor_mix, uncommon, len(mix1), len(mix1)/5272))
    print("in p-tumor, in father, not in p-blood and mother:", len(mix2),len(mix2)/5272)
    print("in p-tumor, in mother, not in p-blood and father:", len(mix3),len(mix3)/5272)
    print("in p-tumor, in brother, not in p-blood:",len(mix4),len(mix4)/5272)  
    print("in p-tumor, in both parents, not in p-blood:",len(mix5),len(mix5)/5272) 
    print("in p-tumor, not in p-blood, in one of parents:",len(mix6),len(mix6)/5272)

    print("in p-tumor, not in p-blood and others (somatic):", len(mix7), len(mix7)/5272)
    print("\n")
    print("mother is carrier, not in p-blood, in p-tumor:", len(mix8),len(mix8)/5272)
    print("mother is carrier, not in p-blood, not in p-tumor:", len(mix9),len(mix9)/5272)
    print("father and mother are carriers, not in p-blood:", len(mix10),len(mix10)/5272)  
    return mix1,mix2,mix3,mix4,mix5,mix6,mix7,mix8, mix9,mix10

main()
