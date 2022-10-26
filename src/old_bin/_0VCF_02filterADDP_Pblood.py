#!/usr/local/bin/python3.6

import glob

#save only "passed" and AF>0.2, DP>50 variants


def main():
    for name in glob.glob("./raw_variants.vcf"): #input is the vcf file 
        inputfile = str(name)
        outputfile = str(name)[:-4]+"_HQ_10K.txt"
        with open(inputfile,'r') as f:
            data_in = f.read().rstrip().split('\n') # Read file and create list split by newline 
        print(len(data_in))
        
        header_line_number, header = find_header(data_in)
        print("header_line_number:", len(header))
        print("Total input SNPs:", len(data_in)-len(header))
        #passed=filter_PASS(data_in,header_line_number, passed)
        
        AF_filtered=filter_AF(data_in,header_line_number)
        print("High quality SNPs:",len(AF_filtered))   
        
        mix_AF=find_AF(AF_filtered)
        
        mix1=find_mix(mix_AF)
        #print=("Pblood AF> 0.2: ", len(mix1))
    
        with open(outputfile,'w') as f:
            for line in header:
                f.write(str(line)+'\n')
            for line in mix_AF:
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
def filter_AF(data_in,header_line_number):
    AF_filtered=[]    
    for line in data_in[(header_line_number + 1):(len(data_in))]:
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
            AF_list.append(str(AF))
        line=line+'\t'+('\t').join(AF_list)
        mix_AF.append(line)      
    return mix_AF 
       
# select line with Pblood AF >0.2   
def find_mix(mix_AF):
    mix1=[]
    multiAD=0
    noAF = 0
    for line in mix_AF:
        field=line.split('\t')
        mother=field[9]
        father=field[10]
        p_blood=field[11]
        p_tumor=field[12]
        brother=field[13]
        pblood_AF=field[16]
        if "0/0" not in p_blood:
            if pblood_AF=="multiAD":
                multiAD += 1
            elif pblood_AF=="":
                noAF += 1
            elif float(pblood_AF)>0.2:
                mix1.append("\t".join(field[0:14]))      
    print("in p-blood:", len(mix1),len(mix1)/102482)
    print("Multi Alleles in p-blood:", multiAD)
    print("No AF in p-blood:", noAF)      
    return mix1

main()
