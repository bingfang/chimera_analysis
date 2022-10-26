#!/usr/local/bin/python3.6

import glob
import unittest

#save only "passed" and AF>0.2, DP>50 variants


def main():
    for name in glob.glob("../VQSR_filtered/snp.recalibrated.vcf"): #input is the vcf file 
        inputfile = str(name)
        outputfile = str(name)[:-17]+"_VQSR_AFDP_filtered.txt"
        print(outputfile )
        
        with open(inputfile,'r') as f:
            data_in = f.read().rstrip().split('\n') # Read file and create list split by newline 
        print(len(data_in))
        
        header_line_number, header = find_header(data_in)
        print("header_line_number:", len(header),header_line_number)
        print("Total input SNPs:", len(data_in)-len(header))
        passed=filter_PASS(data_in,header_line_number)
        
        AF_filtered=filter_AF(passed)
        high_quality_snp=append_AF(AF_filtered)

        
        print("High quality SNPs:",len(high_quality_snp))   
        
        
        
        somatic, both_num, mom_carrier_num, dad_carrier_num, both, mom_carrier, dad_carrier=find_chimera(high_quality_snp)
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
def append_AF(AF_filtered):
    high_quality_snp=[]    
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
        if "multiAD" not in line and "LowDP" not in line:
            high_quality_snp.append(line)      
    return high_quality_snp 



                  
# select line with Pblood AF >0.2   
def find_chimera(high_quality_snp):
    expected_somatic = 0
    somatic = 0
    both_num = 0
    both_invisible =0
    mom_carrier_num = 0
    mom_invisible = 0
    dad_carrier_num = 0
    dad_invisible = 0
    both=[]
    mom_carrier=[] 
    dad_carrier=[]

    for line in high_quality_snp:
        field=line.split('\t')

        mother_AF=float(field[14])
        father_AF=float(field[16])
        p_blood_AF=float(field[18])
        p_tumor_AF=float(field[20])
        brother_AF=float(field[22])
        
        if p_tumor_AF > 0 and p_blood_AF==0:
            expected_somatic += 1
        
        if p_tumor_AF > 0 and p_blood_AF==0 and mother_AF==0 and father_AF==0:
            somatic += 1

        if p_tumor_AF > 0 and p_blood_AF==0 and mother_AF>0 and father_AF>0:
            both_num += 1
            both.append(line+'\tboth carriers)')
        if p_tumor_AF == 0 and p_blood_AF==0 and mother_AF>0 and father_AF>0:
            both_invisible += 1
            
        if p_tumor_AF > 0 and p_blood_AF==0 and mother_AF>0 and father_AF==0:
            mom_carrier_num += 1
            mom_carrier.append(line+'\tmom carriers)')
        if p_tumor_AF == 0 and p_blood_AF==0 and mother_AF>0 and father_AF==0:
            mom_invisible += 1

        elif p_tumor_AF > 0 and p_blood_AF==0 and mother_AF==0 and father_AF>0:
            dad_carrier_num += 1
            dad_carrier.append(line+'\tdad carriers)')  
        if p_tumor_AF == 0 and p_blood_AF==0 and mother_AF==0 and father_AF>0:
            dad_invisible  += 1
 
            
    print("expected somatic variants: ", expected_somatic)
    print("somatic variants: ", somatic)
    print("unexpeceted variants: ", expected_somatic-somatic)
    
    print("Father is carrier, not in p-blood: ", dad_carrier_num + dad_invisible )
    print("in p-tumor, in father, not in p-blood and mother:", dad_carrier_num)
    print("not in p-tumor, in father, not in p-blood and mother:", dad_invisible)
    print ("chimeral percentage: ",dad_carrier_num/dad_invisible*2/(1+dad_carrier_num/dad_invisible)) 
    
    print("mother is carrier, not in p-blood: ", mom_carrier_num + mom_invisible )
    print("in p-tumor, in mother, not in p-blood and father:", mom_carrier_num)
    print("not in p-tumor, in mother, not in p-blood and father:", mom_invisible)
    print ("chimeral percentage: ",mom_carrier_num/mom_invisible*2/(1+mom_carrier_num/mom_invisible)) 
    
    print("both are carriers, not in p-blood: ", both_num + both_invisible )
    print("in p-tumor, in mother and father, not in p-blood :", both_num)
    print(" not in p-tumor, in mother and father, not in p-blood :", both_invisible)
    print ("chimeral percentage: ",both_num/both_invisible*4/(3+3*both_num/both_invisible)) 
    
    return somatic, both_num, mom_carrier_num, dad_carrier_num, both, mom_carrier, dad_carrier
    
    
class TestFilter_append(unittest.TestCase):
        
    def test_find_header(self):
        dat=['##source=ApplyVQSR','##source=CombineGVCFs','#CHROM']
        header_line_number, header = find_header(dat)  
        self.assertEqual(len(header), header_line_number)

main()









