#!/usr/local/bin/python3.6

import glob

#FilterAppend handle vcf files. It filters variants based on AF and DP and appends AF and DP to the end of variant line.
class FilterAppend:


    def __init__(self,vcf_dat):  #vcp_dat is a list from vcf file splited by line
        self.vcf_dat = vcf_dat

    # append header lines to the header list.    
    def find_header(self):
        header=[] 
        for i in range(len(self.vcf_dat)):   
            if '#CHROM' not in self.vcf_dat[i][0:6]:
                header.append(self.vcf_dat[i])
            else:
                header.append(self.vcf_dat[i])
                header_line_number = i
                break 
        return header_line_number, header

    # append lines with quality "PASS" to the filtered list.   
    def filter_PASS(self,header_line_number):
        passed=[]    
        for line in self.vcf_dat[(header_line_number + 1):(len(self.vcf_dat))]:
            field=line.split('\t')
            if field[6]=="PASS":
                passed.append(line)
        return passed
    
    # append lines with AF>20%, DP>50 to the filtered list.    
    def filter_AF(self):
        AF_filtered=[]    
        for line in self.vcf_dat:
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
    def append_AF(self):
        high_quality_snp=[]    
        for line in self.vcf_dat:
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


    # classify and count the variants  
    def find_chimera(self):
        expected_somatic = 0
        somatic = 0
        both_num = 0
        both_invisible =0
        mom_carrier_num = 0
        mom_invisible = 0
        dad_carrier_num = 0
        dad_invisible_num = 0
        both=[]
        mom_carrier=[] 
        dad_carrier=[]
        dad_invisible = []
        for line in self.vcf_dat:
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

            if p_tumor_AF > 0 and p_blood_AF==0 and mother_AF==0 and father_AF>0:
                dad_carrier_num += 1
                dad_carrier.append(line+'\tdad carriers)')  
            if p_tumor_AF == 0 and p_blood_AF==0 and mother_AF==0 and father_AF>0:
                dad_invisible_num  += 1
                dad_invisible.append(line+'\tdad invisible') 
 
            
        print("expected somatic variants: ", expected_somatic)
        print("somatic variants: ", somatic)
        print("unexpeceted variants: ", expected_somatic-somatic)
    
        print("Father is carrier, not in p-blood: ", dad_carrier_num + dad_invisible_num )
        print("in p-tumor, in father, not in p-blood and mother:", dad_carrier_num)
        print("not in p-tumor, in father, not in p-blood and mother:", dad_invisible_num)
        print ("chimeral percentage: ",dad_carrier_num/dad_invisible_num*2/(1+dad_carrier_num/dad_invisible_num)) 
    
        print("mother is carrier, not in p-blood: ", mom_carrier_num + mom_invisible )
        print("in p-tumor, in mother, not in p-blood and father:", mom_carrier_num)
        print("not in p-tumor, in mother, not in p-blood and father:", mom_invisible)
        print ("chimeral percentage: ",mom_carrier_num/mom_invisible*2/(1+mom_carrier_num/mom_invisible)) 
    
        print("both are carriers, not in p-blood: ", both_num + both_invisible )
        print("in p-tumor, in mother and father, not in p-blood :", both_num)
        print(" not in p-tumor, in mother and father, not in p-blood :", both_invisible)
        print ("chimeral percentage: ",both_num/both_invisible*4/(3+3*both_num/both_invisible)) 
    
        return somatic, both_num, mom_carrier_num, dad_carrier_num, dad_invisible_num, both, mom_carrier, dad_carrier,dad_invisible
    
    










