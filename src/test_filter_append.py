#!/usr/local/bin/python3.8

import filter_append as fa
import unittest


class TestFilter_append(unittest.TestCase):

        
    with open("../data/VQSR_filtered/snp.recalibrated.vcf",'r') as f:
        data_in = f.read().rstrip().split('\n') # Read file and create list split by newline 

    filtered=fa.FilterAppend(data_in) # build the FilterAppend instance
        
    header_line_number, header = filtered.find_header()

    passed=filtered.filter_PASS(header_line_number=header_line_number)
    
    filtered1=fa.FilterAppend(passed) # build a new instance
    AF_filtered=filtered1.filter_AF()
    
    filtered2=fa.FilterAppend(AF_filtered) # build a new instance
    high_quality_snp=filtered2.append_AF()
    print("High quality SNPs:",len(high_quality_snp)) 
    
    filtered3=fa.FilterAppend(high_quality_snp) # build a new instance
    somatic, both_num, mom_carrier_num, dad_carrier_num, both, mom_carrier, dad_carrier=filtered3.find_chimera()

        

    def test_find_header(self):

        header_line_number, header = self.filtered.find_header()  
        self.assertEqual(len(header), header_line_number+1)
        
    def test_filter_PASS(self,header_line_number=header_line_number):

        passed = self.filtered.filter_PASS(header_line_number=header_line_number)  
        self.assertIsNotNone(passed)
        self.assertIn("PASS", passed[1])
        
    def test_filter_AF(self):

        AF_filtered = self.filtered1.filter_AF()  
        self.assertIsNotNone(AF_filtered)
        
    def test_append_AF(self):

        high_quality_snp = self.filtered2.append_AF()  
        self.assertIsNotNone(high_quality_snp)
        self.assertEqual(len(high_quality_snp[1].split("\t")),24)
        
    def test_find_chimera(self):

        somatic, both_num, mom_carrier_num, dad_carrier_num, both, mom_carrier, dad_carrier = self.filtered3.find_chimera()  
        self.assertIsNotNone(somatic)
        self.assertEqual(len(mom_carrier[1].split("\t")),25)

if __name__ == '__main__':
    # begin the unittest.main()
    unittest.main()
    
    
## To run unittest
## python -m unittest -v test_file.py    
## python -m unittest test_file.py    



