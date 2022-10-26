#!/usr/local/bin/python3.8

import unittest

import 
for name in glob.glob("./snp.recalibrated.vcf"): #input is the vcf file 
        inputfile = str(name)
        outputfile = str(name)[:-17]+"_VQSR_AFDP_filtered.txt"
        print(outputfile )
        
        with open(inputfile,'r') as f:
            data_in = f.read().rstrip().split('\n') # Read file and create list split by newline 
        print(len(data_in))
       
class TestCuboid(unittest.TestCase):
    def test_volume(self):
        self.assertAlmostEqual(cuboid_volume(2),8)
        self.assertAlmostEqual(cuboid_volume(1),1)
        self.assertAlmostEqual(cuboid_volume(0),0)
        self.assertAlmostEqual(cuboid_volume(5.5),166.375)