#!/usr/bin/env python
"""
Cutflow table information 
"""

__author__ = "Yanxi Gu <yag002@ucsd.edu>"
#"and Luxin Zhang~ mua~"
__copyright__ = "Copyright (c) Yanxi Gu"
__created__ = "[2022-01-26 Fri 01:39]" 

import sys,os,copy
import math
import ROOT 
#from tools import check_outfile_path

#ZZ_Selection = int(sys.argv[1])

def main():

    #tabs = sys.argv[3:]

    signal_sample1 =  ROOT.TFile('../../outfiles/WJetsToQQ_HT-200to400_2_scaled.root')
    signal_sample2 =  ROOT.TFile('../../outfiles/WJetsToQQ_HT-400to600_2_scaled.root')
    signal_sample3 =  ROOT.TFile('../../outfiles/WJetsToQQ_HT-600to800_2_scaled.root')
    signal_sample4 =  ROOT.TFile('../../outfiles/WJetsToQQ_HT-800toInf_2_scaled.root')

    temp1 = signal_sample1.Get('cutflow0')
    temp2 = signal_sample2.Get('cutflow0')
    temp3 = signal_sample3.Get('cutflow0')
    temp4 = signal_sample4.Get('cutflow0')

    for i in range(20):
        exec ("signal1%s = temp1.GetBinContent(%s)"%(i+1,i+1))
        exec ("signal2%s = temp2.GetBinContent(%s)"%(i+1,i+1))
        exec ("signal3%s = temp3.GetBinContent(%s)"%(i+1,i+1))
        exec ("signal4%s = temp4.GetBinContent(%s)"%(i+1,i+1))

    outname = '../../tables/cut_info_WJet.txt'
    table_name = "WJetsToQQ$\_$HT"

    fout_script = open(outname,'w')

    print('GYXILU\n')
    print('mua~\n')
    fout_script.write('GYXILU\n')
    fout_script.write('mua$\\tilde$\n')
    fout_script.write('\\begin{table}[]\n')
    fout_script.write('\\begin{tabular}{c|c|c|c|c}\n')

    print("%-30s%-25s%-25s%-25s%-25s"%('Cut','WJetsToQQ_HT-200to400','WJetsToQQ_HT-400to600','WJetsToQQ_HT-600to800','WJetsToQQ_HT-800toInf'))

    print("%-30s%-25s%-25s%-25s%-25s"%('Total (weighted to 137fb-1)',round(signal11,3),round(signal21,3),round(signal31,3),round(signal41,3)))
    print("%-30s%-25s%-25s%-25s%-25s"%('Selected Hbb at Gen level  ',round(signal12,3),round(signal22,3),round(signal32,3),round(signal42,3)))
    print("%-30s%-25s%-25s%-25s%-25s"%('Trigger                    ',round(signal13,3),round(signal23,3),round(signal33,3),round(signal43,3)))
    print("%-30s%-25s%-25s%-25s%-25s"%('0 lepton                   ',round(signal14,3),round(signal24,3),round(signal34,3),round(signal44,3)))
    print("%-30s%-25s%-25s%-25s%-25s"%('3+ fatjets                 ',round(signal15,3),round(signal25,3),round(signal35,3),round(signal45,3)))
    print("%-30s%-25s%-25s%-25s%-25s"%('2+ jets                    ',round(signal16,3),round(signal26,3),round(signal36,3),round(signal46,3)))
    print("%-30s%-25s%-25s%-25s%-25s"%('VBF cut2                   ',round(signal17,3),round(signal27,3),round(signal37,3),round(signal47,3)))

    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('Cut','200to400','400to600','600to800','800toInf',r'\\'))
    fout_script.write('\hline\n')
    fout_script.write('\hline\n')
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('Total (weighted to 137fb-1)',round(signal11,3),round(signal21,3),round(signal31,3),round(signal41,3),r'\\'))
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('Selected Hbb at Gen level  ',round(signal12,3),round(signal22,3),round(signal32,3),round(signal42,3),r'\\'))
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('Trigger                    ',round(signal13,3),round(signal23,3),round(signal33,3),round(signal43,3),r'\\'))
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('0 lepton                   ',round(signal14,3),round(signal24,3),round(signal34,3),round(signal44,3),r'\\'))
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('3+ fatjets                 ',round(signal15,3),round(signal25,3),round(signal35,3),round(signal45,3),r'\\'))
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('2+ jets                    ',round(signal16,3),round(signal26,3),round(signal36,3),round(signal46,3),r'\\'))
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('VBF cut2                   ',round(signal17,3),round(signal27,3),round(signal37,3),round(signal47,3),r'\\'))

    fout_script.write('\\end{tabular}\n')
    fout_script.write('\\caption{'+table_name+'}\n')
    fout_script.write('\\label{'+table_name+'}\n')
    fout_script.write('\\end{table}\n')

    fout_script.close()
    
if __name__ == '__main__':
    main()
