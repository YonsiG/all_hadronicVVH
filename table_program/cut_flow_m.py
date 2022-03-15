#!/usr/bin/env python
"""
Cutflow table information 
"""

__author__ = "Yanxi Gu <yag002@ucsd.edu>"
__auther_m__ = "Luxin Zhang~ mua~ <zlx2016@mail.ustc.edu.cn>"
__copyright__ = "Copyright (c) Yanxi Gu"
__created__ = "[2022-01-26 Fri 01:39]" 

import sys,os,copy
import math
import ROOT 
#from tools import check_outfile_path

#ZZ_Selection = int(sys.argv[1])

def make_table(category, bkgname, samplename1, samplename2, samplename3, samplename4, outname, table_name):

    signal_samplename1 = '../../outfiles/' + bkgname + samplename1 + '_2_scaled.root'
    signal_samplename2 = '../../outfiles/' + bkgname + samplename2 + '_2_scaled.root'
    signal_samplename3 = '../../outfiles/' + bkgname + samplename3 + '_2_scaled.root'
    signal_samplename4 = '../../outfiles/' + bkgname + samplename4 + '_2_scaled.root'
    
    signal_sample1 =  ROOT.TFile(signal_samplename1)
    signal_sample2 =  ROOT.TFile(signal_samplename2)
    signal_sample3 =  ROOT.TFile(signal_samplename3)
    signal_sample4 =  ROOT.TFile(signal_samplename4)

    temp1 = signal_sample1.Get('cutflow%i'%(category))
    temp2 = signal_sample2.Get('cutflow%i'%(category))
    temp3 = signal_sample3.Get('cutflow%i'%(category))
    temp4 = signal_sample4.Get('cutflow%i'%(category))

    for i in range(20):
        exec ("signal1%s = temp1.GetBinContent(%s)"%(i+1,i+1))
        exec ("signal2%s = temp2.GetBinContent(%s)"%(i+1,i+1))
        exec ("signal3%s = temp3.GetBinContent(%s)"%(i+1,i+1))
        exec ("signal4%s = temp4.GetBinContent(%s)"%(i+1,i+1))

    fout_script = open(outname,'w')

    print('GYXILU\n')
    print('mua~\n')
    fout_script.write('GYXILU\n')
    fout_script.write('mua$\\tilde$\n')
    fout_script.write('\\begin{table}[]\n')
    fout_script.write('\\begin{tabular}{c|c|c|c|c}\n')

    print("%-30s%-25s%-25s%-25s%-25s"%('Cut',samplename1,samplename2,samplename3,samplename4))

    print("%-30s%-25s%-25s%-25s%-25s"%('Total (weighted to 137fb-1)',round(signal11,3),round(signal21,3),round(signal31,3),round(signal41,3)))
    print("%-30s%-25s%-25s%-25s%-25s"%('Selected Hbb at Gen level  ',round(signal12,3),round(signal22,3),round(signal32,3),round(signal42,3)))
    print("%-30s%-25s%-25s%-25s%-25s"%('Trigger                    ',round(signal13,3),round(signal23,3),round(signal33,3),round(signal43,3)))
    print("%-30s%-25s%-25s%-25s%-25s"%('0 lepton                   ',round(signal14,3),round(signal24,3),round(signal34,3),round(signal44,3)))
    if category==0:
        print("%-30s%-25s%-25s%-25s%-25s"%('3+ fatjets                 ',round(signal15,3),round(signal25,3),round(signal35,3),round(signal45,3)))
        print("%-30s%-25s%-25s%-25s%-25s"%('2+ jets                    ',round(signal16,3),round(signal26,3),round(signal36,3),round(signal46,3)))
    if category==3:
        print("%-30s%-25s%-25s%-25s%-25s"%('2+ fatjets                 ',round(signal15,3),round(signal25,3),round(signal35,3),round(signal45,3)))
        print("%-30s%-25s%-25s%-25s%-25s"%('4+ jets                    ',round(signal16,3),round(signal26,3),round(signal36,3),round(signal46,3)))
    if category==4:
        print("%-30s%-25s%-25s%-25s%-25s"%('2+ fatjets                 ',round(signal15,3),round(signal25,3),round(signal35,3),round(signal45,3)))
        print("%-30s%-25s%-25s%-25s%-25s"%('3+ jets                    ',round(signal16,3),round(signal26,3),round(signal36,3),round(signal46,3)))
    print("%-30s%-25s%-25s%-25s%-25s"%('VBF cut2                   ',round(signal17,3),round(signal27,3),round(signal37,3),round(signal47,3)))
#    print("%-30s%-25s%-25s%-25s%-25s"%('fatjet 0 Xbb_modified > 0.9',round(signal18,3),round(signal28,3),round(signal38,3),round(signal48,3)))
    print("%-30s%-25s%-25s%-25s%-25s"%('fatjet XbbMD modified > 0.9',round(signal18,3),round(signal28,3),round(signal38,3),round(signal48,3)))
    print("%-30s%-25s%-25s%-25s%-25s"%('fatjet 1 Xccqq modified > 0.9',round(signal19,3),round(signal29,3),round(signal39,3),round(signal49,3)))
    if category==0:
        print("%-30s%-25s%-25s%-25s%-25s"%('fatjet 2 Xccqq modified > 0.9      ',round(signal110,3),round(signal210,3),round(signal310,3),round(signal410,3)))
    print("%-30s%-25s%-25s%-25s%-25s"%('ST>1800                     ',round(signal111,3),round(signal211,3),round(signal311,3),round(signal411,3)))
    print("%-30s%-25s%-25s%-25s%-25s"%('VBF Etajj$>$4.5               ',round(signal112,3),round(signal212,3),round(signal312,3),round(signal412,3)))
    if category==0:
        print("%-30s%-25s%-25s%-25s%-25s"%('triboson mass $>$ 2200       ',round(signal113,3),round(signal213,3),round(signal313,3),round(signal413,3)))

    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('Cut',samplename1,samplename2,samplename3,samplename4,r'\\'))
    fout_script.write('\hline\n')
    fout_script.write('\hline\n')
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('Total (weighted to 137fb-1)',round(signal11,3),round(signal21,3),round(signal31,3),round(signal41,3),r'\\'))
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('Selected Hbb at Gen level  ',round(signal12,3),round(signal22,3),round(signal32,3),round(signal42,3),r'\\'))
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('Trigger                    ',round(signal13,3),round(signal23,3),round(signal33,3),round(signal43,3),r'\\'))
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('0 lepton                   ',round(signal14,3),round(signal24,3),round(signal34,3),round(signal44,3),r'\\'))
    if category==0:
        fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('3+ fatjets                 ',round(signal15,3),round(signal25,3),round(signal35,3),round(signal45,3),r'\\'))
        fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('2+ jets                    ',round(signal16,3),round(signal26,3),round(signal36,3),round(signal46,3),r'\\'))
    if category==3:
        fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('2+ fatjets                 ',round(signal15,3),round(signal25,3),round(signal35,3),round(signal45,3),r'\\'))
        fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('4+ jets                    ',round(signal16,3),round(signal26,3),round(signal36,3),round(signal46,3),r'\\'))
    if category==4:
        fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('2+ fatjets                 ',round(signal15,3),round(signal25,3),round(signal35,3),round(signal45,3),r'\\'))
        fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('3+ jets                    ',round(signal16,3),round(signal26,3),round(signal36,3),round(signal46,3),r'\\'))
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('VBF cut2                   ',round(signal17,3),round(signal27,3),round(signal37,3),round(signal47,3),r'\\'))
#    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('fatjet 0 Xbb$\_$modified $>$ 0.9',round(signal18,3),round(signal28,3),round(signal38,3),round(signal48,3),r'\\'))
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('fatjet 0 XbbMD modified $>$ 0.9',round(signal18,3),round(signal28,3),round(signal38,3),round(signal48,3),r'\\'))
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('fatjet 1 Xccqq modified $>$ 0.9',round(signal19,3),round(signal29,3),round(signal39,3),round(signal49,3),r'\\'))
    if category==0:
        fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('fatjet 2 Xccqq modified $>$ 0.9  ',round(signal110,3),round(signal210,3),round(signal310,3),round(signal410,3),r'\\'))
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('ST$>$1800                  ',round(signal111,3),round(signal211,3),round(signal311,3),round(signal411,3),r'\\'))
    fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('VBF Etajj$>$4.5              ',round(signal112,3),round(signal212,3),round(signal312,3),round(signal412,3),r'\\'))
    if category==0:
        fout_script.write('%-30s&%-25s&%-25s&%-25s&%-25s%-5s\n'%('triboson mass $>$ 2200  ',round(signal113,3),round(signal213,3),round(signal313,3),round(signal413,3),r'\\'))

    fout_script.write('\\end{tabular}\n')
    fout_script.write('\\caption{'+table_name+'}\n')
    fout_script.write('\\label{'+table_name+'}\n')
    fout_script.write('~\\\\\n')
    fout_script.write('\\end{table}\n')

    fout_script.close()

for icategory in [0,3,4]:
    if icategory == 0:
        Soutname_in = '../../tables/cut_info_signal_3+2.txt'
        Woutname_in = '../../tables/cut_info_WJet_3+2.txt'
        Zoutname_in = '../../tables/cut_info_ZJet_3+2.txt'
        QCDoutname_in = '../../tables/cut_info_QCDJet_3+2.txt'
        ttbaroutname_in = '../../tables/cut_info_ttbarJet_3+2.txt'
    if icategory == 3:
        Soutname_in = '../../tables/cut_info_signal_2+4.txt'
        Woutname_in = '../../tables/cut_info_WJet_2+4.txt'
        Zoutname_in = '../../tables/cut_info_ZJet_2+4.txt'
        QCDoutname_in = '../../tables/cut_info_QCDJet_2+4.txt'
        ttbaroutname_in = '../../tables/cut_info_ttbarJet_2+4.txt'
    if icategory == 4:
        Soutname_in = '../../tables/cut_info_signal_2+3.txt'
        Woutname_in = '../../tables/cut_info_WJet_2+3.txt'
        Zoutname_in = '../../tables/cut_info_ZJet_2+3.txt'
        QCDoutname_in = '../../tables/cut_info_QCDJet_2+3.txt'
        ttbaroutname_in = '../../tables/cut_info_ttbarJet_2+3.txt'
    Stable_name_in = "signal"
    Wtable_name_in = "WJetsToQQ$\_$HT"
    Ztable_name_in = "ZJetsToQQ$\_$HT"
    QCDtable_name_in = "QCD$\_$HT"
    ttbartable_name_in = "ttbar"
    make_table(icategory,"", "WZH", "ZZH", "SSWWH", "OSWWH", Soutname_in, Stable_name_in)
    make_table(icategory,"WJetsToQQ_HT-", "200to400", "400to600", "600to800", "800toInf", Woutname_in, Wtable_name_in)
    make_table(icategory,"ZJetsToQQ_HT-", "200to400", "400to600", "600to800", "800toInf", Zoutname_in, Ztable_name_in)
    make_table(icategory,"QCD_HT", "300to500", "500to700", "700to1000", "1000to1500", QCDoutname_in, QCDtable_name_in)
    make_table(icategory,"", "QCD_HT1500to2000", "QCD_HT2000toInf", "TTToHadronic", "TTToSemiLeptonic", ttbaroutname_in, ttbartable_name_in)



