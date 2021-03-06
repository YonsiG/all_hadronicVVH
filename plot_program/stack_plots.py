from curses import can_change_color
from hashlib import new
from numbers import Integral
from tkinter.tix import Tree
import ROOT
import copy

def draw_plot(plotname="fatjet_msoftdrop", title="myTitle", log=True, compare_data=True, DoRatio=True):
    #open file
    signalfile = ROOT.TFile("../../outfiles/signal_merged_2_scaled.root")
    QCDfile = ROOT.TFile("../../outfiles/QCD_merged_2_scaled.root")
    TTToHadronicfile = ROOT.TFile("../../outfiles/TTToHadronic_2_scaled.root")
    TTToSemiLeptonicfile = ROOT.TFile("../../outfiles/TTToSemiLeptonic_2_scaled.root")
    WJetfile = ROOT.TFile("../../outfiles/WJet_merged_2_scaled.root")
    ZJetfile = ROOT.TFile("../../outfiles/ZJet_merged_2_scaled.root")
    if compare_data:
        datafile = ROOT.TFile("../../outfiles/data_2018_2_selected.root")

    #get historam
    signalplot = signalfile.Get(plotname)
    if compare_data:
        dataplot = datafile.Get(plotname)
        dataplot.Rebin(4)
    QCDplot = QCDfile.Get(plotname)
    TTToHadronicplot = TTToHadronicfile.Get(plotname)
    TTToSemiLeptonicplot = TTToSemiLeptonicfile.Get(plotname)
    WJetplot = WJetfile.Get(plotname)
    ZJetplot = ZJetfile.Get(plotname)

    signalplot.Rebin(4)
    QCDplot.Rebin(4)
    TTToHadronicplot.Rebin(4)
    TTToSemiLeptonicplot.Rebin(4)
    WJetplot.Rebin(4)
    ZJetplot.Rebin(4)

    #draw overflow
    signalplot.SetBinContent(1, signalplot.GetBinContent(1) + signalplot.GetBinContent(0))
    signalplot.SetBinContent(signalplot.GetNbinsX(), signalplot.GetBinContent(signalplot.GetNbinsX() + 1) + signalplot.GetBinContent(signalplot.GetNbinsX()))
    if compare_data:
        dataplot.SetBinContent(1, dataplot.GetBinContent(1) + dataplot.GetBinContent(0))
        dataplot.SetBinContent(dataplot.GetNbinsX(), dataplot.GetBinContent(dataplot.GetNbinsX() + 1) + dataplot.GetBinContent(dataplot.GetNbinsX()))
    QCDplot.SetBinContent(1, QCDplot.GetBinContent(1) + QCDplot.GetBinContent(0))
    QCDplot.SetBinContent(QCDplot.GetNbinsX(), QCDplot.GetBinContent(QCDplot.GetNbinsX() + 1) + QCDplot.GetBinContent(QCDplot.GetNbinsX()))
    TTToHadronicplot.SetBinContent(1, TTToHadronicplot.GetBinContent(1) + TTToHadronicplot.GetBinContent(0))
    TTToHadronicplot.SetBinContent(TTToHadronicplot.GetNbinsX(), TTToHadronicplot.GetBinContent(TTToHadronicplot.GetNbinsX() + 1) + TTToHadronicplot.GetBinContent(TTToHadronicplot.GetNbinsX()))
    TTToSemiLeptonicplot.SetBinContent(1, TTToSemiLeptonicplot.GetBinContent(1) + TTToSemiLeptonicplot.GetBinContent(0))
    TTToSemiLeptonicplot.SetBinContent(TTToSemiLeptonicplot.GetNbinsX(), TTToSemiLeptonicplot.GetBinContent(TTToSemiLeptonicplot.GetNbinsX() + 1) + TTToSemiLeptonicplot.GetBinContent(TTToSemiLeptonicplot.GetNbinsX()))
    WJetplot.SetBinContent(1, WJetplot.GetBinContent(1) + WJetplot.GetBinContent(0))
    WJetplot.SetBinContent(WJetplot.GetNbinsX(), WJetplot.GetBinContent(WJetplot.GetNbinsX() + 1) + WJetplot.GetBinContent(WJetplot.GetNbinsX()))
    ZJetplot.SetBinContent(1, ZJetplot.GetBinContent(1) + ZJetplot.GetBinContent(0))
    ZJetplot.SetBinContent(ZJetplot.GetNbinsX(), ZJetplot.GetBinContent(ZJetplot.GetNbinsX() + 1) + ZJetplot.GetBinContent(ZJetplot.GetNbinsX()))

    #buid stack
    stack = ROOT.THStack("stack","")
    WJetplot.SetFillColor(ROOT.kOrange+3)
    stack.Add(WJetplot)
    ZJetplot.SetFillColor(ROOT.kOrange+4)
    stack.Add(ZJetplot)
    TTToHadronicplot.SetFillColor(ROOT.kOrange+1)
    stack.Add(TTToHadronicplot)
    TTToSemiLeptonicplot.SetFillColor(ROOT.kOrange+2)
    stack.Add(TTToSemiLeptonicplot)
    QCDplot.SetFillColor(ROOT.kOrange)
    stack.Add(QCDplot)
    stack.SetTitle(title)

    #plot legends, ranges
    legend = ROOT.TLegend(0.7,0.83,0.95,0.97)
    legend.SetTextFont(60)
    legend.SetTextSize(0.02)
    legend.AddEntry(signalplot,"signal %.2f"%(signalplot.Integral()))
    if compare_data:
        legend.AddEntry(dataplot,"data %.2f"%(dataplot.Integral()))
    legend.AddEntry(QCDplot, "QCD %.2f"%(QCDplot.Integral()))
    legend.AddEntry(TTToHadronicplot, "TTToHadronic %.2f"%(TTToHadronicplot.Integral()))
    legend.AddEntry(TTToSemiLeptonicplot, "TTToSemiLeptonic %.2f"%(TTToSemiLeptonicplot.Integral()))
    legend.AddEntry(WJetplot,"WJet %.2f"%(WJetplot.Integral()))
    legend.AddEntry(ZJetplot,"ZJet %.2f"%(ZJetplot.Integral()))

    if log==True:
        QCDplot.GetYaxis().SetRangeUser(10e-2,10e8)
        TTToHadronicplot.GetYaxis().SetRangeUser(10e-2,10e8)
        TTToSemiLeptonicplot.GetYaxis().SetRangeUser(10e-2,10e8)
        WJetplot.GetYaxis().SetRangeUser(10e-2,10e8)
        ZJetplot.GetYaxis().SetRangeUser(10e-2,10e8)
        signalplot.GetYaxis().SetRangeUser(10e-2,10e8)
        if compare_data:
            dataplot.GetYaxis().SetRangeUser(10e-2,10e8)

    #define canvas
    canvas = ROOT.TCanvas("canvas","canvas",800,800)

    if DoRatio==True:
        MCplot = copy.deepcopy(QCDplot)
        MCplot.Add(TTToHadronicplot)
        MCplot.Add(TTToSemiLeptonicplot)
        MCplot.Add(WJetplot)
        MCplot.Add(ZJetplot)
        ratioplot=copy.deepcopy(dataplot)
        ratioplot.Divide(MCplot)
        ratioplot.SetTitle(";"+title+";data / MC")
        pad1 = ROOT.TPad("pad1","pad1",0,0.3,1,1)
        pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.3)
        pad1.Draw()
        pad2.Draw()
        pad2.cd()
        ratioplot.Draw("E0")

    if DoRatio==False:
        pad1 = ROOT.TPad("pad1","pad1",0,0,1,1)
        pad1.Draw()

    pad1.cd()
    if log==True:
        pad1.SetLogy()
        #ROOT.gPad.SetLogy(1)

    #plot data,stack, signal, data  
    if log==True:
        signalplot.Draw("Hist")
        stack.Draw("Hist same")
        signalplot.Draw("H same")
        if compare_data:
            dataplot.Draw("E0 same")

    if log==False:
        stack.Draw("Hist")
        signalplot.Draw("H same")
        if compare_data:
            dataplot.Draw("E0 same")
    legend.Draw()

    #print and save
    if log==True:
        if compare_data==False:
            canvas.SaveAs("../../plots/" + plotname + "_s+b_log.png")
        if compare_data==True:
            canvas.SaveAs("../../plots/" + plotname + "_mc+data_log.png")
    if log==False:
        if compare_data==False:
            canvas.SaveAs("../../plots/" + plotname + "_s+b_linear.png")
        if compare_data==True:
            canvas.SaveAs("../../plots/" + plotname + "_mc+data_linear.png")

ROOT.gStyle.SetOptStat(0000)

listofplots1=["fatjet_msoftdrop0_0", "fatjet_msoftdrop0_1", "fatjet_msoftdrop0_2",
                "fatjet_pt0_0", "fatjet_pt0_1", "fatjet_pt0_2",
                "fatjet_eta0_0", "fatjet_eta0_1", "fatjet_eta0_2",
                "fatjet_WvsQCD0_0", "fatjet_WvsQCD0_1", "fatjet_WvsQCD0_2",
                "fatjet_ZvsQCD0_0", "fatjet_ZvsQCD0_1", "fatjet_ZvsQCD0_2",
                #"fatjet_mass0_0", "fatjet_mass0_1", "fatjet_mass0_2",
                "fatjet_Xbb_modified0_0", "fatjet_Xbb_modified0_1", "fatjet_Xbb_modified0_2",
                "fatjet_Xcc0_0", "fatjet_Xcc0_1", "fatjet_Xcc0_2",
                "fatjet_Xqq0_0", "fatjet_Xqq0_1", "fatjet_Xqq0_2",
                "fatjet_QCD0_0", "fatjet_QCD0_1", "fatjet_QCD0_2",
                "fatjet_Xccqq_modified0_0", "fatjet_Xccqq_modified0_1", "fatjet_Xccqq_modified0_2",
                "VBF_max_mass0", "VBF_deltaEta0",
                "ST0", "mVVV0"]

for plot in listofplots1:
    title=plot
    draw_plot(plot, title, False, True, True)
    draw_plot(plot, title, True, True, True)

listofplots2=["fatjet_msoftdrop3_0", "fatjet_msoftdrop3_1", "fatjet_msoftdrop3_2",
                "fatjet_pt3_0", "fatjet_pt3_1", "fatjet_pt3_2",
                "fatjet_eta3_0", "fatjet_eta3_1", "fatjet_eta3_2",
                "fatjet_WvsQCD3_0", "fatjet_WvsQCD3_1", "fatjet_WvsQCD3_2",
                "fatjet_ZvsQCD0_0", "fatjet_ZvsQCD0_1", "fatjet_ZvsQCD0_2",
                #"fatjet_mass3_0", "fatjet_mass3_1", "fatjet_mass3_2",
                "fatjet_Xbb_modified3_0", "fatjet_Xbb_modified3_1", "fatjet_Xbb_modified3_2",
                "fatjet_Xcc3_0", "fatjet_Xcc3_1", "fatjet_Xcc3_2",
                "fatjet_Xqq3_0", "fatjet_Xqq3_1", "fatjet_Xqq3_2",
                "fatjet_QCD3_0", "fatjet_QCD3_1", "fatjet_QCD3_2",
                "fatjet_Xccqq_modified3_0", "fatjet_Xccqq_modified3_1", "fatjet_Xccqq_modified3_2",
                "VBF_max_mass3", "VBF_deltaEta3",
                "ST3", "mVVV3"]

for plot in listofplots2:
    title=plot
    draw_plot(plot, title, True, True, True)
    draw_plot(plot, title, False, True, True)

listofplots3=["fatjet_msoftdrop4_0", "fatjet_msoftdrop4_1", "fatjet_msoftdrop4_2",
                "fatjet_pt4_0", "fatjet_pt4_1", "fatjet_pt4_2",
                "fatjet_eta4_0", "fatjet_eta4_1", "fatjet_eta4_2",
                "fatjet_WvsQCD4_0", "fatjet_WvsQCD4_1", "fatjet_WvsQCD4_2",
                "fatjet_ZvsQCD0_0", "fatjet_ZvsQCD0_1", "fatjet_ZvsQCD0_2",
                #"fatjet_mass4_0", "fatjet_mass4_1", "fatjet_mass4_2",
                "fatjet_Xbb_modified4_0", "fatjet_Xbb_modified4_1", "fatjet_Xbb_modified4_2",
                "fatjet_Xcc4_0", "fatjet_Xcc4_1", "fatjet_Xcc4_2",
                "fatjet_Xqq4_0", "fatjet_Xqq4_1", "fatjet_Xqq4_2",
                "fatjet_QCD4_0", "fatjet_QCD4_1", "fatjet_QCD4_2",
                "fatjet_Xccqq_modified4_0", "fatjet_Xccqq_modified4_1", "fatjet_Xccqq_modified4_2",
                "VBF_max_mass4", "VBF_deltaEta4",
                "ST4", "mVVV4"]

for plot in listofplots3:
    title=plot
    draw_plot(plot, title, True, True, True)
    draw_plot(plot, title, False, True, True)