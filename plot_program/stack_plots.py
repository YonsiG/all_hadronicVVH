import ROOT

def draw_plot(plotname="fatjet_msoftdrop", title="myTitle", log=True):
    #open file
    signalfile = ROOT.TFile("../../outfiles/signal_merged_2_scaled.root")
    QCDfile = ROOT.TFile("../../outfiles/QCD_merged_2_scaled.root")
    TTToHadronicfile = ROOT.TFile("../../outfiles/TTToHadronic_2_scaled.root")
    TTToSemiLeptonicfile = ROOT.TFile("../../outfiles/TTToSemiLeptonic_2_scaled.root")
    WJetfile = ROOT.TFile("../../outfiles/WJet_merged_2_scaled.root")
    ZJetfile = ROOT.TFile("../../outfiles/ZJet_merged_2_scaled.root")

    #get historam
    signalplot = signalfile.Get(plotname)
    QCDplot = QCDfile.Get(plotname)
    TTToHadronicplot = TTToHadronicfile.Get(plotname)
    TTToSemiLeptonicplot = TTToSemiLeptonicfile.Get(plotname)
    WJetplot = WJetfile.Get(plotname)
    ZJetplot = ZJetfile.Get(plotname)

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
    legend.AddEntry(signalplot,"signal")
    legend.AddEntry(QCDplot, "QCD")
    legend.AddEntry(TTToHadronicplot, "TTToHadronic")
    legend.AddEntry(TTToSemiLeptonicplot, "TTToSemiLeptonic")
    legend.AddEntry(WJetplot,"WJet")
    legend.AddEntry(ZJetplot,"ZJet")

    if log==True:
        QCDplot.GetYaxis().SetRangeUser(10e-2,10e8)
        TTToHadronicplot.GetYaxis().SetRangeUser(10e-2,10e8)
        TTToSemiLeptonicplot.GetYaxis().SetRangeUser(10e-2,10e8)
        WJetplot.GetYaxis().SetRangeUser(10e-2,10e8)
        ZJetplot.GetYaxis().SetRangeUser(10e-2,10e8)
        signalplot.GetYaxis().SetRangeUser(10e-2,10e8)

    #define canvas
    canvas = ROOT.TCanvas("canvas","canvas",800,800)
    if log==True:
        ROOT.gPad.SetLogy(1)

    #plot data,stack, signal  
    if log==True:
        signalplot.Draw("Hist")
        stack.Draw("Hist same")
        signalplot.Draw("H same")

    if log==False:
        stack.Draw("Hist")
        signalplot.Draw("H same")
    legend.Draw()

    #print and save
    if log==True:
        canvas.SaveAs("../../plots/" + plotname + "_s+b_log.pdf")
    if log==False:
        canvas.SaveAs("../../plots/" + plotname + "_s+b_linear.pdf")


ROOT.gStyle.SetOptStat(0000)

listofplots=["fatjet_msoftdrop", "fatjet_pt", "fatjet_eta", "VBF_max_mass"]

for plot in listofplots:
    title=plot
    draw_plot(plot, title, True)
    draw_plot(plot, title, False)