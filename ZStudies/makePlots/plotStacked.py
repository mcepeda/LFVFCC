#!/usr/bin/env python3

import argparse
import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# =============================================================================
# USER CONFIGURATION
# =============================================================================

LUMI_AB = 17.0
ECM_GEV = 91.2
BR_LFV = 5.0e-6

XSEC_ZTAUTAU_PB = 1476.58
BR_ZTAUTAU = 0.036
XSEC_Z_TOTAL_PB = XSEC_ZTAUTAU_PB / BR_ZTAUTAU
XSEC_SIGNAL_PB = XSEC_Z_TOTAL_PB * BR_LFV

CHANNELS = {
    "mutauh": {
        "channel_label": "#mu-#tau_{h}",
        "signal_label": "Z #rightarrow #mu#tau",
#        "background_file": "histos_ZTauTau_mutauh_long_sel.root",
#        "background_ngen": 1430000,
        "background_file": "histos_ZTauTau_PY8_mutauh_long_sel.root",
        "background_ngen": 100000000,
        "signal_file": "histos_ZMuTau_mutauh_long_sel.root",
        "signal_ngen": 2e6,
    },

    "emu": {
        "channel_label": "e-#mu",
        "signal_label": "Z #rightarrow #mu#tau",
        "background_file": "histos_ZTauTau_emu_long_sel.root",
        "background_ngen": 1430000,
#        "background_file": "histos_ZTauTau_Py8_emu_long_sel.root",
#        "background_ngen": 20600000,
        "signal_file": "histos_ZMuTau_emu_long_sel.root",
        "signal_ngen": 2e6,
    },
}

HISTOGRAMS = {
    "before": "selection_constrainedMass_before",
    "after": "selection_constrainedMass_after",
    "poverbeam": "obj1_poverebeam",
    "pbefore": "selection_obj1_p_before",
    "pafter": "selection_obj1_p_after",
}   

BACKGROUND_COLOR = ROOT.kAzure - 9
SIGNAL_COLOR = ROOT.kRed + 1

REBIN_FACTOR = 2 

MINX=0.90 # 80 #0.90
MAXX=1.05 # 95 #1.1

SAVE_PNG = True
SAVE_PDF = True


def setup_style():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetPadLeftMargin(0.14)
    ROOT.gStyle.SetPadRightMargin(0.04)
    ROOT.gStyle.SetPadTopMargin(0.12)
    ROOT.gStyle.SetPadBottomMargin(0.14)
    ROOT.gStyle.SetLabelSize(0.045, "XY")
    ROOT.gStyle.SetTitleSize(0.052, "XY")
    ROOT.gStyle.SetTitleOffset(1.18, "Y")
    ROOT.gStyle.SetTitleOffset(1.05, "X")
    ROOT.gStyle.SetTickLength(0.025, "XY")
    ROOT.gStyle.SetNdivisions(510, "XY")


def load_histogram(filename, histname, clone_name):
    f = ROOT.TFile.Open(filename, "READ")

    if not f or f.IsZombie():
        raise RuntimeError("Could not open %s" % filename)

    h = f.Get(histname)

    if not h:
        f.Close()
        raise RuntimeError("Histogram '%s' not found in %s" % (histname, filename))

    h = h.Clone(clone_name)
    h.SetDirectory(0)
    f.Close()

    h.Rebin(REBIN_FACTOR)
   
    return h


def mc_scale(xsec_pb, lumi_ab, ngen):
    if ngen <= 0:
        raise RuntimeError("Ngen must be positive")
    return xsec_pb * lumi_ab * 1000* 1000.0 / float(ngen)


def y_axis_title(h):
    width = h.GetBinWidth(1)

    uniform = True
    for ibin in range(2, h.GetNbinsX() + 1):
        if abs(h.GetBinWidth(ibin) - width) > 1e-9:
            uniform = False
            break

    if not uniform:
        return "Events / bin"

    if width >= 1.0:
        return "Events / %.1f GeV" % width

    return "Events / %.2f GeV" % width


def draw_fcc_labels(channel_label, stage):
    keep = []

    top = ROOT.TLatex()
    top.SetNDC(True)
    top.SetTextAlign(31)
    top.SetTextSize(0.045)

    top.SetTextFont(62)
    top.DrawLatex(0.95, 0.955, "FCC-ee")

    top.SetTextFont(52)
    top.DrawLatex(0.95, 0.910, "(GEN Information)")
    keep.append(top)

    info = ROOT.TLatex()
    info.SetNDC(True)
    info.SetTextAlign(13)
    info.SetTextFont(42)
    info.SetTextSize(0.040)

    info.DrawLatex(0.17, 0.84, "#sqrt{s} = %.1f GeV" % ECM_GEV)
    info.DrawLatex(0.17, 0.79, "L = %.1f ab^{-1}" % LUMI_AB)
    info.DrawLatex(0.17, 0.74, "Channel: %s" % channel_label)
    info.DrawLatex(0.17, 0.69, "Before selection" if stage == "before" else "After selection")
    keep.append(info)

    return keep


def draw_one(channel, stage, config, outdir):
    histname = HISTOGRAMS[stage]

    hbkg = load_histogram(
        config["background_file"],
        histname,
        "hbkg_%s_%s" % (channel, stage)
    )

    hsig = load_histogram(
        config["signal_file"],
        histname,
        "hsig_%s_%s" % (channel, stage)
    )

    bkg_scale = mc_scale(
        XSEC_ZTAUTAU_PB,
        LUMI_AB,
        config["background_ngen"]
    )

    sig_scale = mc_scale(
        XSEC_SIGNAL_PB,
        LUMI_AB,
        config["signal_ngen"]
    )

    hbkg.Scale(bkg_scale)
    hsig.Scale(sig_scale)

    print("")
    print("%s  %s" % (channel, stage))
    print("  background scale/event = %.8g" % bkg_scale)
    print("  signal scale/event     = %.8g" % sig_scale)
    print("  background yield       = %.8g" % hbkg.Integral(0, hbkg.GetNbinsX()+1))
    print("  signal yield           = %.8g" % hsig.Integral(0, hsig.GetNbinsX()+1))

    # Red line is background + LFV signal, matching the visual logic
    # of the FCC-ee reference plot.
    htotal = hbkg.Clone("htotal_%s_%s" % (channel, stage))
    htotal.SetDirectory(0)
    htotal.Add(hsig)

    hbkg.SetFillColor(BACKGROUND_COLOR)
    hbkg.SetLineColor(ROOT.kBlue + 2)
    hbkg.SetLineWidth(1)

    hsig.SetFillStyle(0)
    hsig.SetLineColor(SIGNAL_COLOR)
    hsig.SetLineWidth(3)

    htotal.SetFillStyle(0)
    htotal.SetLineColor(ROOT.kBlack)
    htotal.SetLineWidth(1)

    hbkg.SetTitle("")
    hbkg.GetXaxis().SetTitle("Muon p/E") #"m_{Z} [GeV]")
    hbkg.GetYaxis().SetTitle(y_axis_title(hbkg))

    ymax = max(hbkg.GetMaximum(), htotal.GetMaximum())
    hbkg.SetMinimum(0.0)
    hbkg.SetMaximum(1.35 * ymax if ymax > 0 else 1.0)

    hbkg.GetXaxis().SetRangeUser(MINX,MAXX)

    canvas = ROOT.TCanvas("c_%s_%s" % (channel, stage), "", 900, 700)
    canvas.SetTicks(1, 1)

    hbkg.Draw("HIST")
    hsig.Draw("HIST SAME")
    htotal.Draw("HIST SAME")

    legend = ROOT.TLegend(0.58, 0.69, 0.90, 0.84)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.037)

    legend.AddEntry(hbkg, "Z #rightarrow #tau#tau", "f")
    legend.AddEntry(
        hsig,
        "%s  (BR = %.1e)" %
        (config["signal_label"], BR_LFV),
        "l"
    )
    legend.AddEntry(htotal, "Z #rightarrow #tau#tau + LFV signal", "l")
    legend.Draw()

    labels = draw_fcc_labels(config["channel_label"], stage)

    ROOT.gPad.RedrawAxis()

    output_base = os.path.join(outdir, "ZMass_%s_%s" % (channel, stage))

    if SAVE_PNG:
        canvas.SaveAs(output_base + ".png")

    if SAVE_PDF:
        canvas.SaveAs(output_base + ".pdf")

    canvas.Close()

    _keepalive = [hbkg, hsig, htotal, legend] + labels


def main():
    parser = argparse.ArgumentParser(
        description="FCC-ee style mass plots for LFV studies"
    )

    parser.add_argument(
        "-c", "--channel",
        required=True,
        choices=sorted(CHANNELS.keys())
    )

    parser.add_argument(
        "-v", "--version",
        default="V1"
    )

    args = parser.parse_args()

    setup_style()

    config = CHANNELS[args.channel]

    outdir = "plots_Mass_%s_%s" % (args.channel, args.version)
    os.makedirs(outdir, exist_ok=True)

    print("==================================================")
    print("Visible-mass FCC-ee GEN plot")
    print("Channel            :", args.channel)
    print("sqrt(s)            : %.1f GeV" % ECM_GEV)
    print("Luminosity         : %.1f ab^-1" % LUMI_AB)
    print("BR(LFV)            : %.3g" % BR_LFV)
    print("sigma(Z->tautau)   : %.6g pb" % XSEC_ZTAUTAU_PB)
    print("sigma(Z total)     : %.6g pb" % XSEC_Z_TOTAL_PB)
    print("sigma(signal)      : %.6g pb" % XSEC_SIGNAL_PB)
    print("Output             :", outdir)
    print("==================================================")

    draw_one(args.channel, "before", config, outdir)
    draw_one(args.channel, "after", config, outdir)
    draw_one(args.channel, "poverbeam", config, outdir)
    draw_one(args.channel, "pbefore", config, outdir)
    draw_one(args.channel, "pafter", config, outdir)


    print("")
    print("Done:", outdir)


if __name__ == "__main__":
    main()

