#!/usr/bin/env python3

import argparse
import math
import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

##### python plotLFVHistos.py -c 'mutauh' -v V5

# Shape comparison is usually the useful default for signal/background studies.
# Set to False if you want raw event yields instead.
NORMALIZE = False



# =============================================================================
# USER CONFIGURATION
# =============================================================================
#
# Edit this block to define the samples you want to compare.
#
# "file"       : ROOT histogram file produced by makeLFVHistos.py
# "label"      : text shown in the legend
# "color"      : ROOT color
# "style"      : ROOT line style
# "xsec_pb"    : sample cross section in pb
# "n_original" : number of generated MC events before selections
#
# If NORMALIZE = False, histograms are scaled to expected events / fb^{-1}:
#
#     scale = xsec_pb * 1000 / n_original
#
SAMPLES = [

    {
        "name":  "ZTauTau",
        "label": "Z #rightarrow #tau#tau",
        "file":  "histos_ZTauTau_emu_long_sel.root",
        "color": ROOT.kBlue + 1,
        "style": 2,
        "xsec_pb": 1476.58,       # in pb
        "n_original": 1430000,   # <-- user input
    },

    {
        "name":  "ZMuTau",
        "label": "Z #rightarrow #mu#tau",
        "file":  "histos_ZMuTau_emu_long_sel.root",
        "color": ROOT.kRed + 1,
        "style": 2,
        "xsec_pb": 1476.58/0.0336*1e-5,       # <-- user input
        "n_original": 2e6,   # <-- user input
    },
#

#    {
#        "name":  "ZTauTau",
#        "label": "Z #rightarrow #tau#tau",
#        "file":  "histos_ZTauTau_mutauh_long.root",
#        "color": ROOT.kBlue + 1,
#        "style": 2,
#        "xsec_pb": 1476.58,       # in pb 
#        "n_original": 1430000,   # <-- user input
#    },
#    {
#        "name":  "ZTauTau",
#        "label": "Z #rightarrow #tau#tau",
#        "file":  "histos_ZTauTau_mutauh_long_sel.root",
#        "color": ROOT.kGreen + 1,
#        "style": 1,
#        "xsec_pb": 1476.58,       # in pb
#        "n_original": 1430000,   # <-- user input
#    },

#    {
#        "name":  "ZMuTau",
#        "label": "Z #rightarrow #mu#tau",
#        "file":  "histos_ZMuTau_mutauh_long.root",
#        "color": ROOT.kRed + 1,
#        "style": 2,
#        "xsec_pb": 1476.58/0.0336*1e-5,       # <-- user input
#        "n_original": 2e6,   # <-- user input
#    },
#    {
#        "name":  "ZMuTau",
#        "label": "Z #rightarrow #mu#tau",
#        "file":  "histos_ZMuTau_mutauh_long_sel.root",
#        "color": ROOT.kMagenta + 1,
#        "style": 1,
#        "xsec_pb": 1476.58/0.0336*1e-5,       # <-- user input
#        "n_original": 2e6,   # <-- user input
#    },



    # Example:
    # {
    #     "name":  "HMuTau",
    #     "label": "H #rightarrow #mu#tau",
    #     "file":  "hist_HMuTau.root",
    #     "color": ROOT.kGreen + 2,
    #     "style": 1,
    #     "xsec_pb": 1.0,
    #     "n_original": 2000,
    # },
]

# Save each plot in both formats.
SAVE_PNG = True
SAVE_PDF = True

# Optional global log scale for the individual plots and summary pads.
LOGY = False

# ROOT drawing configuration.
LINE_WIDTH = 3
LEGEND_X1 = 0.60
LEGEND_Y1 = 0.72
LEGEND_X2 = 0.88
LEGEND_Y2 = 0.88

# =============================================================================


CHANNEL_TITLES = {
    "mutauh":   "#mu-#tau_{h}",
    "etauh":    "e-#tau_{h}",
    "emu":      "e-#mu",
    "tauhtauh": "#tau_{h}-#tau_{h}",
    "mumu":     "#mu-#mu",
    "ee":       "e-e",
}


# Summary canvases.
# The histogram names are those written by makeLFVHistos.py.
SUMMARY_GROUPS = {
    "obj1_kinematics": {
        "title": "Selected object 1 kinematics",
        "histograms": [
            "obj1_p",
            "obj1_pt",
            "obj1_e",
            "obj1_theta",
            "obj1_costheta",
            "obj1_phi",
        ],
        "cols": 3,
        "rows": 2,
    },

    "obj2_kinematics": {
        "title": "Selected object 2 kinematics",
        "histograms": [
            "obj2_p",
            "obj2_pt",
            "obj2_e",
            "obj2_theta",
            "obj2_costheta",
            "obj2_phi",
        ],
        "cols": 3,
        "rows": 2,
    },

    "neutrinos": {
        "title": "Summed status-1 neutrino system",
        "histograms": [
            "nuSum_p",
            "nuSum_pt",
            "nuSum_e",
            "nuSum_mass",
            "nuSum_theta",
            "nuSum_costheta",
            "nuSum_phi",
            "nNeutrinos",
        ],
        "cols": 4,
        "rows": 2,
    },

    "event_multiplicities": {
        "title": "Generator-level object multiplicities",
        "histograms": [
            "nMuons",
            "nElectrons",
            "nTaus",
            "nPhotons",
            "nNeutrinos",
        ],
        "cols": 3,
        "rows": 2,
    },

    "angles": {
        "title": "Angular relations between selected objects and neutrino system",
        "histograms": [
            "obj1_obj2_deltaPhi",
            "obj1_obj2_deltaTheta",
            "obj1_obj2_dRThetaPhi",
            "obj1_obj2_angle3D",
            "obj1_nuSum_deltaPhi",
            "obj1_nuSum_deltaTheta",
            "obj1_nuSum_dRThetaPhi",
            "obj2_nuSum_deltaPhi",
            "obj2_nuSum_deltaTheta",
            "obj2_nuSum_dRThetaPhi",
        ],
        "cols": 4,
        "rows": 3,
    },

    "visible_mother": {
        "title": "Mother reconstructed from visible selected four-vectors",
        "histograms": [
            "visMother_mass",
            "visMother_p",
            "visMother_pt",
            "visMother_e",
            "visMother_theta",
            "visMother_costheta",
            "visMother_phi",
        ],
        "cols": 4,
        "rows": 2,
    },

    "constrained_mass_ltauh": {
      "title": "Tau-mass constrained reconstruction",
      "histograms": [
        "visMother_mass", 
        "constrainedMass_ltauh",
        "constrainedNuPz_ltauh",
        "tauConstraintPhysical_ltauh",
      ],
      "cols": 3,
      "rows": 1,
    },

    "constrained_mass_ll": {
      "title": "Tau-mass constrained reconstruction",
      "histograms": [
        "visMother_mass",
        "collinearMass_ll",
        "coll_alpha_ll",
        "nuPerp_ll",
      ],
      "cols": 2,
      "rows": 2,
    },



    "real_mother": {
        "title": "Truth-level mother reconstructed using true parent tau four-vectors",
        "histograms": [
            "realMother_mass",
            "realMother_p",
            "realMother_pt",
            "realMother_e",
            "realMother_theta",
            "realMother_costheta",
            "realMother_phi",
        ],
        "cols": 4,
        "rows": 2,
    },

    "gen_mother": {
        "title": "Generated Z/H before decay",
        "histograms": [
            "genMother_mass",
            "genMother_p",
            "genMother_pt",
            "genMother_e",
            "genMother_theta",
            "genMother_costheta",
            "genMother_phi",
            "genMother_ID",
        ],
        "cols": 4,
        "rows": 2,
    },

    "tau_information": {
        "title": "Selected-object tau information",
        "histograms": [
            "obj1FromTau",
            "obj2FromTau",
            "obj1_TauID",
            "obj2_TauID",
        ],
        "cols": 2,
        "rows": 2,
    },

    "mass_comparison": {
        "title": "Visible, true-parent and generated mother mass information",
        "histograms": [
            "visMother_mass",
            "realMother_mass",
            "genMother_mass",
            "visMinusGenMass",
            "realMinusGenMass",
        ],
        "cols": 3,
        "rows": 2,
    },
}


def set_root_style():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetPadLeftMargin(0.14)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadTopMargin(0.10)
    ROOT.gStyle.SetPadRightMargin(0.05)


def open_samples():
    """Open all configured sample files."""
    opened = []

    for sample in SAMPLES:
        f = ROOT.TFile.Open(sample["file"], "READ")

        if not f or f.IsZombie():
            raise RuntimeError(
                "Could not open sample '%s': %s" %
                (sample["name"], sample["file"])
            )

        opened.append({
            **sample,
            "rootfile": f,
        })

    return opened


def list_histograms(samples):
    """Return the union of all TH1 histogram names found in the sample files."""
    names = set()

    for sample in samples:
        rootfile = sample["rootfile"]

        for key in rootfile.GetListOfKeys():
            obj = rootfile.Get(key.GetName())

            if obj and obj.InheritsFrom("TH1"):
                names.add(key.GetName())

    return sorted(names)


def get_histogram(sample, hist_name):
    """
    Get a detached clone of one histogram and apply the sample style.
    Returns None if this sample does not contain the histogram.
    """
    original = sample["rootfile"].Get(hist_name)

    if not original or not original.InheritsFrom("TH1"):
        return None

    h = original.Clone("%s__%s" % (hist_name, sample["name"]))
    h.SetDirectory(0)

    if h.GetEntries() == 0:
      return None

    h.SetLineColor(sample["color"])
    h.SetMarkerColor(sample["color"])
    h.SetLineStyle(sample["style"])
    h.SetLineWidth(LINE_WIDTH)

    if NORMALIZE:
        integral = h.Integral(0, h.GetNbinsX() + 1)
        if integral > 0:
            h.Scale(1.0 / integral)
        h.GetYaxis().SetTitle("Normalized events")
    else:
        xsec_pb = float(sample["xsec_pb"])
        n_original = float(sample["n_original"])

        if n_original <= 0:
            raise RuntimeError(
                "Sample '%s' has invalid n_original = %s" %
                (sample["name"], sample["n_original"])
            )

        # 1 pb = 1000 fb.
        # Expected events for 1 fb^{-1}:
        #   N = sigma[pb] * 1000[fb/pb]
        scale = xsec_pb * 1000.0 / n_original
        h.Scale(scale)

        h.GetYaxis().SetTitle("Events / fb^{-1}")

    return h


def make_legend():
    legend = ROOT.TLegend(
        LEGEND_X1, LEGEND_Y1,
        LEGEND_X2, LEGEND_Y2
    )
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    return legend


def draw_histogram_comparison(samples, hist_name, pad=None, legend=True):
    """
    Draw one histogram for all samples on the current pad.
    Returns objects that must be kept alive until the canvas is saved.
    """
    if pad:
        pad.cd()

    if LOGY:
        ROOT.gPad.SetLogy(True)

    histograms = []
    used_samples = []

    for sample in samples:
        h = get_histogram(sample, hist_name)

        if h is None:
            continue

        histograms.append(h)
        used_samples.append(sample)

    if len(histograms) == 0:
        return []

    ymax = max(h.GetMaximum() for h in histograms)

    if LOGY:
        positive = [
            h.GetMinimum(0.0)
            for h in histograms
            if h.GetMinimum(0.0) > 0
        ]
        ymin = min(positive) * 0.5 if positive else 1e-6
        for h in histograms:
            h.SetMinimum(ymin)
            h.SetMaximum(max(ymax * 20.0, ymin * 100.0))
    else:
        for h in histograms:
            h.SetMinimum(0.0)
            h.SetMaximum(ymax * 1.30 if ymax > 0 else 1.0)

    # Draw first histogram with axes, then overlay the others.
    histograms[0].Draw("HIST")

    for h in histograms[1:]:
        h.Draw("HIST SAME")

    keepalive = list(histograms)

    if legend:
        leg = make_legend()

        for h, sample in zip(histograms, used_samples):
            leg.AddEntry(h, sample["label"], "l")

        leg.Draw()
        keepalive.append(leg)

    ROOT.gPad.RedrawAxis()
    return keepalive


def save_canvas(canvas, output_base):
    if SAVE_PNG:
        canvas.SaveAs(output_base + ".png")

    if SAVE_PDF:
        canvas.SaveAs(output_base + ".pdf")


def make_individual_plots(samples, hist_names, output_dir, channel_title):
    individual_dir = os.path.join(output_dir, "individual")
    os.makedirs(individual_dir, exist_ok=True)

    for hist_name in hist_names:
        canvas = ROOT.TCanvas(
            "c_%s" % hist_name,
            "%s: %s" % (channel_title, hist_name),
            850, 700
        )

        keepalive = draw_histogram_comparison(
            samples,
            hist_name,
            canvas,
            legend=True
        )

        if len(keepalive) == 0:
            canvas.Close()
            continue

        # Put the channel clearly on every individual plot.
        latex = ROOT.TLatex()
        latex.SetNDC(True)
        latex.SetTextSize(0.035)
        latex.SetTextAlign(13)
        latex.DrawLatex(
            0.15, 0.94,
            "Channel: %s" % channel_title
        )
        keepalive.append(latex)

        save_canvas(
            canvas,
            os.path.join(individual_dir, hist_name)
        )

        canvas.Close()


def make_summary_canvas(samples, group_name, group, output_dir, channel_title):
    """
    Build one canvas with a separate pad for each variable.
    A dedicated header pad carries the summary title.
    """

    available = []

    for hist_name in group["histograms"]:

      has_entries = False

      for sample in samples:
        h = sample["rootfile"].Get(hist_name)

        if h and h.InheritsFrom("TH1") and h.GetEntries() > 0:
            has_entries = True
            break

      if has_entries:
        available.append(hist_name)

    if len(available) == 0:
      print("Skipping empty summary:", group_name)
      return

    cols = group["cols"]
    rows = group["rows"]

    width = 420 * cols
    height = 330 * rows + 70

    canvas = ROOT.TCanvas(
        "summary_%s" % group_name,
        "%s - %s" % (channel_title, group["title"]),
        width, height
    )

    # Header with a clear overall title.
    header = ROOT.TPad(
        "header_%s" % group_name,
        "",
        0.0, 0.94, 1.0, 1.0
    )
    header.SetFillStyle(0)
    header.SetMargin(0, 0, 0, 0)
    header.Draw()
    header.cd()

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextAlign(22)
    latex.SetTextSize(0.45)
    latex.DrawLatex(
        0.5, 0.50,
        "%s  |  %s" % (channel_title, group["title"])
    )

    # Main pad containing the variable plots.
    canvas.cd()
    body = ROOT.TPad(
        "body_%s" % group_name,
        "",
        0.0, 0.0, 1.0, 0.94
    )
    body.SetFillStyle(0)
    body.Draw()
    body.cd()
    body.Divide(cols, rows, 0.001, 0.001)

    keepalive = [header, body, latex]

    for i, hist_name in enumerate(available):
        body.cd(i + 1)

        ROOT.gPad.SetTicks(1, 1)

        drawn = draw_histogram_comparison(
            samples,
            hist_name,
            ROOT.gPad,
            legend=True
        )

        keepalive.extend(drawn)

    canvas.cd()
    canvas.Update()

    save_canvas(
        canvas,
        os.path.join(output_dir, "summary_%s" % group_name)
    )

    canvas.Close()


def main():
    parser = argparse.ArgumentParser(
        description="Compare LFV histogram files and make individual and summary plots."
    )

    parser.add_argument(
        "-c", "--channel",
        required=True,
        choices=sorted(CHANNEL_TITLES.keys()),
        help="Selected analysis channel"
    )

    parser.add_argument(
        "-v", "--version",
        default="V1",
        help="Version tag used in the output directory, e.g. V1, V2"
    )

    parser.add_argument(
        "-d", "--outdir",
        default="",
        help="Optional output directory. Default: plots_<channel>_<version>"
    )

    args = parser.parse_args()

    set_root_style()

    channel_title = CHANNEL_TITLES[args.channel]

    output_dir = args.outdir
    if not output_dir:
        output_dir = "plots_%s_%s" % (args.channel, args.version)

    os.makedirs(output_dir, exist_ok=True)

    samples = open_samples()

    print("Samples:")
    for sample in samples:
        if NORMALIZE:
            scale_text = "unit-area normalization"
        else:
            scale = sample["xsec_pb"] * 1000.0 / float(sample["n_original"])
            scale_text = "scale = %.6g per generated event" % scale

        print(
            "  %-15s  %-30s  %s  |  xsec = %.6g pb, Ngen = %d, %s" %
            (
                sample["name"],
                sample["label"],
                sample["file"],
                sample["xsec_pb"],
                sample["n_original"],
                scale_text,
            )
        )

    print("Channel :", args.channel)
    print("Version :", args.version)
    print("Output  :", output_dir)
    print("Normalize:", NORMALIZE)

    hist_names = list_histograms(samples)

    print("Found %d histogram variables" % len(hist_names))

    # -------------------------------------------------------------------------
    # One plot per histogram variable.
    # -------------------------------------------------------------------------
    make_individual_plots(
        samples,
        hist_names,
        output_dir,
        channel_title
    )

    # -------------------------------------------------------------------------
    # Grouped summary canvases.
    # -------------------------------------------------------------------------
    for group_name, group in SUMMARY_GROUPS.items():
        print("Making summary:", group_name)

        make_summary_canvas(
            samples,
            group_name,
            group,
            output_dir,
            channel_title
        )

    for sample in samples:
        sample["rootfile"].Close()

    print("Done. Plots stored in:", output_dir)


if __name__ == "__main__":
    main()

