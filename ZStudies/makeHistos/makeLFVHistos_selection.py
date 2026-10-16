#!/usr/bin/env python3

import argparse
import glob
import math
import ROOT

import reconstructZMass

###### HELP
##### to run:
###### "python makeLFVHistos_selection  -f INPUTFILE.root -o OUTPUTFILE.root -c CHANNEL (eg: 'mutauh') --pmax 45 (max P for plots) --mmax 160 (max M) --selection "

ROOT.TH1.AddDirectory(False)

CHANNELS = {
    "mutauh":   ("MuTauh",   "Muon",     "Tau"),
    "etauh":    ("ETauh",    "Electron", "Tau"),
    "emu":      ("EMu",      "Electron", "Muon"),
    "tauhtauh": ("TauhTauh", "Tau",      "Tau"),
    "mumu":     ("MuMu",     "Muon",     "Muon"),
    "ee":       ("EE",       "Electron", "Electron"),
}


# =============================================================================
# OPTIONAL EVENT SELECTION
# =============================================================================
#
# Applied only when running with:
#
#     --selection
#
# The normal histograms are filled AFTER these cuts when the selection is on.
# A small set of selection_*_before/after histograms and selection_cutflow
# are also written for diagnostics.
#
# For MuTauh:
#     object 1 = muon
#     object 2 = visible hadronic tau
#
# Set a cut to None to disable it.
#
SELECTIONS = {
    "mutauh": {
        "obj1_p_min": 42.0,
        "obj1_p_max": None,
        "obj2_p_min": None,
        "obj2_p_max": None,

        "obj1_nu_dphi_min": 2.5,
        "obj1_nu_dphi_max": None,
        "obj2_nu_dphi_min": None,
        "obj2_nu_dphi_max": 0.5,

        "vis_mass_min": None,
        "vis_mass_max": None,

        # True  -> require a real solution from the constrained-mass reconstruction
        # False -> do not apply this cut
        "hasRealSolution": True,

    },

    # Add channel-specific selections here when needed.
    "etauh": {},
    "emu": {
        "obj1_p_min": None,
        "obj1_p_max": None,
        "obj2_p_min": 40,
        "obj2_p_max": None,

        "obj1_nu_dphi_min": None,
        "obj1_nu_dphi_max": 0.5,
        "obj2_nu_dphi_min": 2.5,
        "obj2_nu_dphi_max": None,

        "vis_mass_min": None,
        "vis_mass_max": None,

        # True  -> require a real solution from the constrained-mass reconstruction
        # False -> do not apply this cut
        "hasRealSolution": True,
    },
    "tauhtauh": {},
    "mumu": {},
    "ee": {},
}


def make_p4(px, py, pz, e):
    p4 = ROOT.TLorentzVector()
    p4.SetPxPyPzE(px, py, pz, e)
    return p4


def object_p4(tree, obj, idx):
    """Visible four-vector of the selected object."""
    prefix = "Tau_vis" if obj == "Tau" else obj
    return make_p4(
        getattr(tree, prefix + "_px")[idx],
        getattr(tree, prefix + "_py")[idx],
        getattr(tree, prefix + "_pz")[idx],
        getattr(tree, prefix + "_e")[idx],
    )


def pair_p4(tree, channel, kind):
    """Four-vector of the selected pair. kind = vis or real."""
    prefix = channel + "_" + kind
    return make_p4(
        getattr(tree, prefix + "_px"),
        getattr(tree, prefix + "_py"),
        getattr(tree, prefix + "_pz"),
        getattr(tree, prefix + "_e"),
    )


def neutrino_sum_p4(tree):
    """Four-vector made from the sum of all status-1 neutrinos in the event."""
    return make_p4(
        tree.NeutrinoSum_px,
        tree.NeutrinoSum_py,
        tree.NeutrinoSum_pz,
        tree.NeutrinoSum_e,
    )


def delta_phi(p41, p42):
    return abs(ROOT.TVector2.Phi_mpi_pi(p41.Phi() - p42.Phi()))


def delta_theta(p41, p42):
    return abs(p41.Theta() - p42.Theta())


def delta_r_theta_phi(p41, p42):
    """
    Angular distance using theta and phi, NOT eta:
        sqrt(DeltaTheta^2 + DeltaPhi^2)
    """
    dphi = delta_phi(p41, p42)
    dtheta = delta_theta(p41, p42)
    return math.sqrt(dtheta * dtheta + dphi * dphi)


def book1(name, title, nbins, xmin, xmax):
    return ROOT.TH1F(name, title, nbins, xmin, xmax)


def book_object_histos(h, label, description, pmax):
    h[label + "_p"] = book1(
        label + "_p", "%s: momentum;p [GeV];Events" % description, 100, 0, pmax)
    h[label + "_pt"] = book1(
        label + "_pt", "%s: transverse momentum;p_{T} [GeV];Events" % description, 100, 0, pmax)
    h[label + "_e"] = book1(
        label + "_e", "%s: energy;E [GeV];Events" % description, 100, 0, pmax)
    h[label + "_theta"] = book1(
        label + "_theta", "%s: polar angle;#theta [rad];Events" % description, 64, 0, math.pi)
    h[label + "_costheta"] = book1(
        label + "_costheta", "%s: cosine of polar angle;cos(#theta);Events" % description, 60, -1, 1)
    h[label + "_phi"] = book1(
        label + "_phi", "%s: azimuthal angle;#phi [rad];Events" % description, 64, -math.pi, math.pi)
    h[label + "_poverebeam"] = book1(
        label + "_poverebeam", "%s: momentum fraction ;x=p/p_{beam};Events" % description, 1000, 0, 2)



def fill_object_histos(h, label, p4,ebeam):
    h[label + "_p"].Fill(p4.P())
    h[label + "_pt"].Fill(p4.Pt())
    h[label + "_e"].Fill(p4.E())
    h[label + "_theta"].Fill(p4.Theta())
    h[label + "_costheta"].Fill(math.cos(p4.Theta()))
    h[label + "_phi"].Fill(p4.Phi())
    h[label + "_poverebeam"].Fill(p4.P()/ebeam)


def book_mother_histos(h, label, description, pmax, mmax):
    h[label + "_mass"] = book1(
        label + "_mass", "%s: invariant mass;m [GeV];Events" % description, 120, 0, mmax)
    h[label + "_p"] = book1(
        label + "_p", "%s: momentum;p [GeV];Events" % description, 100, 0, pmax)
    h[label + "_pt"] = book1(
        label + "_pt", "%s: transverse momentum;p_{T} [GeV];Events" % description, 100, 0, pmax)
    h[label + "_e"] = book1(
        label + "_e", "%s: energy;E [GeV];Events" % description, 100, 0, pmax)
    h[label + "_theta"] = book1(
        label + "_theta", "%s: polar angle;#theta [rad];Events" % description, 64, 0, math.pi)
    h[label + "_costheta"] = book1(
        label + "_costheta", "%s: cosine of polar angle;cos(#theta);Events" % description, 60, -1, 1)
    h[label + "_phi"] = book1(
        label + "_phi", "%s: azimuthal angle;#phi [rad];Events" % description, 64, -math.pi, math.pi)


def fill_mother_histos(h, label, p4):
    h[label + "_mass"].Fill(p4.M())
    h[label + "_p"].Fill(p4.P())
    h[label + "_pt"].Fill(p4.Pt())
    h[label + "_e"].Fill(p4.E())
    h[label + "_theta"].Fill(p4.Theta())
    h[label + "_costheta"].Fill(math.cos(p4.Theta()))
    h[label + "_phi"].Fill(p4.Phi())


def book_pair_angles(h, prefix, description):
    h[prefix + "_deltaPhi"] = book1(
        prefix + "_deltaPhi",
        "%s: |#Delta#phi|;|#Delta#phi| [rad];Events" % description,
        64, 0, math.pi)

    h[prefix + "_deltaTheta"] = book1(
        prefix + "_deltaTheta",
        "%s: |#Delta#theta|;|#Delta#theta| [rad];Events" % description,
        64, 0, math.pi)

    h[prefix + "_dRThetaPhi"] = book1(
        prefix + "_dRThetaPhi",
        "%s: #DeltaR_{#theta#phi};#DeltaR_{#theta#phi};Events" % description,
        80, 0, math.sqrt(2.0) * math.pi)


def fill_pair_angles(h, prefix, p41, p42):
    h[prefix + "_deltaPhi"].Fill(delta_phi(p41, p42))
    h[prefix + "_deltaTheta"].Fill(delta_theta(p41, p42))
    h[prefix + "_dRThetaPhi"].Fill(delta_r_theta_phi(p41, p42))



def build_selection_cuts(config, obj1_type, obj2_type):
    """
    Convert the user configuration into an ordered list of active cuts.

    Each entry is:
        (configuration key, threshold/value, human-readable label)

    For hasRealSolution:
        True  -> require a real constrained-mass solution
        False -> do not apply this cut
    """
    order = [
        "obj1_p_min", "obj1_p_max",
        "obj2_p_min", "obj2_p_max",
        "obj1_nu_dphi_min", "obj1_nu_dphi_max",
        "obj2_nu_dphi_min", "obj2_nu_dphi_max",
        "vis_mass_min", "vis_mass_max",
    ]

    cuts = []

    for key in order:
        if key not in config or config[key] is None:
            continue

        value = float(config[key])

        labels = {
            "obj1_p_min": "p(%s) > %.3g GeV" % (obj1_type, value),
            "obj1_p_max": "p(%s) < %.3g GeV" % (obj1_type, value),
            "obj2_p_min": "p(%s) > %.3g GeV" % (obj2_type, value),
            "obj2_p_max": "p(%s) < %.3g GeV" % (obj2_type, value),

            "obj1_nu_dphi_min":
                "|DeltaPhi(%s,nu sum)| > %.3g" % (obj1_type, value),
            "obj1_nu_dphi_max":
                "|DeltaPhi(%s,nu sum)| < %.3g" % (obj1_type, value),

            "obj2_nu_dphi_min":
                "|DeltaPhi(%s,nu sum)| > %.3g" % (obj2_type, value),
            "obj2_nu_dphi_max":
                "|DeltaPhi(%s,nu sum)| < %.3g" % (obj2_type, value),

            "vis_mass_min": "m(vis) > %.3g GeV" % value,
            "vis_mass_max": "m(vis) < %.3g GeV" % value,
        }

        cuts.append((key, value, labels[key]))

    # Boolean switch: False means the cut is disabled.
    if config.get("hasRealSolution", False):
        cuts.append((
            "hasRealSolution",
            True,
            "Constrained mass: real solution"
        ))

    return cuts

def evaluate_selection_cut(
    key, threshold, p41, p42, nu_sum, vis, hasRealSolution
):
    """Evaluate one configured selection cut."""

    if key == "obj1_p_min":
        return p41.P() > threshold
    if key == "obj1_p_max":
        return p41.P() < threshold

    if key == "obj2_p_min":
        return p42.P() > threshold
    if key == "obj2_p_max":
        return p42.P() < threshold

    if key == "obj1_nu_dphi_min":
        return nu_sum.P() > 0 and delta_phi(p41, nu_sum) > threshold
    if key == "obj1_nu_dphi_max":
        return nu_sum.P() > 0 and delta_phi(p41, nu_sum) < threshold

    if key == "obj2_nu_dphi_min":
        return nu_sum.P() > 0 and delta_phi(p42, nu_sum) > threshold
    if key == "obj2_nu_dphi_max":
        return nu_sum.P() > 0 and delta_phi(p42, nu_sum) < threshold

    if key == "vis_mass_min":
        return vis.M() > threshold
    if key == "vis_mass_max":
        return vis.M() < threshold

    if key == "hasRealSolution":
        return bool(hasRealSolution)

    raise RuntimeError("Unknown selection cut: %s" % key)

def book_selection_diagnostics(h, pmax, mmax, active_cuts):
    """Before/after plots and cumulative cutflow for the event selection."""

    h["selection_obj1_p_before"] = book1(
        "selection_obj1_p_before",
        "Before selection: object 1 momentum;p(object 1) [GeV];Events",
        100, 0, pmax)
    h["selection_obj1_p_after"] = book1(
        "selection_obj1_p_after",
        "After selection: object 1 momentum;p(object 1) [GeV];Events",
        100, 0, pmax)

    h["selection_obj2_p_before"] = book1(
        "selection_obj2_p_before",
        "Before selection: object 2 momentum;p(object 2) [GeV];Events",
        100, 0, pmax)
    h["selection_obj2_p_after"] = book1(
        "selection_obj2_p_after",
        "After selection: object 2 momentum;p(object 2) [GeV];Events",
        100, 0, pmax)

    h["selection_obj1_nu_deltaPhi_before"] = book1(
        "selection_obj1_nu_deltaPhi_before",
        "Before selection: |#Delta#phi|(object 1, neutrino sum);|#Delta#phi| [rad];Events",
        64, 0, math.pi)
    h["selection_obj1_nu_deltaPhi_after"] = book1(
        "selection_obj1_nu_deltaPhi_after",
        "After selection: |#Delta#phi|(object 1, neutrino sum);|#Delta#phi| [rad];Events",
        64, 0, math.pi)

    h["selection_obj2_nu_deltaPhi_before"] = book1(
        "selection_obj2_nu_deltaPhi_before",
        "Before selection: |#Delta#phi|(object 2, neutrino sum);|#Delta#phi| [rad];Events",
        64, 0, math.pi)
    h["selection_obj2_nu_deltaPhi_after"] = book1(
        "selection_obj2_nu_deltaPhi_after",
        "After selection: |#Delta#phi|(object 2, neutrino sum);|#Delta#phi| [rad];Events",
        64, 0, math.pi)

    h["selection_visMass_before"] = book1(
        "selection_visMass_before",
        "Before selection: visible mother mass;m_{vis} [GeV];Events",
        120, 0, mmax)
    h["selection_visMass_after"] = book1(
        "selection_visMass_after",
        "After selection: visible mother mass;m_{vis} [GeV];Events",
        120, 0, mmax)

    # Constrained and collinear masses before/after the full selection.
    # Keep the same useful Z-region ranges as the main histograms.
    h["selection_constrainedMass_before"] = book1(
        "selection_constrainedMass_before",
        "Before selection: constrained mass;m_{constr} [GeV];Events",
        240, 70, 100)
    h["selection_constrainedMass_after"] = book1(
        "selection_constrainedMass_after",
        "After selection: constrained mass;m_{constr} [GeV];Events",
        240, 70, 100)

    h["selection_collinearMass_before"] = book1(
        "selection_collinearMass_before",
        "Before selection: collinear mass;m_{coll} [GeV];Events",
        200, 60, 120)
    h["selection_collinearMass_after"] = book1(
        "selection_collinearMass_after",
        "After selection: collinear mass;m_{coll} [GeV];Events",
        200, 60, 120)

    # Cumulative cutflow using NORMAL ROOT bins.
    #
    # bin 1                 = all processed events
    # bin 2                 = valid channel pairing
    # bins 3 ... N+2        = cumulative active cuts
    # final bin N+3         = final selection
    #
    # Do not use ROOT bin 0: it is the underflow bin and is not displayed.
    nbins = len(active_cuts) + 3

    cutflow = ROOT.TH1F(
        "selection_cutflow",
        "Cumulative event selection;Selection step;Events",
        nbins, 0.5, nbins + 0.5
    )

    cutflow.GetXaxis().SetBinLabel(1, "All events")
    cutflow.GetXaxis().SetBinLabel(2, "Channel Selection")

    for i, (_, _, label) in enumerate(active_cuts):
        cutflow.GetXaxis().SetBinLabel(i + 3, label)

    cutflow.GetXaxis().SetBinLabel(nbins, "Final selection")

    # Long cut labels are easier to read vertically when drawn directly in ROOT.
    cutflow.GetXaxis().LabelsOption("v")
    cutflow.GetXaxis().SetLabelSize(0.030)

    h["selection_cutflow"] = cutflow

def fill_selection_diagnostics(
    h, suffix, p41, p42, nu_sum, vis, mConstrained, mCollinear
):
    h["selection_obj1_p_" + suffix].Fill(p41.P())
    h["selection_obj2_p_" + suffix].Fill(p42.P())
    h["selection_visMass_" + suffix].Fill(vis.M())

    if mConstrained > -900:
        h["selection_constrainedMass_" + suffix].Fill(mConstrained)

    if mCollinear > -900:
        h["selection_collinearMass_" + suffix].Fill(mCollinear)

    if nu_sum.P() > 0:
        h["selection_obj1_nu_deltaPhi_" + suffix].Fill(delta_phi(p41, nu_sum))
        h["selection_obj2_nu_deltaPhi_" + suffix].Fill(delta_phi(p42, nu_sum))

def passes_selection(
    h, active_cuts, p41, p42, nu_sum, vis, hasRealSolution
):
    """
    Apply cuts sequentially and fill a cumulative cutflow.

    Bin 1 (All events) and bin 2 (Channel Selection) are filled
    in the event loop. Here we fill the active-cut bins and the
    final-selection bin.
    """

    for i, (key, threshold, _) in enumerate(active_cuts):
        if not evaluate_selection_cut(
            key, threshold, p41, p42, nu_sum, vis, hasRealSolution
        ):
            return False

        # Cumulative count after this cut.
        h["selection_cutflow"].Fill(i + 3)

    # Final selection is a separate visible bin, even though its content
    # equals the final active-cut count.
    h["selection_cutflow"].Fill(len(active_cuts) + 3)

    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--input", required=True,
                        help="Input ROOT file, glob, or comma-separated list")
    parser.add_argument("-o", "--output", default="",
                        help="Output ROOT histogram file")
    parser.add_argument("-c", "--channel", required=True,
                        choices=sorted(CHANNELS.keys()))
    parser.add_argument("--pmax", type=float, default=150.0)
    parser.add_argument("--mmax", type=float, default=160.0)
    parser.add_argument(
        "--selection",
        action="store_true",
        help="Apply the channel selection configured in SELECTIONS"
    )
    args = parser.parse_args()

    channel, obj1_type, obj2_type = CHANNELS[args.channel]

    selection_config = SELECTIONS.get(args.channel, {})
    active_cuts = build_selection_cuts(
        selection_config, obj1_type, obj2_type
    )

    if args.selection and len(active_cuts) == 0:
        raise RuntimeError(
            "Selection requested, but no active cuts are configured for channel '%s'"
            % args.channel
        )

    output = args.output if args.output else "hist_%s.root" % args.channel

    chain = ROOT.TChain("outtree")

    inputs = []
    for item in args.input.split(","):
        matches = glob.glob(item)
        inputs.extend(matches if matches else [item])

    for filename in inputs:
#        print("Adding", filename)
        chain.Add(filename)

    if chain.GetEntries() == 0:
        raise RuntimeError("No entries found in input tree")

    print ("I will run over ",chain.GetEntries()," events")

    required = [
        channel + "_valid",
        channel + "_idx1",
        channel + "_idx2",
        channel + "_vis_px",
        channel + "_vis_py",
        channel + "_vis_pz",
        channel + "_vis_e",
        channel + "_real_px",
        channel + "_real_py",
        channel + "_real_pz",
        channel + "_real_e",
        "NeutrinoSum_px",
        "NeutrinoSum_py",
        "NeutrinoSum_pz",
        "NeutrinoSum_e",
    ]

    for name in required:
        if not chain.GetBranch(name):
            raise RuntimeError("Missing branch: %s" % name)

    h = {}

    obj1_description = "Selected object 1 (%s)" % obj1_type
    obj2_description = "Selected object 2 (%s)" % obj2_type

    # Event multiplicities.
    h["nMuons"] = book1("nMuons", "Number of status-1 muons;N_{#mu};Events", 8, -0.5, 7.5)
    h["nElectrons"] = book1("nElectrons", "Number of status-1 electrons;N_{e};Events", 8, -0.5, 7.5)
    h["nTaus"] = book1("nTaus", "Number of generated taus;N_{#tau};Events", 8, -0.5, 7.5)
    h["nPhotons"] = book1("nPhotons", "Number of status-1 photons;N_{#gamma};Events", 15, -0.5, 14.5)
    if chain.GetBranch("nNeutrinos"):
        h["nNeutrinos"] = book1("nNeutrinos", "Number of status-1 neutrinos;N_{#nu};Events", 12, -0.5, 11.5)

    # Selected visible objects.
    book_object_histos(h, "obj1", obj1_description, args.pmax)
    book_object_histos(h, "obj2", obj2_description, args.pmax)

    # Angular separation between the two selected objects.
    book_pair_angles(h, "obj1_obj2", "Object 1 - object 2")

    # Full 3D angle between the two momentum vectors. This is distinct from DeltaTheta.
    h["obj1_obj2_angle3D"] = book1(
        "obj1_obj2_angle3D",
        "3D angle between selected object 1 and object 2 momentum vectors;3D angle [rad];Events",
        64, 0, math.pi)

    # Sum of all status-1 neutrinos in the event.
    book_object_histos(
        h, "nuSum",
        "Vector sum of all status-1 neutrinos",
        args.pmax)

    # Angular relation of each selected object to the summed-neutrino vector.
    book_pair_angles(h, "obj1_nuSum", "Object 1 - neutrino sum")
    book_pair_angles(h, "obj2_nuSum", "Object 2 - neutrino sum")

    # Reconstructed mother from selected pair.
    book_mother_histos(
        h, "visMother",
        "Mother reconstructed from the visible selected objects",
        args.pmax, args.mmax)

    book_mother_histos(
        h, "realMother",
        "Mother reconstructed using the true parent tau when an object comes from tau decay",
        args.pmax, args.mmax)

    # Generator Z/H before decay.
    book_mother_histos(
        h, "genMother",
        "Generated Z/H before decay",
        args.pmax, args.mmax)

    h["genMother_ID"] = book1(
        "genMother_ID", "PDG ID of generated mother;PDG ID;Events", 61, -30.5, 30.5)

    # Mass comparisons.
    h["visMinusGenMass"] = book1(
        "visMinusGenMass",
        "Visible-pair mass minus generated-mother mass;m_{vis}-m_{gen} [GeV];Events",
        120, -120, 30)

    h["realMinusGenMass"] = book1(
        "realMinusGenMass",
        "True-parent pair mass minus generated-mother mass;m_{real}-m_{gen} [GeV];Events",
        120, -30, 30)

    # Prompt-vs-tau provenance of selected objects.
    if chain.GetBranch(channel + "_obj1FromTau"):
        h["obj1FromTau"] = book1(
            "obj1FromTau", "Selected object 1 originates from tau decay;from #tau (0/1);Events", 2, -0.5, 1.5)
    if chain.GetBranch(channel + "_obj2FromTau"):
        h["obj2FromTau"] = book1(
            "obj2FromTau", "Selected object 2 originates from tau decay;from #tau (0/1);Events", 2, -0.5, 1.5)

    # Tau decay mode for selected tau candidates.
    if obj1_type == "Tau":
        h["obj1_TauID"] = book1(
            "obj1_TauID", "TauID of selected tau object 1;TauID;Events", 40, -20.5, 19.5)
    if obj2_type == "Tau":
        h["obj2_TauID"] = book1(
            "obj2_TauID", "TauID of selected tau object 2;TauID;Events", 40, -20.5, 19.5)

    # mass histograms

    h["constrainedMass"] = book1(
       "constrainedMass",
       "Tau-mass constrained mass); m_{#mu#tau}^{constr} [GeV];Events",
         240, 70, 100
        )

    h["constrainedNuPz"] = book1(
       "constrainedNuPz",
       "Neutrino p_{z} from tau-mass constraint ;p_{z}^{#nu,constr} [GeV];Events",
        120, -60, 60
       )

    h["tauConstraintPhysical"] = book1(
        "tauConstraintPhysical",
        "Tau-mass constraint has real solution;real solution;Events",
         2, -0.5, 1.5
    )
#    h["constrainedMass_ll"] = book1(  # this did not work well
#       "constrainedMass_ll",
#       "Tau-mass constrained visible mass (emu); m_{#mu#tau}^{constr} [GeV];Events",
#         40, 80, 100
#        )

    h["collinearMass"] = book1(
       "collinearMass",
       "Collinear Mass; m_{#mu#tau}^{constr} [GeV];Events",
         200, 60, 120
        )
    h["coll_alpha"] = book1(
       "coll_alpha",
       "Alpha, used in the construction of coll mass; alpha ; Events",
        100, 0, 2
       )
    h["nuPerp"] = book1(
        "nuPerp",
        "Neutrino Perp; nu perp;Events",
         100, -50, 50
    )

    if args.selection:
        book_selection_diagnostics(
            h, args.pmax, args.mmax, active_cuts
        )

        print("Applying selection for channel", args.channel)
        for _, _, label in active_cuts:
            print("   ", label)
    else:
        print("Selection is OFF")

    nchannelvalid = 0
    nselected = 0
    ntotal=0

    for event in chain:
        ntotal+=1
        if ntotal%100000==0 :
           print ("....",ntotal)

        if args.selection:
            h["selection_cutflow"].Fill(1)   # All events

        if getattr(event, channel + "_valid") != 1:
            continue

        idx1 = int(getattr(event, channel + "_idx1"))
        idx2 = int(getattr(event, channel + "_idx2"))

        if idx1 < 0 or idx2 < 0:
            continue

        if args.selection:
            h["selection_cutflow"].Fill(2)   # Channel Selection

        nchannelvalid += 1

        # Build the quantities needed both for the selection and for plotting.
        p41 = object_p4(event, obj1_type, idx1)
        p42 = object_p4(event, obj2_type, idx2)
        nu_sum = neutrino_sum_p4(event)
        vis = pair_p4(event, channel, "vis")
        real = pair_p4(event, channel, "real")

        mConstrained=0
        nuPzConstrained=0
        hasRealSolution=0
        mCollinear=0
        alpha=0
        nuPerp=0


        if (channel=="MuTauh" or channel=="ETauh"):

         mConstrained, nuPzConstrained, hasRealSolution = \
           reconstructZMass.reconstruct_mass_ltauh(
            p41,             # muon in mutauh, or electron in etauh
            p42,             # visible tau
            event.NeutrinoSum_px,
            event.NeutrinoSum_py,
            91.2
           )

         mCollinear, alpha, nuPerp = \
           reconstructZMass.reconstruct_mass_collinear(
              p41,             # prompt muon 
              p42,             # had tau
              event.NeutrinoSum_px,
              event.NeutrinoSum_py
           )
     

        if (channel == "EMu"):

         mConstrained, nuPzConstrained, hasRealSolution = \
          reconstructZMass.reconstruct_mass_ll_constrained(
           p42,   # prompt muon in ZMuTau, change for ZETau. Adaptfor mumu or ee if wanted
           p41,   # electron from tau
           event.NeutrinoSum_px,
           event.NeutrinoSum_py,
           91.2
         )

         mCollinear, alpha, nuPerp = \
           reconstructZMass.reconstruct_mass_collinear(
              p42,             # prompt muon # if you want to change to ETau_mu, change these
              p41,             # electron from tau
              event.NeutrinoSum_px,
              event.NeutrinoSum_py
         )


        # Optional exploratory selection.
        if args.selection:
            fill_selection_diagnostics(
                h, "before",
                p41, p42, nu_sum, vis,
                mConstrained, mCollinear
            )

            if not passes_selection(
                h, active_cuts,
                p41, p42, nu_sum, vis,
                hasRealSolution
            ):
                continue

            fill_selection_diagnostics(
                h, "after",
                p41, p42, nu_sum, vis,
                mConstrained, mCollinear
            )


        h["tauConstraintPhysical"].Fill(hasRealSolution)

        nselected += 1



        # From here onwards, normal histograms contain:
        #   - all valid channel events when --selection is OFF
        #   - only events passing the configured cuts when --selection is ON

        # Multiplicities.
        h["nMuons"].Fill(event.nMuons)
        h["nElectrons"].Fill(event.nElectrons)
        h["nTaus"].Fill(event.nTaus)
        h["nPhotons"].Fill(event.nPhotons)
        if "nNeutrinos" in h:
            h["nNeutrinos"].Fill(event.nNeutrinos)

        # Selected visible objects.
        fill_object_histos(h, "obj1", p41,event.beamE_0)
        fill_object_histos(h, "obj2", p42,event.beamE_0)

        fill_pair_angles(h, "obj1_obj2", p41, p42)

        if p41.P() > 0 and p42.P() > 0:
            h["obj1_obj2_angle3D"].Fill(p41.Vect().Angle(p42.Vect()))

        # Sum of all status-1 neutrinos.
        fill_object_histos(h, "nuSum", nu_sum,event.beamE_0)

        # Angles to the neutrino system only make sense if its momentum is nonzero.
        if nu_sum.P() > 0:
            fill_pair_angles(h, "obj1_nuSum", p41, nu_sum)
            fill_pair_angles(h, "obj2_nuSum", p42, nu_sum)

        # Pair / reconstructed mother.
        fill_mother_histos(h, "visMother", vis)
        fill_mother_histos(h, "realMother", real)

        h["constrainedMass"].Fill(mConstrained)
        h["constrainedNuPz"].Fill(nuPzConstrained) 

        #print(mCollinear)
        h["collinearMass"].Fill(mCollinear)
        h["coll_alpha"].Fill(alpha)
        h["nuPerp"].Fill(nuPerp)


        # Exact generated Z/H before decay.
        if hasattr(event, "GenMother_ID") and abs(event.GenMother_ID) in (23, 25):
            gen = make_p4(
                event.GenMother_px,
                event.GenMother_py,
                event.GenMother_pz,
                event.GenMother_e,
            )

            fill_mother_histos(h, "genMother", gen)
            h["genMother_ID"].Fill(event.GenMother_ID)

            h["visMinusGenMass"].Fill(vis.M() - gen.M())
            h["realMinusGenMass"].Fill(real.M() - gen.M())

        # Prompt / tau-decay provenance.
        if "obj1FromTau" in h:
            h["obj1FromTau"].Fill(getattr(event, channel + "_obj1FromTau"))
        if "obj2FromTau" in h:
            h["obj2FromTau"].Fill(getattr(event, channel + "_obj2FromTau"))

        # TauID for selected tau objects.
        if obj1_type == "Tau":
            h["obj1_TauID"].Fill(event.Tau_TauID[idx1])
        if obj2_type == "Tau":
            h["obj2_TauID"].Fill(event.Tau_TauID[idx2])

    fout = ROOT.TFile(output, "RECREATE")
    for hist in h.values():
        hist.Write()
    fout.Close()

    print("Total events:", ntotal)
    print("Channel-valid events:", nchannelvalid, " ( %.2f",nchannelvalid/ntotal*100,"% of total")

    if args.selection:
        efficiency = (
            float(nselected) / nchannelvalid
            if nchannelvalid > 0 else 0.0
        )
        print(
            "After selection: %d / %d  efficiency = %.6g"
            % (nselected, nchannelvalid, efficiency)
        )
    else:
        print("Selected events:", nselected)

    print("Wrote", output)


if __name__ == "__main__":
    main()


