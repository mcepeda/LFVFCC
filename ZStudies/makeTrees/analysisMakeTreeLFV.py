import sys, os, math
from array import array
from pathlib import Path

import ROOT
from podio import root_io

import reco.tauReco

import argparse
parser = argparse.ArgumentParser(
    description="Build a GEN-level LFV Z/H -> mutau, etau, emu analysis tree",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

parser.add_argument(
    "-f", "--inputfile",
#     default="/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/Py8_ONLYSIM_2M_LFV/out_sim_edm4hep_1.root",
    default="/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/KKMCee_ONLYSIM_10M_90GeV/out_sim_edm4hep_1.root",
)
parser.add_argument("-s", "--sample", default="Py8_ONLYSIM_2M_LFV")
parser.add_argument("-fn", "--filename", default="out_sim_edm4hep")
parser.add_argument("-o", "--outfile", default="test_GEN_LFV_tree.root")
parser.add_argument("-w", "--Winter23", default="False")
parser.add_argument(
    "-c", "--channel",
    default="all",
    choices=["all", "tauhtauh","mutauh", "etauh", "emu","mumu","ee"],
    help="Event selection. 'all' keeps events passing at least one LFV channel.",
)

args = parser.parse_args()
print(vars(args))

fileOutName = args.outfile
winter23 = (args.Winter23 == "True")
print(args.Winter23, winter23)

if winter23:
    import cppyy
    edm4hep = cppyy.gbl.edm4hep
    genparts = "Particle"
else:
    import edm4hep
    genparts = "MCParticles"


# ----------------------------------------------------------------------
# Input files
# ----------------------------------------------------------------------

filenames = []

if args.inputfile != "":
    filenames.append(args.inputfile)
    print("Reading file", args.inputfile)
else:
    path = "/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/"
    dir_path = os.path.join(path, args.sample)
    badfiles = [-1]
    nfiles = 100

    print(dir_path)

    for i in range(1, nfiles + 1):
        if i in badfiles:
            continue

        filename = os.path.join(dir_path, "%s_%d.root" % (args.filename, i))
        if not Path(filename).is_file():
            continue

        root_file = reco.tauReco.open_root_file(filename)
        if not root_file or root_file.IsZombie():
            continue
        root_file.Close()

        filenames.append(filename)

    print("Read %d files" % len(filenames))

if len(filenames) == 0:
    raise RuntimeError("No input files found")

reader = root_io.Reader(filenames)


# ----------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------

def particle_p4(part):
    p4 = ROOT.TLorentzVector()
    p4.SetXYZM(
        part.getMomentum().x,
        part.getMomentum().y,
        part.getMomentum().z,
        part.getMass(),
    )
    return p4


def particle_key(part):
    """Stable-enough key for matching a particle returned through two podio handles."""
    try:
        oid = part.getObjectID()
        collection_id = getattr(oid, "collectionID", -1)
        return (int(collection_id), int(oid.index))
    except Exception:
        return None


def tau_is_hadronic(tau_id):
    # tauReco.visTauGen uses -13/-11 for leptonic tau decays and -999 for
    # an intermediate tau copy. Other valid modes are treated as tau_h.
    return tau_id not in (-13, -11, -999)


def find_tau_ancestor(part):
    """Return the first tau found in the parent chain, if any.

    Tau construction/decay classification itself is *not* redone here; this
    is only a fallback for tagging final-state e/mu after generator FSR.
    """
    seen = set()
    stack = list(part.getParents())

    while stack:
        parent = stack.pop(0)
        key = particle_key(parent)
        if key is not None and key in seen:
            continue
        if key is not None:
            seen.add(key)

        if abs(parent.getPDG()) == 15:
            # Use the same tauReco definition here too: a tau with TauID=-999
            # does not have a decay we are contemplating here 
            try:
                tau_info = reco.tauReco.visTauGen(parent, False)
                if int(tau_info[1]) != -999:
                    return parent
            except Exception:
                pass

        try:
            stack.extend(list(parent.getParents()))
        except Exception:
            pass

    return None


def match_tau_index_from_p4(tau_particle, taus):
    # not clear we need this?º
    """Match a tau particle from the ancestry chain to a tauReco tau record."""
    if tau_particle is None:
        return -1

    p4 = particle_p4(tau_particle)
    q = tau_particle.getCharge()

    best_idx = -1
    best_metric = 1.0e30
    for i, tau in enumerate(taus):
        if abs(tau["charge"] - q) > 0.1:
            continue
        tp4 = tau["true_p4"]
        metric = (
            abs(tp4.Px() - p4.Px())
            + abs(tp4.Py() - p4.Py())
            + abs(tp4.Pz() - p4.Pz())
            + abs(tp4.E() - p4.E())
        )
        if metric < best_metric:
            best_metric = metric
            best_idx = i

    # Exact copies should agree to numerical precision. The loose bound only
    # protects against accidentally matching an unrelated tau.
    if best_metric > 1.0e-3:
        return -1
    return best_idx


def build_gen_taus(mc_particles):
    """Use tauReco.findAllGenTaus/visTauGen as the single tau definition.

    Returned tau tuple from tauReco is:
      [0] visible P4
      [1] TauID
      [2] charge
      [3] true tau P4
      [4] max constituent angle
      [5] n constituents
      [6] constituent dictionary
      [7] helicity
    """
    raw_taus = reco.tauReco.findAllGenTaus(mc_particles, False)
    taus = []

    for i in range(len(raw_taus)):
        tau = raw_taus[i]
        tau_id = int(tau[1])

        if tau_id == -999:
            continue

        record = {
            "vis_p4": tau[0],
            "TauID": tau_id,
            "charge": float(tau[2]),
            "true_p4": tau[3],
            "constituents": tau[6],
        }
        taus.append(record)

    return taus


def build_tau_constituent_map(taus):
    """Map visible tau constituents onto the tau index using tauReco output."""
    constituent_to_tau = {}
    for itau, tau in enumerate(taus):
        const = tau["constituents"]
        for iconst in range(len(const)):
            key = particle_key(const[iconst])
            if key is not None:
                constituent_to_tau[key] = itau
    return constituent_to_tau

def get_photon_info(photon, constituent_to_tau, taus):
    """
    Return direct mother information and tau-decay provenance.

    photon_from_tau = 1 if the photon belongs to the decay chain of one
    of the physical taus returned by tauReco.
    """

    mothers = list(photon.getParents())

    mother_id = -999
    mother_status = -999
    mother_idx = -1

    if len(mothers) > 0:
        mother = mothers[0]

        mother_id = int(mother.getPDG())
        mother_status = int(mother.getGeneratorStatus())

        try:
            mother_idx = int(mother.getObjectID().index)
        except Exception:
            mother_idx = -1

    # First use exactly the same tau constituent information as tauReco.
    key = particle_key(photon)
    parent_tau_idx = (
        constituent_to_tau.get(key, -1)
        if key is not None
        else -1
    )

    # Fallback for cases in which showering/FSR has inserted an
    # intermediate particle between the tau and the final photon.
    if parent_tau_idx < 0:
        tau_ancestor = find_tau_ancestor(photon)
        parent_tau_idx = match_tau_index_from_p4(tau_ancestor, taus)

    from_tau_decay = int(parent_tau_idx >= 0)

    return {
        "motherID": mother_id,
        "motherStatus": mother_status,
        "motherIdx": mother_idx,
        "nMothers": len(mothers),
        "fromTauDecay": from_tau_decay,
        "parentTauIdx": parent_tau_idx,
    }

def sum_status1_neutrinos(mc_particles):
    """Sum px, py, pz (and E for convenience) over all status-1 neutrinos."""
    px = 0.0
    py = 0.0
    pz = 0.0
    e = 0.0

    for part in mc_particles:
        if part.getGeneratorStatus() != 1:
            continue
        if abs(part.getPDG()) not in (12, 14, 16):
            continue

        px += part.getMomentum().x
        py += part.getMomentum().y
        pz += part.getMomentum().z
        e += part.getEnergy()

    return px, py, pz, e


def find_decay_boson(mc_particles):
    """Find the Z/H copy immediately upstream of its lepton/tau decay."""
    best = None
    best_score = -1

    for part in mc_particles:
        if abs(part.getPDG()) not in (23, 25):
            continue

        daughters = list(part.getDaughters())
        if len(daughters) == 0:
            continue

        daughter_ids = [abs(d.getPDG()) for d in daughters]
        n_leptonic_daughters = sum(pid in (11, 13, 15) for pid in daughter_ids)

        # Strongly prefer the boson whose direct daughters are the LFV/tau pair.
        score = 10 * n_leptonic_daughters + len(daughters)
        if score > best_score:
            best = part
            best_score = score

    return best


def original_lepton_p4(lepton, taus):
    """For tau->l use the parent true tau P4; otherwise use the status-1 lepton P4."""
    itau = lepton["parentTauIdx"]
    if itau >= 0:
        return taus[itau]["true_p4"]
    return lepton["p4"]


def choose_best_pair(candidates):
    if len(candidates) == 0:
        return None
    # No mass-based ranking: keep the highest visible momentum-sum candidate.
    return max(candidates, key=lambda x: x["rank"])

def make_tauhtauh_pair(taus):
    candidates = []
    for itau1, tau1 in enumerate(taus):
        for itau2, tau2 in enumerate(taus):
           if itau2>itau1:
            if not tau_is_hadronic(tau1["TauID"]):
                continue
            if not tau_is_hadronic(tau2["TauID"]):
                continue
            if tau1["charge"] * tau2["charge"] >= 0:
                continue

            vis = tau1["vis_p4"] + tau2["vis_p4"]
            real = tau1["true_p4"] + tau2["true_p4"]
            candidates.append({
                "idx1": itau1,
                "idx2": itau2,
                "vis": vis,
                "real": real,
                "fromTau1": 1,
                "fromTau2": 1,
                "rank": tau1["vis_p4"].P() + tau2["vis_p4"].P(),
            })
    return choose_best_pair(candidates)



def make_mutauh_pair(muons, taus):
    candidates = []
    for imu, mu in enumerate(muons):
        for itau, tau in enumerate(taus):
            if not tau_is_hadronic(tau["TauID"]):
                continue
            if mu["charge"] * tau["charge"] >= 0:
                continue

            vis = mu["p4"] + tau["vis_p4"]
            real = original_lepton_p4(mu, taus) + tau["true_p4"]
            candidates.append({
                "idx1": imu,
                "idx2": itau,
                "vis": vis,
                "real": real,
                "fromTau1": mu["isTauDecay"],
                "fromTau2": 1,
                "rank": mu["p4"].P() + tau["vis_p4"].P(),
            })
    return choose_best_pair(candidates)


def make_etauh_pair(electrons, taus):
    candidates = []
    for ie, ele in enumerate(electrons):
        for itau, tau in enumerate(taus):
            if not tau_is_hadronic(tau["TauID"]):
                continue
            if ele["charge"] * tau["charge"] >= 0:
                continue

            vis = ele["p4"] + tau["vis_p4"]
            real = original_lepton_p4(ele, taus) + tau["true_p4"]
            candidates.append({
                "idx1": ie,
                "idx2": itau,
                "vis": vis,
                "real": real,
                "fromTau1": ele["isTauDecay"],
                "fromTau2": 1,
                "rank": ele["p4"].P() + tau["vis_p4"].P(),
            })
    return choose_best_pair(candidates)


def make_emu_pair(electrons, muons, taus):
    candidates = []
    for ie, ele in enumerate(electrons):
        for imu, mu in enumerate(muons):
            if ele["charge"] * mu["charge"] >= 0:
                continue

            # Avoid the pathological case where both visible leptons are mapped
            # to the same tau copy.
            if ele["parentTauIdx"] >= 0 and ele["parentTauIdx"] == mu["parentTauIdx"]:
                continue

            vis = ele["p4"] + mu["p4"]
            real = original_lepton_p4(ele, taus) + original_lepton_p4(mu, taus)
            candidates.append({
                "idx1": ie,
                "idx2": imu,
                "vis": vis,
                "real": real,
                "fromTau1": ele["isTauDecay"],
                "fromTau2": mu["isTauDecay"],
                "rank": ele["p4"].P() + mu["p4"].P(),
            })
    return choose_best_pair(candidates)

def make_mumu_pair(muons, taus):
    candidates = []
    for imu1, mu1 in enumerate(muons):
        for imu2, mu2 in enumerate(muons):
           if imu2>imu1:
            if mu1["charge"] * mu2["charge"] >= 0:
                continue

            # Avoid the pathological case where both visible leptons are mapped
            # to the same tau copy.
            if mu1["parentTauIdx"] >= 0 and mu1["parentTauIdx"] == mu2["parentTauIdx"]:
                continue

            vis = mu1["p4"] + mu2["p4"]
            real = original_lepton_p4(mu1, taus) + original_lepton_p4(mu2, taus)
            candidates.append({
                "idx1": imu1,
                "idx2": imu2,
                "vis": vis,
                "real": real,
                "fromTau1": mu1["isTauDecay"],
                "fromTau2": mu2["isTauDecay"],
                "rank": mu1["p4"].P() + mu2["p4"].P(),
            })
    return choose_best_pair(candidates)

def make_ee_pair(electrons, taus):
    candidates = []
    for ie1, e1 in enumerate(electrons):
        for ie2, e2 in enumerate(electrons):
           if ie2>ie1:
            if e1["charge"] * e2["charge"] >= 0:
                continue

            # Avoid the pathological case where both visible leptons are mapped
            # to the same tau copy.
            if e1["parentTauIdx"] >= 0 and e1["parentTauIdx"] == e2["parentTauIdx"]:
                continue

            vis = e1["p4"] + e2["p4"]
            real = original_lepton_p4(e1, taus) + original_lepton_p4(e2, taus)
            candidates.append({
                "idx1": ie1,
                "idx2": ie2,
                "vis": vis,
                "real": real,
                "fromTau1": e1["isTauDecay"],
                "fromTau2": e2["isTauDecay"],
                "rank": e1["p4"].P() + e2["p4"].P(),
            })
    return choose_best_pair(candidates)

# ----------------------------------------------------------------------
# Output tree and branches
# ----------------------------------------------------------------------

treeName = "outtree"
outfile = ROOT.TFile(fileOutName, "RECREATE")
new_tree = ROOT.TTree(treeName, "GEN-level LFV variables")

scalar_f = {}
scalar_i = {}
vector_f = {}
vector_i = {}


def add_float(name, default=-999.0):
    scalar_f[name] = array("f", [default])
    new_tree.Branch(name, scalar_f[name], "%s/F" % name)


def add_int(name, default=-999):
    scalar_i[name] = array("i", [default])
    new_tree.Branch(name, scalar_i[name], "%s/I" % name)


def add_vfloat(name):
    vector_f[name] = ROOT.std.vector("float")()
    new_tree.Branch(name, vector_f[name])


def add_vint(name):
    vector_i[name] = ROOT.std.vector("int")()
    new_tree.Branch(name, vector_i[name])


def reset_branches():
    for b in scalar_f.values():
        b[0] = -999.0
    for b in scalar_i.values():
        b[0] = -999
    for b in vector_f.values():
        b.clear()
    for b in vector_i.values():
        b.clear()


# Event/object counts
for name in ["nMuons", "nElectrons", "nTaus", "nNeutrinos", "nPhotons"]:
    add_int(name, 0)

# Beam information
for name in ["beamE_0", "beamE_1", "beamE", "sqrtS"]:
    add_float(name)

# Final-state muons/electrons.
for obj in ["Muon", "Electron"]:
    for comp in ["px", "py", "pz", "e"]:
        add_vfloat("%s_%s" % (obj, comp))
    add_vint("%s_ID" % obj)
    add_vfloat("%s_charge" % obj)
    add_vint("%s_isTauDecay" % obj)
    add_vint("%s_parentTauIdx" % obj)

# Generator taus: both visible and true P4 are retained from tauReco.visTauGen.
for comp in ["px", "py", "pz", "e"]:
    add_vfloat("Tau_vis_%s" % comp)
    add_vfloat("Tau_true_%s" % comp)
add_vint("Tau_ID")
add_vint("Tau_TauID")
add_vfloat("Tau_charge")
add_vint("Tau_isHadronic")

# Status-1 neutrinos.
for comp in ["px", "py", "pz", "e"]:
    add_vfloat("Neutrino_%s" % comp)
add_vint("Neutrino_ID")
add_vfloat("Neutrino_charge")

for name in ["NeutrinoSum_px", "NeutrinoSum_py", "NeutrinoSum_pz", "NeutrinoSum_e", "NeutrinoSum_pt"]:
    add_float(name)

# Status-1 photons.
for comp in ["px", "py", "pz", "e"]:
    add_vfloat("Photon_%s" % comp)

add_vint("Photon_ID")

# Direct generator mother information.
add_vint("Photon_motherID")
add_vint("Photon_motherStatus")
add_vint("Photon_motherIdx")
add_vint("Photon_nMothers")

# Tau provenance.
add_vint("Photon_fromTauDecay")
add_vint("Photon_parentTauIdx")


# Z/H before it decays.
add_int("GenMother_ID")
for comp in ["px", "py", "pz", "e", "m"]:
    add_float("GenMother_%s" % comp)

# One best candidate per LFV channel. idx1/idx2 are:
#   MuTauh : Muon index, Tau index
#   ETauh  : Electron index, Tau index
#   EMu    : Electron index, Muon index
for channel in ["MuTauh", "ETauh", "EMu","TauhTauh","MuMu","EE"]:
    add_int(channel + "_valid", 0)
    add_int(channel + "_idx1")
    add_int(channel + "_idx2")
    add_int(channel + "_obj1FromTau", 0)
    add_int(channel + "_obj2FromTau", 0)
    for kind in ["vis", "real"]:
        for comp in ["px", "py", "pz", "e", "m"]:
            add_float("%s_%s_%s" % (channel, kind, comp))

# Which selection retained the event. 0=none, 1=MuTauh, 2=ETauh, 3=EMu, 4=TauhTauh, 5=MuMu, 6=EE
# 7=more than one channel valid (possible for inclusive/background events).
add_int("selectedChannel", 0)


# ----------------------------------------------------------------------
# Fill helpers
# ----------------------------------------------------------------------

def push_p4(prefix, p4):
    vector_f[prefix + "_px"].push_back(p4.Px())
    vector_f[prefix + "_py"].push_back(p4.Py())
    vector_f[prefix + "_pz"].push_back(p4.Pz())
    vector_f[prefix + "_e"].push_back(p4.E())


def fill_pair(channel, pair):
    if pair is None:
        scalar_i[channel + "_valid"][0] = 0
        return

    scalar_i[channel + "_valid"][0] = 1
    scalar_i[channel + "_idx1"][0] = pair["idx1"]
    scalar_i[channel + "_idx2"][0] = pair["idx2"]
    scalar_i[channel + "_obj1FromTau"][0] = int(pair["fromTau1"])
    scalar_i[channel + "_obj2FromTau"][0] = int(pair["fromTau2"])

    for kind in ["vis", "real"]:
        p4 = pair[kind]
        scalar_f["%s_%s_px" % (channel, kind)][0] = p4.Px()
        scalar_f["%s_%s_py" % (channel, kind)][0] = p4.Py()
        scalar_f["%s_%s_pz" % (channel, kind)][0] = p4.Pz()
        scalar_f["%s_%s_e" % (channel, kind)][0] = p4.E()
        scalar_f["%s_%s_m" % (channel, kind)][0] = p4.M()


# ----------------------------------------------------------------------
# Event loop
# ----------------------------------------------------------------------

totalEvents = 0
selectedEvents = 0

for event in reader.get("events"):
    if totalEvents % 10000 == 0:
        print(totalEvents)

    totalEvents += 1
    reset_branches()

    mc_particles = event.get(genparts)
    if len(mc_particles) < 2:
        continue

    # Beam energy. Keep both beam entries and the average beam energy.
    beamE_0 = float(mc_particles[0].getEnergy())
    beamE_1 = float(mc_particles[1].getEnergy())
    scalar_f["beamE_0"][0] = beamE_0
    scalar_f["beamE_1"][0] = beamE_1
    scalar_f["beamE"][0] = 0.5 * (beamE_0 + beamE_1)
    scalar_f["sqrtS"][0] = beamE_0 + beamE_1

    # Taus are built exclusively with the existing tauReco helper.
    taus = build_gen_taus(mc_particles)
    constituent_to_tau = build_tau_constituent_map(taus)

    # Store taus.
    for tau in taus:
        push_p4("Tau_vis", tau["vis_p4"])
        push_p4("Tau_true", tau["true_p4"])

        q = tau["charge"]
        pdgid = -15 if q > 0 else 15 if q < 0 else 0
        vector_i["Tau_ID"].push_back(pdgid)
        vector_i["Tau_TauID"].push_back(tau["TauID"])
        vector_f["Tau_charge"].push_back(q)
        vector_i["Tau_isHadronic"].push_back(int(tau_is_hadronic(tau["TauID"])))

    # Final-state particle collections.
    muons = []
    electrons = []
    neutrinos = []
    photons = []

    for part in mc_particles:
        if part.getGeneratorStatus() != 1:
            continue

        pdgid = int(part.getPDG())
        apid = abs(pdgid)

        if apid == 22:
            p4 = particle_p4(part)
        
            pho_info = get_photon_info(
                part,
                constituent_to_tau,
                taus
            )
        
            photons.append({
                "p4": p4,
                "ID": pdgid,
                "motherID": pho_info["motherID"],
                "motherStatus": pho_info["motherStatus"],
                "motherIdx": pho_info["motherIdx"],
                "nMothers": pho_info["nMothers"],
                "fromTauDecay": pho_info["fromTauDecay"],
                "parentTauIdx": pho_info["parentTauIdx"],
            })
        
            continue


        if apid not in (11, 12, 13, 14, 16):
            continue

        p4 = particle_p4(part)
        charge = float(part.getCharge())

        if apid in (11, 13):
            key = particle_key(part)
            parent_tau_idx = constituent_to_tau.get(key, -1) if key is not None else -1

            # Usually the status-1 lepton is directly among visTauGen constituents.
            # If generator FSR inserts another lepton copy, use ancestry only to
            # recover the link; the tau P4/TauID still come from tauReco.
            if parent_tau_idx < 0:
                tau_ancestor = find_tau_ancestor(part)
                parent_tau_idx = match_tau_index_from_p4(tau_ancestor, taus)

            record = {
                "p4": p4,
                "ID": pdgid,
                "charge": charge,
                "isTauDecay": int(parent_tau_idx >= 0),
                "parentTauIdx": parent_tau_idx,
            }

            if apid == 13:
                muons.append(record)
            else:
                electrons.append(record)

        else:
            neutrinos.append({"p4": p4, "ID": pdgid, "charge": charge})

    # Store status-1 muons and electrons, including tau provenance.
    for obj_name, objects in [("Muon", muons), ("Electron", electrons)]:
        for obj in objects:
            push_p4(obj_name, obj["p4"])
            vector_i[obj_name + "_ID"].push_back(obj["ID"])
            vector_f[obj_name + "_charge"].push_back(obj["charge"])
            vector_i[obj_name + "_isTauDecay"].push_back(obj["isTauDecay"])
            vector_i[obj_name + "_parentTauIdx"].push_back(obj["parentTauIdx"])

    # Store status-1 photons and their generator provenance.
    for pho in photons:
        push_p4("Photon", pho["p4"])
    
        vector_i["Photon_ID"].push_back(pho["ID"])
    
        vector_i["Photon_motherID"].push_back(pho["motherID"])
        vector_i["Photon_motherStatus"].push_back(pho["motherStatus"])
        vector_i["Photon_motherIdx"].push_back(pho["motherIdx"])
        vector_i["Photon_nMothers"].push_back(pho["nMothers"])
    
        vector_i["Photon_fromTauDecay"].push_back(pho["fromTauDecay"])
        vector_i["Photon_parentTauIdx"].push_back(pho["parentTauIdx"])



    # Store status-1 neutrinos and their vector sum.
    for nu in neutrinos:
        push_p4("Neutrino", nu["p4"])
        vector_i["Neutrino_ID"].push_back(nu["ID"])
        vector_f["Neutrino_charge"].push_back(nu["charge"])

    nu_px, nu_py, nu_pz, nu_e = sum_status1_neutrinos(mc_particles)
    scalar_f["NeutrinoSum_px"][0] = nu_px
    scalar_f["NeutrinoSum_py"][0] = nu_py
    scalar_f["NeutrinoSum_pz"][0] = nu_pz
    scalar_f["NeutrinoSum_e"][0] = nu_e
    scalar_f["NeutrinoSum_pt"][0] = math.sqrt(nu_px * nu_px + nu_py * nu_py)

    scalar_i["nMuons"][0] = len(muons)
    scalar_i["nElectrons"][0] = len(electrons)
    scalar_i["nTaus"][0] = len(taus)
    scalar_i["nNeutrinos"][0] = len(neutrinos)
    scalar_i["nPhotons"][0] = len(photons)

    # Exact generated Z/H four-vector before decay.
    boson = find_decay_boson(mc_particles)
    if boson is not None:
        bp4 = particle_p4(boson)
        scalar_i["GenMother_ID"][0] = int(boson.getPDG())
        scalar_f["GenMother_px"][0] = bp4.Px()
        scalar_f["GenMother_py"][0] = bp4.Py()
        scalar_f["GenMother_pz"][0] = bp4.Pz()
        scalar_f["GenMother_e"][0] = bp4.E()
        scalar_f["GenMother_m"][0] = bp4.M()

    # Build all three possible LFV pairs. "real" replaces a visible tau decay
    # product by the true parent tau P4; "vis" uses only visible objects.
    mutauh = make_mutauh_pair(muons, taus)
    etauh = make_etauh_pair(electrons, taus)
    emu = make_emu_pair(electrons, muons, taus)
    tauhtauh = make_tauhtauh_pair(taus)
    mumu = make_mumu_pair(muons,taus)
    ee = make_ee_pair(electrons, taus)

    fill_pair("MuTauh", mutauh)
    fill_pair("ETauh", etauh)
    fill_pair("EMu", emu)
    fill_pair("TauhTauh", tauhtauh)
    fill_pair("MuMu", mumu)
    fill_pair("EE", ee)


    valid = {
        "mutauh": mutauh is not None,
        "etauh": etauh is not None,
        "emu": emu is not None,
        "tauhtauh": tauhtauh is not None,
        "mumu": mumu is not None,
        "ee": ee is not None,
    }

    nvalid = sum(int(v) for v in valid.values())
    if nvalid == 1:
        if valid["mutauh"]:
            scalar_i["selectedChannel"][0] = 1
        elif valid["etauh"]:
            scalar_i["selectedChannel"][0] = 2
        elif valid["emu"]:
            scalar_i["selectedChannel"][0] = 3
        elif valid["tauhtauh"]:
            scalar_i["selectedChannel"][0] = 4
        elif valid["mumu"]:
            scalar_i["selectedChannel"][0] = 5
        elif valid["ee"]:
            scalar_i["selectedChannel"][0] = 6
    elif nvalid > 1:
        scalar_i["selectedChannel"][0] = 7
    else:
        scalar_i["selectedChannel"][0] = 0

    if args.channel == "all":
        keep_event = (nvalid > 0)
    else:
        keep_event = valid[args.channel]

#    if not keep_event:
#        continue

    new_tree.Fill()
    selectedEvents += 1


print("Total events:", totalEvents)
print("Selected events:", selectedEvents)

outfile.cd()
new_tree.Write()
outfile.Close()

print("Wrote:", fileOutName)

