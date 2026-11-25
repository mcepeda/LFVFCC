#!/usr/bin/env python3
"""Print decay trees and measure generator-level sample composition."""

import sys
import math
from collections import Counter
from podio.root_io import Reader

filename = sys.argv[1] if len(sys.argv) > 1 else "out_sim_edm4hep.root"
max_events = int(sys.argv[2]) if len(sys.argv) > 2 else 3

names = {
    23: "Z", 25: "H", 24: "W+", -24: "W-",
    13: "mu-", -13: "mu+", 15: "tau-", -15: "tau+",
    11: "e-", -11: "e+", 22: "gamma", 21: "g",
    12: "nu_e", -12: "anti_nu_e",
    14: "nu_mu", -14: "anti_nu_mu",
    16: "nu_tau", -16: "anti_nu_tau",
    1: "d", -1: "dbar", 2: "u", -2: "ubar",
    3: "s", -3: "sbar", 4: "c", -4: "cbar",
    5: "b", -5: "bbar", 6: "t", -6: "tbar",
    111: "pi0", 211: "pi+", -211: "pi-",
    113: "rho0", 213: "rho+", -213: "rho-",
    20213: "a1+", -20213: "a1-",
    221: "eta", 223: "omega", 331: "eta_prime",
    321: "K+", -321: "K-", 311: "K0", -311: "K0bar",
    310: "KS", 130: "KL", 313: "K*0", -313: "K*0bar",
    323: "K*+", -323: "K*-", 333: "phi",
}
self_conjugate = {21, 22, 23, 25, 111, 113, 221, 223, 331, 333, 130, 310}
neutrinos = {12, 14, 16}
tau_decay_products = {
    11, 12, 13, 14, 16, 22, 111, 211,
    311, 321, 310, 130, 221, 2112, 2212,
}


def pdg_label(pdg):
    return names.get(pdg, str(pdg))


def particle_id(p):
    oid = p.getObjectID()
    return (getattr(oid, "collectionID", -1), oid.index)


def daughters(p):
    return [d for d in p.getDaughters() if not d.isCreatedInSimulation()]


def tree(p, depth=0, seen=None):
    if seen is None:
        seen = set()
    idx = particle_id(p)
    if idx in seen:
        return
    seen.add(idx)
    print(f"{'  '*depth}{pdg_label(p.getPDG()):>10} "
          f"[{p.getObjectID().index:3d}]  PDG={p.getPDG():>6}  "
          f"status={p.getGeneratorStatus():2d}  E={p.getEnergy():7.2f} GeV")
    if depth >= 12:
        return
    for d in daughters(p):
        tree(d, depth+1, seen)


def terminal_copies(mc, pdg):
    return [p for p in mc if abs(p.getPDG()) == abs(pdg)
            and not p.isCreatedInSimulation()
            and not any(d.getPDG() == p.getPDG() for d in daughters(p))]


def canonical_channel(pdgs, conjugate=False):
    if conjugate:
        pdgs = [x if abs(x) in self_conjugate else -x for x in pdgs]
    return tuple(sorted(pdgs, key=lambda x: (abs(x), -x)))


def direct_decay(p):
    ds = daughters(p)
    if not ds:
        return None
    return canonical_channel([d.getPDG() for d in ds], p.getPDG() == -15)


def channel_label(channel):
    return " ".join(f"{n} {pdg_label(pdg)}" if n > 1 else pdg_label(pdg)
                    for pdg, n in Counter(channel).items())


def tau_final_products(p):
    products, seen = [], set()

    def visit(q):
        idx = particle_id(q)
        if idx in seen:
            return
        seen.add(idx)
        ds = daughters(q)
        if abs(q.getPDG()) in tau_decay_products or not ds:
            products.append(q)
            return
        for d in ds:
            visit(d)

    for d in daughters(p):
        visit(d)
    return products


def tau_decay_mode(products):
    visible = [p for p in products if abs(p.getPDG()) not in neutrinos
               and abs(p.getPDG()) != 22]
    pdgs = [abs(p.getPDG()) for p in visible]
    if pdgs == [11]:
        return "Electron"
    if pdgs == [13]:
        return "Muon"
    if not visible:
        return "Unclassified"
    return f"Hadronic {sum(abs(p.getCharge()) > 0.1 for p in visible)}-prong"


def tau_decay_mode_pi0(products):
    mode = tau_decay_mode(products)
    if not mode.startswith("Hadronic"):
        return mode
    return f"{mode}, {sum(p.getPDG() == 111 for p in products)} pi0"


def print_fractions(title, counts, formatter=str):
    total = sum(counts.values())
    print("\n" + "="*85)
    print(title)
    print(f"Total decays: {total:,}")
    print("="*85)
    if total == 0:
        print("No decays found")
        return
    print(f'{"Decay channel":<45}{"N":>12}{"BR [%]":>12}{"Stat. [%]":>12}')
    print("-"*85)
    for channel, n in counts.most_common():
        br = n / total
        err = math.sqrt(br * (1.0-br) / total)
        print(f"{formatter(channel):<45}{n:>12,}{100*br:>12.5f}{100*err:>12.5f}")
    print("-"*85)
    print(f'{"Total":<45}{total:>12,}{100.:>12.5f}')


def main():
    n_events = 0
    n_particles, n_undecayed, n_events_with_particle = Counter(), Counter(), Counter()
    decay_counts = {pdg: Counter() for pdg in (23, 25, 15)}
    tau_exclusive_counts, tau_mode_counts = Counter(), Counter()
    tau_mode_pi0_counts, tau_charge_counts = Counter(), Counter()
    reader = Reader(filename)

    for ievt, event in enumerate(reader.get("events")):
        n_events += 1
        mc = event.get("MCParticles")
        if ievt < max_events:
            counts = Counter(p.getPDG() for p in mc)
            print(f"\n========== Event {ievt} | {len(mc)} MCParticles ==========")
            print("PDG counts (all MCParticles):", counts.most_common(20))
            for pdg in (23, 25, 15):
                roots = [p for p in mc if abs(p.getPDG()) == pdg
                         and not p.isCreatedInSimulation()
                         and not any(parent.getPDG() == p.getPDG()
                                     for parent in p.getParents())]
                print(f"\n{pdg_label(pdg)} decay tree:")
                if not roots:
                    print("  Not found in MCParticles")
                for p in roots:
                    tree(p)

        for pdg in (23, 25, 15):
            particles = terminal_copies(mc, pdg)
            if particles:
                n_events_with_particle[pdg] += 1
            for p in particles:
                n_particles[pdg] += 1
                if pdg == 15:
                    tau_charge_counts[p.getPDG()] += 1
                channel = direct_decay(p)
                if channel is None:
                    n_undecayed[pdg] += 1
                    continue
                decay_counts[pdg][channel] += 1
                if pdg == 15:
                    products = tau_final_products(p)
                    tau_exclusive_counts[canonical_channel(
                        [d.getPDG() for d in products], p.getPDG() == -15)] += 1
                    tau_mode_counts[tau_decay_mode(products)] += 1
                    tau_mode_pi0_counts[tau_decay_mode_pi0(products)] += 1
        if (ievt+1) % 100000 == 0:
            print(f"Processed {ievt+1:,} events")

    print("\n" + "="*85)
    print("FULL SAMPLE SUMMARY")
    print("="*85)
    print(f"Total events processed: {n_events:,}")
    for pdg in (23, 25, 15):
        print(f"\nParticle: {pdg_label(pdg) if pdg != 15 else 'tau (both charges)'}")
        print(f"  Events containing particle: {n_events_with_particle[pdg]:,}")
        print(f"  Total terminal copies: {n_particles[pdg]:,}")
        print(f"  Decaying particles: {sum(decay_counts[pdg].values()):,}")
        print(f"  Without generator-level daughters: {n_undecayed[pdg]:,}")
    print(f"\nTotal tau-: {tau_charge_counts[15]:,}")
    print(f"Total tau+: {tau_charge_counts[-15]:,}")

    for pdg, title in [(23, "Z"), (25, "HIGGS"), (15, "TAU")]:
        print_fractions(title + " DIRECT DECAY BRANCHING FRACTIONS",
                        decay_counts[pdg], channel_label)
    print_fractions("TAU EXCLUSIVE DECAYS (RESONANCES RESOLVED)",
                    tau_exclusive_counts, channel_label)
    print_fractions("TAU INCLUSIVE DECAY MODES", tau_mode_counts)
    print_fractions("TAU DECAY MODES BY PRONG AND PI0 MULTIPLICITY",
                    tau_mode_pi0_counts)


if __name__ == "__main__":
    main()
