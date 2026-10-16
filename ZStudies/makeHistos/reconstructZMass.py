import math
import ROOT

MTAU = 1.77686


def reconstruct_mass_ltauh(promptP4, tauVisP4,
                           nu_px, nu_py,
                           ecm=91.2):
    """
    prompt lepton + tau_h.

    The tau invisible system is one massless neutrino.
    Tau mass gives pz(nu).
    sqrt(s) is used only to choose between the two pz solutions.
    """

    pt2_nu = nu_px*nu_px + nu_py*nu_py

    A = (
        0.5 * (MTAU*MTAU - tauVisP4.M2())
        + tauVisP4.Px()*nu_px
        + tauVisP4.Py()*nu_py
    )

    a = (
        tauVisP4.E()*tauVisP4.E()
        - tauVisP4.Pz()*tauVisP4.Pz()
    )

    if a <= 0:
        return -999., -999., 0

    disc = A*A - a*pt2_nu

    real_solution = 1
    if disc < 0:
        disc = 0.
        real_solution = 0

    sqrt_disc = math.sqrt(disc)

    pz_solutions = [
        (A*tauVisP4.Pz() + tauVisP4.E()*sqrt_disc) / a,
        (A*tauVisP4.Pz() - tauVisP4.E()*sqrt_disc) / a,
    ]

    best_mass = -999.
    best_pz = -999.
    best_score = 1.e30

    for nu_pz in pz_solutions:

        nu_e = math.sqrt(pt2_nu + nu_pz*nu_pz)

        nuP4 = ROOT.TLorentzVector()
        nuP4.SetPxPyPzE(
            nu_px,
            nu_py,
            nu_pz,
            nu_e
        )

        mother = promptP4 + tauVisP4 + nuP4

        # Beam information only resolves the quadratic ambiguity.
        score = (
            (mother.E() - ecm)**2
            + mother.Px()**2
            + mother.Py()**2
            + mother.Pz()**2
        )

        if score < best_score:
            best_score = score
            best_mass = mother.M()
            best_pz = nu_pz

    return best_mass, best_pz, real_solution


def reconstruct_mass_ll_constrained(promptP4, tauLepP4,
                                    nu_px, nu_py,
                                    ecm=91.2):
    """
    prompt lepton + tau -> prompt lepton + lepton + nu + nu.

    Treat nu+nu as ONE effective massless invisible object.

        q = p_nu1 + p_nu2
        q^2 ~= 0

    Then impose

        (tauLep + q)^2 = m_tau^2

    to solve qz.

    sqrt(s) is used only to choose between the two qz solutions.

    Returns:
        constrained_mass, qz, real_solution
    """

    qt2 = nu_px*nu_px + nu_py*nu_py

    A = (
        0.5 * (MTAU*MTAU - tauLepP4.M2())
        + tauLepP4.Px()*nu_px
        + tauLepP4.Py()*nu_py
    )

    a = (
        tauLepP4.E()*tauLepP4.E()
        - tauLepP4.Pz()*tauLepP4.Pz()
    )

    if a <= 0:
        return -999., -999., 0

    disc = A*A - a*qt2

    real_solution = 1
    if disc < 0:
        disc = 0.
        real_solution = 0

    sqrt_disc = math.sqrt(disc)

    qz_solutions = [
        (A*tauLepP4.Pz() + tauLepP4.E()*sqrt_disc) / a,
        (A*tauLepP4.Pz() - tauLepP4.E()*sqrt_disc) / a,
    ]

    best_mass = -999.
    best_qz = -999.
    best_score = 1.e30

    for qz in qz_solutions:

        # Effective nu-nu system approximated as massless.
        qE = math.sqrt(qt2 + qz*qz)

        qP4 = ROOT.TLorentzVector()
        qP4.SetPxPyPzE(
            nu_px,
            nu_py,
            qz,
            qE
        )

        mother = promptP4 + tauLepP4 + qP4

        # Again, ecm only chooses which quadratic solution is preferable.
        score = (
            (mother.E() - ecm)**2
            + mother.Px()**2
            + mother.Py()**2
            + mother.Pz()**2
        )

        if score < best_score:
            best_score = score
            best_mass = mother.M()
            best_qz = qz

    return best_mass, best_qz, real_solution


def reconstruct_mass_collinear(promptP4, tauVisP4,
                               nu_px, nu_py):
    """
    Generic collinear tau estimator.

    Works for both:
      prompt lepton + tau_h
      prompt lepton + tau -> lepton nu nu

    Assumes the missing transverse momentum associated with the tau
    is parallel to the visible tau decay product/system.

    Returns:
        mass, alpha, nuPerp
    """

    pt2 = (
        tauVisP4.Px()*tauVisP4.Px()
        + tauVisP4.Py()*tauVisP4.Py()
    )

    if pt2 <= 0:
        return -999., -999., -999.

    alpha = (
        nu_px*tauVisP4.Px()
        + nu_py*tauVisP4.Py()
    ) / pt2

    nu_parallel_x = alpha * tauVisP4.Px()
    nu_parallel_y = alpha * tauVisP4.Py()

    nu_perp = math.sqrt(
        (nu_px - nu_parallel_x)**2
        + (nu_py - nu_parallel_y)**2
    )

    if alpha < 0:
        return -999., alpha, nu_perp

    scale = 1.0 + alpha

    tau_px = scale * tauVisP4.Px()
    tau_py = scale * tauVisP4.Py()
    tau_pz = scale * tauVisP4.Pz()

    tau_e = math.sqrt(
        tau_px*tau_px
        + tau_py*tau_py
        + tau_pz*tau_pz
        + MTAU*MTAU
    )

    tauP4 = ROOT.TLorentzVector()
    tauP4.SetPxPyPzE(
        tau_px,
        tau_py,
        tau_pz,
        tau_e
    )

    mother = promptP4 + tauP4

    return mother.M(), alpha, nu_perp





