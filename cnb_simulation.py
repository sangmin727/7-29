import argparse
import time

import numpy as np


GRID_K = 100
GRID_COD = 150


def J_cr(s_inf, l_cr, R_0, E):
    """Energy release rate at the crack tip."""
    norm_a = l_cr / R_0
    K_cr_value = (
        s_inf
        * np.sqrt(np.pi * l_cr)
        * (1.0 / ((1 - norm_a) ** 1.5))
        * (1.122 - 1.302 * norm_a + 0.988 * norm_a**2 - 0.308 * norm_a**3)
    )
    return (K_cr_value**2) / E


def G_SIF_CNB(a, R_0, x):
    x11 = a - x
    R_l = R_0 - a
    return G_SIF_INF(R_l, x11) * f(R_l, R_0, x11)


def G_SIF_INF(R_l, x):
    numerator = 2 * np.pi * (R_l + x)
    denominator = (np.pi * R_l) ** 1.5

    ratio = R_l / (R_l + x)
    ratio = np.clip(ratio, -1.0, 1.0)
    acos_term = np.arccos(ratio)

    radicand = ((R_l + x) ** 2) - (R_l**2)
    radicand = np.maximum(radicand, np.finfo(float).eps)
    sqrt_term = R_l / np.sqrt(radicand)

    return (numerator / denominator) * (acos_term + sqrt_term)


def K_cr(s_inf, l_cr, R_0):
    """SIF at the crack tip."""
    norm_a = l_cr / R_0
    return (
        s_inf
        * np.sqrt(np.pi * l_cr)
        * (1.0 / ((1 - norm_a) ** 1.5))
        * (1.122 - 1.302 * norm_a + 0.988 * norm_a**2 - 0.308 * norm_a**3)
    )


def K_dr(s_dr, l_cr, l_pz, R_0):
    """SIF due to the drawing stress s_dr."""
    if l_cr == l_pz:
        return 0

    n = GRID_K
    h = (l_pz - l_cr) / n
    x_h = np.arange(l_cr + h / 2, l_pz, h)
    Kx = G_SIF_CNB(l_pz, R_0, x_h) * h
    return s_dr * np.sum(Kx)


def K_inf(s_inf, l_pz, R_0):
    """SIF due to the remote stress s_inf."""
    n = GRID_K
    h = l_pz / n
    x_h = np.arange(h / 2, l_pz, h)
    Kx = G_SIF_CNB(l_pz, R_0, x_h) * h
    return s_inf * np.sum(Kx)


def K_tot(s_inf, s_dr, l_cr, l_pz, R_0):
    K_tot_value = K_inf(s_inf, l_pz, R_0) - K_dr(s_dr, l_cr, l_pz, R_0)
    return max(K_tot_value, 0)


def V_AZ_simp(s_inf, s_dr, l_cr, l_pz, R_0, E, _lambda):
    a = l_cr
    b = l_pz
    COD_a = COD_tot(s_inf, s_dr, a, l_cr, l_pz, R_0, E)
    COD_mid1 = COD_tot(s_inf, s_dr, ((2 * a + b) / 3), l_cr, l_pz, R_0, E)
    COD_mid2 = COD_tot(s_inf, s_dr, ((a + 2 * b) / 3), l_cr, l_pz, R_0, E)
    COD_b = COD_tot(s_inf, s_dr, b, l_cr, l_pz, R_0, E)
    integral = ((b - a) / 8) * (
        ((R_0 - a) * COD_a)
        + 3 * ((R_0 - ((2 * a + b) / 3)) * COD_mid1)
        + 3 * ((R_0 - ((a + 2 * b) / 3)) * COD_mid2)
        + ((R_0 - b) * COD_b)
    )
    return (2 * np.pi / (_lambda - 1)) * integral


def COD_tot(s_inf, s_dr, x1, l_cr, l_pz, R_0, E):
    """Compute COD_tot at a position x1."""
    n = GRID_COD
    dxi = (l_pz - x1) / n
    ds1 = x1 / n
    ds2 = (l_pz - x1) / n

    if l_cr == l_pz or x1 >= l_pz:
        cod_total = 0
    else:
        I = np.arange(1, n + 1)
        J = np.arange(1, n + 1)
        i, j = np.meshgrid(I, J)

        xi = x1 + dxi * (i - 0.5)
        s1 = ds1 * (j - 0.5)

        COD1_segment = (
            ((2 * s_inf) / E)
            * G_SIF_CNB(xi, R_0, x1)
            * G_SIF_CNB(xi, R_0, s1)
            * dxi
            * ds1
        )
        COD1 = np.sum(COD1_segment)

        K = np.arange(1, n)
        M = np.arange(1, n)
        m, k = np.meshgrid(K, M)

        xi_k = x1 + dxi * (m + 0.5)
        s2 = x1 + ds2 * (k - 0.5)

        COD2_segment = np.triu(
            ((2 * s_inf) / E)
            * G_SIF_CNB(xi_k, R_0, x1)
            * G_SIF_CNB(xi_k, R_0, s2)
            * dxi
            * ds2
        )
        COD2 = np.sum(COD2_segment)

        if x1 == l_cr:
            COD_dr_segment = np.triu(
                ((2 * s_dr) / E)
                * G_SIF_CNB(xi_k, R_0, x1)
                * G_SIF_CNB(xi_k, R_0, s2)
                * dxi
                * ds2
            )
            COD_dr = np.sum(COD_dr_segment)
        else:
            ds_dr1 = (x1 - l_cr) / n
            ds_dr2 = (l_pz - x1) / n
            dxi_dr = (l_pz - x1) / n

            H = np.arange(1, n + 1)
            G_vals = np.arange(1, n + 1)
            g, h = np.meshgrid(H, G_vals)

            xi_dr = x1 + dxi_dr * (h - 0.5)
            s_dr1 = ds_dr1 * (g - 0.5)

            COD_dr_1_segment = (
                ((2 * s_dr) / E)
                * G_SIF_CNB(xi_dr, R_0, x1)
                * G_SIF_CNB(xi_dr, R_0, s_dr1)
                * dxi_dr
                * ds_dr1
            )
            COD_dr_1 = np.sum(COD_dr_1_segment)

            Q = np.arange(1, n)
            W = np.arange(1, n)
            w, q = np.meshgrid(Q, W)

            xi_dr2 = x1 + dxi_dr * (w + 0.5)
            s_dr2 = x1 + ds_dr2 * (q - 0.5)

            COD_dr_2_segment = np.triu(
                ((2 * s_dr) / E)
                * G_SIF_CNB(xi_dr2, R_0, x1)
                * G_SIF_CNB(xi_dr2, R_0, s_dr2)
                * dxi_dr
                * ds_dr2
            )
            COD_dr_2 = np.sum(COD_dr_2_segment)

            COD_dr = COD_dr_1 + COD_dr_2

        cod_total = COD1 + COD2 - COD_dr

    return max(cod_total, 0)


def f(R_l, R_0, x):
    alpha = x / (R_0 - R_l)
    gamma_val = R_l / R_0
    f3 = -0.2007 * gamma_val + 0.06840
    f2 = -0.1001 + 0.2671 * gamma_val + 0.1260 * (gamma_val**2)
    f1 = 0.2651 * gamma_val - 1.210 * (gamma_val**2) + 1.590 * (gamma_val**3)
    f0 = 0.002101 * gamma_val + 1
    return (alpha**3) * f3 + (alpha**2) * f2 + alpha * f1 + f0


def R_AZ(s_inf, s_dr, l_cr, l_pz, R_0, E, _lambda):
    hh = l_pz / 10000
    V_plus = V_AZ_simp(s_inf, s_dr, l_cr, l_pz + hh, R_0, E, _lambda)
    V_minus = V_AZ_simp(s_inf, s_dr, l_cr, l_pz - hh, R_0, E, _lambda)
    return (1 / (2 * np.pi * (R_0 - l_pz))) * (V_plus - V_minus) / (2 * hh)


def run_simulation(use_dt_current_for_growth=False, n_t=1000):
    """
    Run the CNB simulation.

    If use_dt_current_for_growth=False, this matches the original provided script,
    where growth updates use `dt`.
    """
    dt_initial = 10
    dt_middle = 10
    dt = 10
    s_inf = 10
    s_dr = 16
    l = 1.5
    R_0 = 5.0
    E = 900
    gamma_0 = 15.0
    gamma_tr = 35.0
    _lambda = 5.0
    t_star = 260
    r = 1.0
    k_cr = 0.002
    k_pz = 0.002
    k1c = 2.75

    l_cr = [l]
    l_pz = [l]
    gamma = [gamma_0]
    t = [0]
    step = []
    SIF = []

    for i in range(n_t):
        if t[i] < 200:
            dt_current = dt_initial
        elif t[i] < 1000:
            dt_current = dt_middle
        else:
            dt_current = dt

        t_new = t[i] + dt_current

        K_tot_value = K_tot(s_inf, s_dr, l_cr[i], l_pz[i], R_0)
        term1 = (K_tot_value**2) / E
        R_AZ_value = R_AZ(s_inf, s_dr, l_cr[i], l_pz[i], R_0, E, _lambda)
        term2 = gamma_tr * R_AZ_value

        growth_dt = dt_current if use_dt_current_for_growth else dt

        if (term1 - term2) > 0:
            l_pz_new = l_pz[i] + growth_dt * k_pz * (term1 - term2)
        else:
            l_pz_new = l_pz[i]

        J_cr_value = J_cr(s_inf, l_cr[i], R_0, E)
        if (J_cr_value - 2 * gamma[i]) > 0:
            l_cr_new = l_cr[i] + growth_dt * k_cr * (J_cr_value - 2 * gamma[i])
        else:
            l_cr_new = l_cr[i]

        if l_cr_new > l_pz_new:
            l_cr_new = l_pz_new

        gamma_new = gamma[i]
        if l_cr_new == l_pz_new:
            gamma_new = gamma_0
        else:
            for j in range(i + 1):
                idx1 = i - j
                idx2 = i - j + 1
                if idx1 >= 0 and idx2 < len(l_pz):
                    if l_pz[idx1] <= l_cr_new < l_pz[idx2]:
                        gamma_new = gamma_0 / (1 + (((t_new - t[idx1]) / t_star) ** r))
                        break

        t.append(t_new)
        l_pz.append(l_pz_new)
        l_cr.append(l_cr_new)
        gamma.append(gamma_new)

        if (
            i >= 11
            and l_cr[i - 9] == l_cr[i]
            and l_cr[i - 3] == l_cr[i]
            and l_cr[i - 1] == l_cr[i]
            and l_cr_new > l_cr[i]
            and (l_cr_new - l_cr[i]) > (l_cr[i] - l_cr[i - 1])
        ):
            step.append(i)

        K_cr_value = K_cr(s_inf, l_cr[i], R_0)
        if K_cr_value >= (k1c * (10**1.5)):
            l_pz[-1] = R_0
            l_cr[-1] = R_0
            break
        if l_pz_new >= R_0:
            l_pz[-1] = R_0
            break
        if l_cr_new >= R_0 or l_pz_new >= R_0:
            break

    t = np.array(t)
    l_cr = np.array(l_cr)
    l_pz = np.array(l_pz)

    m = len(step)
    if m > 1:
        lengthstep = l_pz[step] - l_cr[step]
        stepcrack = l_cr[step]
        timestep = np.zeros(m)
        timestep[1:] = np.diff(t[step])
        timestep[0] = t[step][0]

        dadt = lengthstep / timestep
        logdadt = np.log10(dadt)

        for idx in range(m):
            SIF.append(K_cr(s_inf, stepcrack[idx], R_0))

        SIF = np.array(SIF)
        logSIF = np.log10(SIF)
        fit = np.polyfit(logSIF, logdadt, 1)
        slopefit = fit[0] * logSIF + fit[1]
    else:
        logdadt = np.array([])
        logSIF = np.array([])
        slopefit = np.array([])
        fit = None

    return {
        "t": t,
        "l_cr": l_cr,
        "l_pz": l_pz,
        "m": m,
        "fit": fit,
        "logdadt": logdadt,
        "logSIF": logSIF,
        "slopefit": slopefit,
    }


def plot_results(results):
    import matplotlib.pyplot as plt

    t = results["t"]
    l_cr = results["l_cr"]
    l_pz = results["l_pz"]
    m = results["m"]

    plt.figure(figsize=(14, 6))

    plt.subplot(1, 2, 1)
    plt.plot(t, l_cr, "-k", label="Crack length (mm)")
    plt.plot(t, l_pz, "-r", label="Crack layer length (mm)")
    plt.axis([0, t[-1] + 1000, 0, 7])
    plt.legend()
    plt.xlabel("Time (sec)", fontsize=15)
    plt.ylabel("Length (mm)", fontsize=15)
    plt.title(
        "CNB specimen, Surface energy "
        "$\\gamma= \\frac{{\\gamma}_{o}}{1 + {\\left(\\frac{{t_i - t_o}}{{t^*}}\\right)}^{r}}$",
        fontsize=14,
    )
    plt.grid(True)

    if m > 1:
        plt.subplot(1, 2, 2)
        plt.scatter(results["logSIF"], results["logdadt"])
        plt.plot(results["logSIF"], results["slopefit"], "r-.")
        plt.xlabel("log(SIF), MPa·m$^{0.5}$", fontsize=15)
        plt.ylabel("log(da/dt), mm/s", fontsize=15)
        fit = results["fit"]
        plt.text(
            results["logSIF"][0],
            results["logdadt"][-1] + 0.1,
            f"$\\log(da/dt) = {fit[0]:.3f}\\; \\log(SIF) + {fit[1]:.3f}$",
            fontsize=14,
        )
        plt.grid(True)
    else:
        plt.subplot(1, 2, 2)
        plt.text(
            0.5,
            0.5,
            "Insufficient data for plot",
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=15,
        )
        plt.axis("off")

    plt.tight_layout()
    plt.show()


def compare_original_and_revised(n_t=1000):
    original = run_simulation(use_dt_current_for_growth=False, n_t=n_t)
    revised = run_simulation(use_dt_current_for_growth=True, n_t=n_t)

    print("Comparison of original vs revised growth-time-step handling")
    print(
        f"- original: final t={original['t'][-1]:.2f}, "
        f"l_cr={original['l_cr'][-1]:.6f}, l_pz={original['l_pz'][-1]:.6f}, "
        f"steps={len(original['t'])}"
    )
    print(
        f"- revised : final t={revised['t'][-1]:.2f}, "
        f"l_cr={revised['l_cr'][-1]:.6f}, l_pz={revised['l_pz'][-1]:.6f}, "
        f"steps={len(revised['t'])}"
    )
    print(
        f"- deltas  : Δl_cr={revised['l_cr'][-1] - original['l_cr'][-1]:.6e}, "
        f"Δl_pz={revised['l_pz'][-1] - original['l_pz'][-1]:.6e}, "
        f"Δt={revised['t'][-1] - original['t'][-1]:.6e}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=["original", "revised", "compare"],
        default="compare",
        help="Run the original version, revised version, or compare both.",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Plot results for the selected single mode (requires matplotlib).",
    )
    parser.add_argument("--n-t", type=int, default=1000, help="Number of time iterations.")
    parser.add_argument(
        "--fast",
        action="store_true",
        help="Use reduced integration grids for faster approximate comparisons.",
    )
    args = parser.parse_args()

    start_time = time.time()

    global GRID_K, GRID_COD
    if args.fast:
        GRID_K = 40
        GRID_COD = 60

    if args.mode == "compare":
        compare_original_and_revised(n_t=args.n_t)
    else:
        use_dt_current = args.mode == "revised"
        results = run_simulation(use_dt_current_for_growth=use_dt_current, n_t=args.n_t)
        print(
            f"mode={args.mode}, final t={results['t'][-1]:.2f}, "
            f"l_cr={results['l_cr'][-1]:.6f}, l_pz={results['l_pz'][-1]:.6f}, "
            f"steps={len(results['t'])}"
        )
        if args.plot:
            plot_results(results)

    end_time = time.time()
    print(f"Elapsed time: {end_time - start_time:.2f} seconds")


if __name__ == "__main__":
    main()
