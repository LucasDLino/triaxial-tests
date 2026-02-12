"""
VALIDAÇÃO E VERIFICAÇÃO DO MODELO MOHR-COULOMB
===============================================

Arquivo unificado contendo:

  PARTE 1 — Verificação numérica (PASS/FAIL)
    Testes automatizados com critério analítico/qualitativo.
    Garante que o modelo está correto após qualquer alteração de código.

  PARTE 2 — Validação gráfica
    Gráficos de q vs εa, εv vs εa, poropressão, círculos de Mohr.
    Permite inspeção visual do comportamento do modelo.

Parâmetros de referência:
  E = 30000 kPa, ν = 0.3, φ = 30°, c = 20 kPa
  σ3 = 100 kPa (referência para PASS/FAIL)
  σ3_list = [50, 100, 200] kPa (para gráficos)

Uso:
  python validacao.py              → roda tudo (verificação + gráficos salvos)
  python validacao.py --check      → apenas verificação PASS/FAIL
  python validacao.py --plot       → apenas gráficos (salva PNG)
  python validacao.py --plot --show → gráficos salvos + exibe na tela
"""

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mohr_courlomb import MohrCoulombModel
from triaxial import TriaxialTest


# =====================================================================
# PARÂMETROS GLOBAIS
# =====================================================================
E = 30000.0        # kPa
NU = 0.3
PHI_DEG = 30.0     # graus
C0 = 20.0          # kPa  (coesão)
SIGMA3_REF = 100.0 # kPa  (confinamento de referência)
SIGMA3_LIST = [50, 100, 200]  # kPa  (para gráficos multi-σ3)

SIN_PHI = np.sin(np.radians(PHI_DEG))
COS_PHI = np.cos(np.radians(PHI_DEG))
KP = (1 + SIN_PHI) / (1 - SIN_PHI)
Q_ANALITICO = SIGMA3_REF * (KP - 1) + 2 * C0 * np.sqrt(KP)


# =====================================================================
# HELPERS
# =====================================================================
_PASS = 0
_FAIL = 0
_TOTAL = 0


def _check(cond, label, detail=""):
    global _PASS, _FAIL, _TOTAL
    _TOTAL += 1
    if cond:
        _PASS += 1
        print(f"  [PASS] {label}")
    else:
        _FAIL += 1
        print(f"  [FAIL] {label}  {detail}")


def _make_model(psi_deg=0, cohesion=C0, sampling_pairs=None):
    return MohrCoulombModel(E, NU, PHI_DEG, cohesion,
                            psi_deg=psi_deg,
                            sampling_pairs=sampling_pairs,
                            compression_positive=True)


def _run(test_type, psi_deg=0, cohesion=C0, sampling_pairs=None,
         sigma3=SIGMA3_REF, steps=200, eps_max=0.15):
    model = _make_model(psi_deg, cohesion, sampling_pairs)
    test = TriaxialTest(model, sigma3, test_type=test_type)
    return test.run(eps_max=eps_max, steps=steps), model


# =====================================================================
#  PARTE 1 — VERIFICAÇÃO NUMÉRICA  (PASS / FAIL)
# =====================================================================

def _test_cd_psi0():
    """CD ψ=0: q ≈ q_analítico, εv > 0 (contração elástica)."""
    print("\n" + "=" * 70)
    print("TESTE 1: CD ψ=0 — Validação analítica")
    print("=" * 70)
    print(f"  q_analítico = σ3·(Kp−1) + 2c·√Kp = {Q_ANALITICO:.2f} kPa")

    res, _ = _run("CD", psi_deg=0)
    q_max = max(res['q'])
    err = abs(q_max - Q_ANALITICO) / Q_ANALITICO * 100

    print(f"  q_simulado  = {q_max:.2f} kPa  (erro {err:.2f}%)")
    _check(err < 1.0, "q erro < 1%", f"err={err:.2f}%")

    ev = res['volumetric_strain'][-1]
    print(f"  εv_final    = {ev:.6f}")
    _check(ev > 0, "εv > 0 (contração elástica em drenado)", f"εv={ev:.6f}")


def _test_cu_psi0():
    """CU ψ=0: Terzaghi consistente, εv=0, q_CU < q_CD."""
    print("\n" + "=" * 70)
    print("TESTE 2: CU ψ=0 — Consistência de Terzaghi")
    print("=" * 70)

    res, _ = _run("CU", psi_deg=0)

    # u ≈ σ3_total − σ3'
    max_err = max(abs(res['pore_pressure'][i] - (SIGMA3_REF - res['sigma3'][i]))
                  for i in range(len(res['pore_pressure'])))
    print(f"  max|u − (σ3_total−σ3')| = {max_err:.4f} kPa")
    _check(max_err < 1.0, "Terzaghi u ≈ σ3_total − σ3'")

    # εv = 0
    max_ev = max(abs(v) for v in res['volumetric_strain'])
    _check(max_ev < 1e-6, "εv = 0 (não drenado)", f"|εv|={max_ev:.2e}")

    # q_CU < q_CD
    q_cu = max(res['q'])
    print(f"  q_CU = {q_cu:.2f},  q_CD_analítico = {Q_ANALITICO:.2f}")
    _check(q_cu > 0, "q_CU > 0")
    _check(q_cu < Q_ANALITICO, "q_CU < q_CD (correto para φ>0, NC)")


def _test_uu():
    """UU: cu independe de σ3."""
    print("\n" + "=" * 70)
    print("TESTE 3: UU — cu independente de σ3")
    print("=" * 70)

    sigmas = [50, 100, 200, 400]
    q_vals = []
    for s3 in sigmas:
        r, _ = _run("UU", sigma3=s3)
        q = max(r['q'])
        q_vals.append(q)
        print(f"  σ3={s3:4d}: q_max = {q:.2f} kPa")

    var = (max(q_vals) - min(q_vals)) / q_vals[0] * 100
    _check(var < 1.0, "variação < 1%", f"var={var:.2f}%")


def _test_cd_cu_relation():
    """CD vs CU: q_CU < q_CD para φ>0 (NC)."""
    print("\n" + "=" * 70)
    print("TESTE 4: CD vs CU — relação qualitativa")
    print("=" * 70)

    r_cd, _ = _run("CD", psi_deg=0)
    r_cu, _ = _run("CU", psi_deg=0)
    q_cd = max(r_cd['q'])
    q_cu = max(r_cu['q'])
    ratio = q_cu / q_cd * 100
    print(f"  q_CD={q_cd:.2f}, q_CU={q_cu:.2f}, ratio={ratio:.1f}%")
    _check(q_cu < q_cd, "q_CU < q_CD", f"ratio={ratio:.1f}%")
    _check(q_cu > 0.3 * q_cd, "q_CU > 30% q_CD (não colapsou)")


def _test_cd_psi_neq0():
    """CD ψ≠0: q ≈ constante, εv dilata com ψ crescente."""
    print("\n" + "=" * 70)
    print("TESTE 5: CD ψ≠0 — q constante, εv dilata com ψ")
    print("=" * 70)

    psi_list = [0, 5, 10, 15]
    q_res = {}
    ev_res = {}
    for psi in psi_list:
        r, _ = _run("CD", psi_deg=psi, steps=300)
        q_res[psi] = max(r['q'])
        ev_res[psi] = r['volumetric_strain'][-1]
        err = abs(q_res[psi] - Q_ANALITICO) / Q_ANALITICO * 100
        print(f"  ψ={psi:2d}°: q={q_res[psi]:.2f} (err={err:.1f}%), εv={ev_res[psi]:+.6f}")

    # q constante
    q_vals = list(q_res.values())
    spread = (max(q_vals) - min(q_vals)) / np.mean(q_vals) * 100
    _check(spread < 2.0, "q spread < 2%", f"spread={spread:.2f}%")

    for psi in psi_list:
        err = abs(q_res[psi] - Q_ANALITICO) / Q_ANALITICO * 100
        _check(err < 2.0, f"q(ψ={psi}°) err < 2%", f"err={err:.2f}%")

    # εv dilata monotonicamente
    for psi in [5, 10, 15]:
        _check(ev_res[psi] < ev_res[0], f"εv(ψ={psi}°) < εv(ψ=0°)")
    _check(ev_res[10] < ev_res[5], "εv(10°) < εv(5°)")
    _check(ev_res[15] < ev_res[10], "εv(15°) < εv(10°)")


def _test_yield_surface():
    """Φ(σ_return) = 0 para vários ψ (incluindo ψ=0)."""
    print("\n" + "=" * 70)
    print("TESTE 6: Yield surface Φ(σ_return) = 0")
    print("=" * 70)

    for psi in [0, 5, 10, 15, 20]:
        r, model = _run("CD", psi_deg=psi, steps=200, eps_max=0.10)
        eig = np.sort(np.linalg.eigvalsh(model.stress))[::-1]
        c_f = model.plfun(model.equivalent_plastic_strain)
        Phi = (eig[0] - eig[2] + model.sign_conv * (eig[0] + eig[2]) * SIN_PHI
               - 2 * c_f * COS_PHI)
        print(f"  ψ={psi:2d}°: Φ = {Phi:.8f}")
        _check(abs(Phi) < 0.1, f"Φ≈0 para ψ={psi}°", f"Φ={Phi:.8f}")


def _test_hardening():
    """Hardening/softening com ψ=0 e ψ=10°."""
    print("\n" + "=" * 70)
    print("TESTE 7: Hardening / Softening")
    print("=" * 70)

    pairs = [[0.0, 20.0], [0.02, 40.0], [0.06, 10.0]]
    q_peak_exp = SIGMA3_REF * (KP - 1) + 2 * 40 * np.sqrt(KP)
    q_res_exp = SIGMA3_REF * (KP - 1) + 2 * 10 * np.sqrt(KP)

    for psi in [0, 10]:
        r, _ = _run("CD", psi_deg=psi, sampling_pairs=pairs, steps=300, eps_max=0.12)
        q = r['q']
        q_peak = max(q)
        q_final = q[-1]
        idx = q.index(q_peak)
        print(f"\n  ψ={psi}°: peak={q_peak:.1f} (exp≈{q_peak_exp:.1f}), "
              f"res={q_final:.1f} (exp≈{q_res_exp:.1f}), idx={idx}/{len(q)}")
        _check(q_peak > q_final, f"ψ={psi}°: peak > residual")
        _check(idx < len(q) - 10, f"ψ={psi}°: pico no meio")
        _check(abs(q_peak - q_peak_exp) / q_peak_exp < 0.10, f"ψ={psi}°: peak err<10%")
        _check(abs(q_final - q_res_exp) / q_res_exp < 0.15, f"ψ={psi}°: res err<15%")


def _test_cu_psi_neq0():
    """CU ψ≠0: verificação qualitativa (K_coupling empírico)."""
    print("\n" + "=" * 70)
    print("TESTE 8: CU ψ≠0 — qualitativo")
    print("=" * 70)
    print("  NOTA: K_coupling empírico — verificamos estabilidade, não valores exatos.\n")

    for psi in [0, 5, 10]:
        r, _ = _run("CU", psi_deg=psi)
        print(f"  ψ={psi:2d}°: q={max(r['q']):.1f}, u_final={r['pore_pressure'][-1]:.1f}")

    r0, _ = _run("CU", psi_deg=0)
    r10, _ = _run("CU", psi_deg=10)
    q0 = max(r0['q'])
    q10 = max(r10['q'])
    _check(q10 > q0 * 0.8, "q(ψ=10) > 80% q(ψ=0)")
    _check(q10 < 5000, "q(ψ=10) < 5000 (sem divergência)")
    _check(not np.isnan(q10), "q(ψ=10) finito")


def verificacao_numerica():
    """Executa todos os testes PASS/FAIL. Retorna True se todos passaram."""
    global _PASS, _FAIL, _TOTAL
    _PASS = _FAIL = _TOTAL = 0

    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  PARTE 1 — VERIFICAÇÃO NUMÉRICA                                ║")
    print("╠══════════════════════════════════════════════════════════════════╣")
    print(f"║  E={E:.0f}, ν={NU}, φ={PHI_DEG}°, c={C0}, σ3={SIGMA3_REF} kPa          ║")
    print(f"║  q_analítico = {Q_ANALITICO:.2f} kPa                                  ║")
    print("╚══════════════════════════════════════════════════════════════════╝")

    _test_cd_psi0()
    _test_cu_psi0()
    _test_uu()
    _test_cd_cu_relation()
    _test_cd_psi_neq0()
    _test_yield_surface()
    _test_hardening()
    _test_cu_psi_neq0()

    print("\n" + "=" * 70)
    print(f"RESULTADO: {_PASS}/{_TOTAL} PASS, {_FAIL}/{_TOTAL} FAIL")
    if _FAIL == 0:
        print(">>> TODOS OS TESTES PASSARAM <<<")
    else:
        print(f">>> {_FAIL} TESTE(S) FALHARAM <<<")
    print("=" * 70)
    return _FAIL == 0


# =====================================================================
#  PARTE 2 — VALIDAÇÃO GRÁFICA
# =====================================================================

def _plot_comparacao_ensaios(r_cd, r_cu, r_uu_list, q_teo):
    """Gráfico 2×3: q, u, UU, hardening, softening, Mohr."""
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    fig.suptitle(f"Comparação CD / CU / UU  (φ={PHI_DEG}°, c={C0} kPa)",
                 fontsize=13)

    # 1. q vs εa — CD, CU, UU
    ax = axes[0, 0]
    ax.plot(np.array(r_cd['axial_strain']) * 100, r_cd['q'], 'b-', lw=2, label='CD')
    ax.plot(np.array(r_cu['axial_strain']) * 100, r_cu['q'], 'r-', lw=2, label='CU')
    ax.plot(np.array(r_uu_list[0][1]['axial_strain']) * 100,
            r_uu_list[0][1]['q'], 'g-', lw=2, label='UU')
    ax.axhline(q_teo, color='b', ls='--', alpha=.5, label=f'CD teórico ({q_teo:.0f})')
    ax.set_xlabel('εa (%)')
    ax.set_ylabel('q (kPa)')
    ax.set_title('CD / CU / UU')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=.3)

    # 2. Poropressão CU
    ax = axes[0, 1]
    ax.plot(np.array(r_cu['axial_strain']) * 100, r_cu['pore_pressure'], 'r-', lw=2)
    ax.set_xlabel('εa (%)')
    ax.set_ylabel('u (kPa)')
    ax.set_title('CU — Poropressão')
    ax.grid(True, alpha=.3)

    # 3. UU multi-σ3
    ax = axes[0, 2]
    for s3, r in r_uu_list:
        ax.plot(np.array(r['axial_strain']) * 100, r['q'], lw=2,
                label=f'σ3={s3}')
    ax.set_xlabel('εa (%)')
    ax.set_ylabel('q (kPa)')
    ax.set_title('UU — cu independe de σ3')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=.3)

    # 4. Hardening
    ax = axes[1, 0]
    pairs_h = [[0.0, 20.0], [0.05, 50.0]]
    r_h, _ = _run("CD", psi_deg=0, sampling_pairs=pairs_h, steps=150, eps_max=0.15)
    ax.plot(np.array(r_h['axial_strain']) * 100, r_h['q'], 'b-', lw=2)
    ax.set_xlabel('εa (%)')
    ax.set_ylabel('q (kPa)')
    ax.set_title('Hardening (c: 20→50)')
    ax.grid(True, alpha=.3)

    # 5. Softening
    ax = axes[1, 1]
    pairs_s = [[0.0, 50.0], [0.01, 50.0], [0.05, 15.0]]
    r_s, _ = _run("CD", psi_deg=0, cohesion=50, sampling_pairs=pairs_s,
                   steps=150, eps_max=0.15)
    ax.plot(np.array(r_s['axial_strain']) * 100, r_s['q'], 'r-', lw=2)
    ax.set_xlabel('εa (%)')
    ax.set_ylabel('q (kPa)')
    ax.set_title('Softening (c: 50→15)')
    ax.grid(True, alpha=.3)

    # 6. Círculo de Mohr simplificado (CD)
    ax = axes[1, 2]
    theta = np.linspace(0, 2 * np.pi, 200)
    idx_pk = np.argmax(r_cd['q'])
    s1 = r_cd['sigma1'][idx_pk]
    s3 = r_cd['sigma3'][idx_pk]
    ct = (s1 + s3) / 2
    rd = (s1 - s3) / 2
    ax.plot(ct + rd * np.cos(theta), rd * np.sin(theta), 'b-', lw=2, label='Ruptura')
    sigma_env = np.linspace(0, s1 * 1.2, 100)
    tau_env = C0 + sigma_env * np.tan(np.radians(PHI_DEG))
    ax.plot(sigma_env, tau_env, 'k-', lw=2, label='Envoltória MC')
    ax.set_xlabel("σ' (kPa)")
    ax.set_ylabel('τ (kPa)')
    ax.set_title('Círculo de Mohr — CD')
    ax.set_aspect('equal')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=.3)
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)

    plt.tight_layout()
    plt.savefig('validacao_ensaios.png', dpi=150, bbox_inches='tight')
    print("  Salvo: validacao_ensaios.png")


def _plot_circulos_mohr(r_cd_list, r_cu_list, r_uu_list):
    """Círculos de Mohr: CD, CU (efetivo+total), UU (φu=0)."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    fig.suptitle("Círculos de Mohr Comparativos", fontsize=13)
    theta = np.linspace(0, 2 * np.pi, 200)
    colors = ['blue', 'green', 'red', 'orange', 'purple']

    def circle(ax, s1, s3, color, ls, label, lw=2, alpha=1.0):
        ct = (s1 + s3) / 2
        rd = (s1 - s3) / 2
        ax.plot(ct + rd * np.cos(theta), rd * np.sin(theta),
                color=color, ls=ls, lw=lw, label=label, alpha=alpha)
        return ct, rd

    def envelope(ax, smax, lbl=True, ls='-'):
        s = np.linspace(0, smax * 1.3, 100)
        tau = C0 + s * np.tan(np.radians(PHI_DEG))
        ax.plot(s, tau, color='k', ls=ls, lw=2.5,
                label=f"φ'={PHI_DEG}°, c'={C0}" if lbl else None)

    # --- CD ---
    ax = axes[0]
    ax.set_title('CD (Drenado)', fontsize=11)
    smax = 0
    for i, (s3v, r) in enumerate(r_cd_list):
        ip = np.argmax(r['q'])
        circle(ax, r['sigma1'][ip], r['sigma3'][ip], colors[i], '-',
               f'σ₃={s3v}')
        smax = max(smax, r['sigma1'][ip])
    envelope(ax, smax)
    ax.set_xlabel("σ' (kPa)")
    ax.set_ylabel('τ (kPa)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=.3)
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)
    ax.set_aspect('equal')

    # --- CU (efetivo + total) ---
    ax = axes[1]
    ax.set_title('CU (Não-Drenado)', fontsize=11)
    smax_eff = smax_tot = 0
    sm_tot, rd_tot = [], []
    for i, (s3v, r) in enumerate(r_cu_list):
        ip = np.argmax(r['q'])
        u = r['pore_pressure'][ip]
        s1e = r['sigma1'][ip]
        s3e = r['sigma3'][ip]
        s1t = s1e + u
        s3t = s3e + u
        smax_eff = max(smax_eff, s1e)
        smax_tot = max(smax_tot, s1t)
        sm_tot.append((s1t + s3t) / 2)
        rd_tot.append((s1t - s3t) / 2)
        circle(ax, s1t, s3t, colors[i], '-', f"σ₃={s3v} (total)", lw=2)
        circle(ax, s1e, s3e, colors[i], '--', f"σ₃={s3v} (σ')", lw=1.5, alpha=.7)
    # Envoltória efetiva (tracejada)
    envelope(ax, max(smax_eff, smax_tot), ls='--')
    # Envoltória de tensões totais (regressão linear: R = c_cu·cosφ + σm·sinφ)
    sm_arr = np.array(sm_tot)
    rd_arr = np.array(rd_tot)
    n_pts = len(sm_arr)
    b_reg = (n_pts * np.sum(sm_arr * rd_arr) - np.sum(sm_arr) * np.sum(rd_arr)) / \
            (n_pts * np.sum(sm_arr**2) - np.sum(sm_arr)**2)
    a_reg = (np.sum(rd_arr) - b_reg * np.sum(sm_arr)) / n_pts
    phi_cu_rad = np.arcsin(np.clip(b_reg, -1, 1))
    c_cu = a_reg / np.cos(phi_cu_rad) if abs(np.cos(phi_cu_rad)) > 1e-10 else 0
    phi_cu_deg = np.degrees(phi_cu_rad)
    s_env = np.linspace(0, max(smax_eff, smax_tot) * 1.3, 100)
    tau_tot = c_cu + s_env * np.tan(phi_cu_rad)
    ax.plot(s_env, tau_tot, 'k-', lw=2,
            label=f"Env. total (φ_cu={phi_cu_deg:.0f}°, c_cu={c_cu:.0f})")
    ax.set_xlabel("σ (kPa)")
    ax.set_ylabel('τ (kPa)')
    ax.legend(fontsize=7, loc='upper left')
    ax.grid(True, alpha=.3)
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)
    ax.set_aspect('equal')

    # --- UU ---
    ax = axes[2]
    ax.set_title('UU  (φᵤ = 0)', fontsize=11)
    smax = 0
    for i, (s3v, r) in enumerate(r_uu_list):
        ip = np.argmax(r['q'])
        u = r['pore_pressure'][ip]
        s1t = r['sigma1'][ip] + u
        s3t = r['sigma3'][ip] + u
        smax = max(smax, s1t)
        circle(ax, s1t, s3t, colors[i], '-', f'σ₃={s3v}', lw=2)
        if i == 0:
            circle(ax, r['sigma1'][ip], r['sigma3'][ip], 'black', '--',
                   "σ' (efetivo)", lw=1.5)
    cu = max(r_uu_list[0][1]['q']) / 2
    ax.axhline(cu, color='k', lw=2.5, label=f'cu = {cu:.0f} kPa')
    envelope(ax, smax)
    ax.set_xlabel("σ (kPa)")
    ax.set_ylabel('τ (kPa)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=.3)
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig('circulos_mohr.png', dpi=150, bbox_inches='tight')
    print("  Salvo: circulos_mohr.png")


def _plot_dilatancia():
    """Efeito puro da dilatância (ψ) em CD e CU.

    Linha 1 — CD: q ≈ constante (independe de ψ), εv dilata com ψ.
    Linha 2 — CU: q cresce e u fica negativo com ψ (sem estado crítico
    no MC com ψ constante — σ₃' cresce sem limite).
    """
    psi_vals = [0, 2, 5]
    labels = ['ψ=0° (NC)', 'ψ=2° (OC leve)', 'ψ=5° (OC moderado)']
    colors = ['blue', 'green', 'red']

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f'Efeito da Dilatância (ψ) — CD vs CU  —  σ₃={SIGMA3_REF} kPa',
                 fontsize=12)

    # --- CD: ψ varia ---
    cd_data = []
    for psi, lbl, col in zip(psi_vals, labels, colors):
        r_cd, _ = _run("CD", psi_deg=psi, steps=200, eps_max=0.10)
        cd_data.append((lbl, col, r_cd))

    ax = axes[0, 0]
    for lbl, col, r in cd_data:
        ax.plot(np.array(r['axial_strain']) * 100, r['q'], color=col, lw=2, label=lbl)
    ax.set_xlabel('εa (%)')
    ax.set_ylabel('q (kPa)')
    ax.set_title('CD: Tensão desviadora (q ≈ const)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=.3)

    ax = axes[0, 1]
    for lbl, col, r in cd_data:
        ax.plot(np.array(r['axial_strain']) * 100,
                np.array(r['volumetric_strain']) * 100,
                color=col, lw=2, label=lbl)
    ax.axhline(0, color='k', ls='--', alpha=.5)
    ax.set_xlabel('εa (%)')
    ax.set_ylabel('εv (%)')
    ax.set_title('CD: Deformação volumétrica (ψ → dilatação)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=.3)

    # --- CU: ψ varia (sem softening) ---
    # NOTA: MC + ψ constante + CU não tem estado crítico.
    # σ₃' cresce sem limite → u cada vez mais negativo → q cresce.
    # Isso é uma LIMITAÇÃO do modelo, não comportamento real.
    cu_data = []
    for psi, lbl, col in zip(psi_vals, labels, colors):
        r_cu, _ = _run("CU", psi_deg=psi, steps=200, eps_max=0.10)
        cu_data.append((lbl, col, r_cu))

    ax = axes[1, 0]
    for lbl, col, r in cu_data:
        ax.plot(np.array(r['axial_strain']) * 100, r['q'], color=col, lw=2, label=lbl)
    ax.set_xlabel('εa (%)')
    ax.set_ylabel('q (kPa)')
    ax.set_title('CU: Tensão desviadora (q cresce sem limite*)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=.3)

    ax = axes[1, 1]
    for lbl, col, r in cu_data:
        ax.plot(np.array(r['axial_strain']) * 100, r['pore_pressure'],
                color=col, lw=2, label=lbl)
    ax.axhline(0, color='k', ls='--', lw=2, alpha=.7)
    ax.set_xlabel('εa (%)')
    ax.set_ylabel('u (kPa)')
    ax.set_title('CU: Poropressão (u → negativo com ψ*)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=.3)

    # nota de rodapé
    fig.text(0.5, -0.01,
             '* Limitação do MC com ψ constante em CU: sem estado crítico '
             '(σ₃\' cresce sem limite). Ver comportamento_OC.png para pico+queda via softening.',
             ha='center', fontsize=8, style='italic', color='gray')

    plt.tight_layout()
    plt.savefig('comparacao_dilatancia.png', dpi=150, bbox_inches='tight')
    print("  Salvo: comparacao_dilatancia.png")


def _plot_comportamento_OC():
    """CU: Comportamento NC vs OC via softening de coesão (ψ=0).

    O modelo MC com ψ constante não tem estado crítico em CU (σ₃' cresce
    sem limite). O comportamento realista de argila OC (pico + queda) é
    capturado via softening de coesão c(εp): coesão alta inicial que
    decai até o valor residual.

    Combinar softening + ψ>0 em CU NÃO funciona: o efeito de ψ domina
    completamente (q e u idênticos com ou sem softening quando ψ>0).
    """
    cu_configs = [
        ('NC (c=20 const)',       'blue',  None, C0),
        ('OC leve (c: 30→20)',    'green',
         [[0.0, 30.0], [0.01, 30.0], [0.05, 20.0]], 30.0),
        ('OC forte (c: 40→20)',   'red',
         [[0.0, 40.0], [0.01, 40.0], [0.05, 20.0]], 40.0),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f'CU: Comportamento NC vs OC (Softening de Coesão, ψ=0)  —  '
                 f'σ₃={SIGMA3_REF} kPa', fontsize=12)

    oc_data = []
    for lbl, col, pairs, c0 in cu_configs:
        r, _ = _run("CU", psi_deg=0, sampling_pairs=pairs,
                     cohesion=c0, steps=200, eps_max=0.10)
        oc_data.append((lbl, col, r))

    # q vs εa — pico + queda
    ax = axes[0]
    for lbl, col, r in oc_data:
        ax.plot(np.array(r['axial_strain']) * 100, r['q'],
                color=col, lw=2, label=lbl)
    ax.set_xlabel('εa (%)')
    ax.set_ylabel('q (kPa)')
    ax.set_title('Tensão desviadora (pico + queda em OC)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=.3)

    # u vs εa
    ax = axes[1]
    for lbl, col, r in oc_data:
        ax.plot(np.array(r['axial_strain']) * 100, r['pore_pressure'],
                color=col, lw=2, label=lbl)
    ax.axhline(0, color='k', ls='--', lw=2, alpha=.7)
    ax.set_xlabel('εa (%)')
    ax.set_ylabel('u (kPa)')
    ax.set_title('Poropressão (u convergem ao residual)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=.3)

    plt.tight_layout()
    plt.savefig('comportamento_OC.png', dpi=150, bbox_inches='tight')
    print("  Salvo: comportamento_OC.png")


def validacao_grafica(show=False):
    """Gera todos os gráficos de validação."""
    if not show:
        matplotlib.use('Agg')
    print("\n╔══════════════════════════════════════════════════════════════════╗")
    print("║  PARTE 2 — VALIDAÇÃO GRÁFICA                                    ║")
    print("╚══════════════════════════════════════════════════════════════════╝\n")

    # CD, CU com múltiplos σ3
    cd_list = []
    cu_list = []
    for s3 in SIGMA3_LIST:
        r_cd, _ = _run("CD", sigma3=s3, steps=150, eps_max=0.10)
        r_cu, _ = _run("CU", sigma3=s3, steps=150, eps_max=0.10)
        cd_list.append((s3, r_cd))
        cu_list.append((s3, r_cu))

    # UU multi-σ3
    uu_list = []
    for s3 in [50, 100, 200]:
        r_uu, _ = _run("UU", sigma3=s3, steps=100, eps_max=0.10)
        uu_list.append((s3, r_uu))

    # σ3 de referência para comparação geral
    idx_ref = SIGMA3_LIST.index(SIGMA3_REF) if SIGMA3_REF in SIGMA3_LIST else 1
    r_cd_ref = cd_list[idx_ref][1]
    r_cu_ref = cu_list[idx_ref][1]
    q_teo_ref = SIGMA3_LIST[idx_ref] * (KP - 1) + 2 * C0 * np.sqrt(KP)

    print("Gerando gráficos...")
    _plot_comparacao_ensaios(r_cd_ref, r_cu_ref, uu_list, q_teo_ref)
    _plot_circulos_mohr(cd_list, cu_list, uu_list)
    _plot_dilatancia()
    _plot_comportamento_OC()
    print("Gráficos concluídos.\n")
    if show:
        plt.show()


# =====================================================================
#  MAIN
# =====================================================================

def main():
    """Executa validação completa (numérica + gráfica)."""

    args = sys.argv[1:]

    do_check = True
    do_plot = True
    do_show = '--show' in args

    if '--check' in args:
        do_plot = False
    elif '--plot' in args:
        do_check = False

    print("\n" + "=" * 70)
    print("  VALIDAÇÃO DO MODELO MOHR-COULOMB — ENSAIOS TRIAXIAIS")
    print("=" * 70)
    print(f"  E={E:.0f} kPa, ν={NU}, φ={PHI_DEG}°, c={C0} kPa")
    print(f"  σ3_ref={SIGMA3_REF} kPa,  q_analítico={Q_ANALITICO:.2f} kPa")
    print()

    ok = True
    if do_check:
        ok = verificacao_numerica()

    if do_plot:
        validacao_grafica(show=do_show)

    print("\n" + "=" * 70)
    if do_check and ok:
        print("  VALIDAÇÃO CONCLUÍDA COM SUCESSO!")
    elif do_check:
        print("  VALIDAÇÃO CONCLUÍDA — VERIFICAR FALHAS ACIMA")
    else:
        print("  GRÁFICOS GERADOS COM SUCESSO!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
