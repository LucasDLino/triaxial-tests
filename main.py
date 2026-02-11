"""
Main - Ensaios Triaxiais com Modelo de Mohr-Coulomb

Demonstra os três tipos principais de ensaio com validação analítica.

MÓDULOS DISPONÍVEIS:
-------------------
- main.py (este arquivo)     : Exemplos básicos CD, CU, UU
- exemplos_hardening.py      : Hardening e softening
- analise_mohr.py           : Círculos de Mohr e ângulo de atrito
- validacao_analitica.py    : Validações numéricas

Para executar:
  python main.py                    # Básicos
  python exemplos_hardening.py      # Hardening/Softening
  python analise_mohr.py            # Círculos de Mohr
"""

import numpy as np
import matplotlib.pyplot as plt
from mohr_courlomb import MohrCoulombModel
from triaxial import TriaxialTest
from validacao_analitica import validar_ensaio_CD, validar_ensaio_CU, validar_ensaio_UU

plt.style.use('seaborn-v0_8-darkgrid')

# ============================================================================
# PARÂMETROS GLOBAIS
# ============================================================================
PARAMS = {
    'E': 25000,         # kPa - Módulo de Young
    'nu': 0.30,         # Coeficiente de Poisson
    'phi_deg': 25,      # graus - Ângulo de atrito
    'cohesion': 35,     # kPa - Coesão
    'psi_deg': 0,       # graus - Dilatância (NÃO-ASSOCIADO: ψ = 0)
    'sigma3': 30,       # kPa - Tensão confinante
    'eps_max': 0.15,    # Deformação máxima
    'steps': 200        # Número de passos
}


def criar_modelo():
    """Cria modelo Mohr-Coulomb com parâmetros globais"""
    return MohrCoulombModel(
        E=PARAMS['E'], nu=PARAMS['nu'], 
        phi_deg=PARAMS['phi_deg'], cohesion=PARAMS['cohesion'],
        psi_deg=PARAMS['psi_deg']
    )


def rodar_ensaio(test_type):
    """Executa ensaio triaxial"""
    model = criar_modelo()
    test = TriaxialTest(model, sigma3=PARAMS['sigma3'], test_type=test_type)
    return test.run(eps_max=PARAMS['eps_max'], steps=PARAMS['steps'])


def print_parametros():
    """Imprime parâmetros do modelo"""
    phi = PARAMS['phi_deg']
    c = PARAMS['cohesion']
    s3 = PARAMS['sigma3']
    
    Kp = (1 + np.sin(np.radians(phi))) / (1 - np.sin(np.radians(phi)))
    q_teo = s3 * (Kp - 1) + 2 * c * np.sqrt(Kp)
    p_apex = c / np.tan(np.radians(phi))
    
    print(f"  E={PARAMS['E']} kPa, ν={PARAMS['nu']}, φ={phi}°, c={c} kPa, ψ={PARAMS['psi_deg']}°")
    print(f"  σ₃={s3} kPa, q_teórico={q_teo:.1f} kPa, p_apex={p_apex:.1f} kPa")


def plot_resultado(results, title, save_path):
    """Plota resultado de ensaio triaxial"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(title, fontsize=13, fontweight='bold')
    
    eps = np.array(results['axial_strain']) * 100
    q = np.array(results['q'])
    p = np.array(results['p'])
    
    # q vs εa
    axes[0].plot(eps, q, 'b-', lw=2.5)
    idx = np.argmax(q)
    axes[0].plot(eps[idx], q[idx], 'r*', ms=15, label=f'Pico: {q[idx]:.1f} kPa')
    axes[0].set_xlabel('Deformação Axial εₐ (%)')
    axes[0].set_ylabel('Tensão Desviadora q (kPa)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Trajetória p-q
    axes[1].plot(p, q, 'g-', lw=2.5)
    axes[1].plot(p[0], q[0], 'go', ms=10, label='Início')
    axes[1].plot(p[-1], q[-1], 'rs', ms=10, label='Fim')
    axes[1].set_xlabel("Tensão Média p' (kPa)")
    axes[1].set_ylabel('Tensão Desviadora q (kPa)')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"  ✓ {save_path}")
    plt.close()


def exemplo_CD():
    """Ensaio CD - Consolidated Drained"""
    print("\n" + "="*60)
    print("ENSAIO CD (Consolidated Drained)")
    print("="*60)
    print_parametros()
    
    results = rodar_ensaio('CD')
    q_max = max(results['q'])
    print(f"\n  Resultado: q_max = {q_max:.2f} kPa")
    
    # Validação
    idx = np.argmax(results['q'])
    val = validar_ensaio_CD(PARAMS['sigma3'], PARAMS['phi_deg'], PARAMS['cohesion'],
                            results['q'], results['p'], results['sigma1'][idx])
    
    plot_resultado(results, f"Ensaio CD: σ₃={PARAMS['sigma3']} kPa", "ensaio_CD.png")
    return results, val


def exemplo_CU():
    """Ensaio CU - Consolidated Undrained"""
    print("\n" + "="*60)
    print("ENSAIO CU (Consolidated Undrained)")
    print("="*60)
    print_parametros()
    
    results = rodar_ensaio('CU')
    q_max = max(results['q'])
    u_max = max(results['pore_pressure'])
    epsv_max = max(abs(np.array(results['volumetric_strain'])))
    
    print(f"\n  Resultados: q_max={q_max:.2f} kPa, u_max={u_max:.2f} kPa")
    print(f"  Verificação: |εv|_max = {epsv_max:.2e} (deve ser ≈0)")
    
    val = validar_ensaio_CU(results['volumetric_strain'])
    
    plot_resultado(results, f"Ensaio CU: σ₃_total={PARAMS['sigma3']} kPa", "ensaio_CU.png")
    
    # Gráfico de poropressão
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(np.array(results['axial_strain'])*100, results['pore_pressure'], 'c-', lw=2)
    ax.set_xlabel('Deformação Axial εₐ (%)')
    ax.set_ylabel('Poropressão u (kPa)')
    ax.set_title(f'Poropressão CU (u_max = {u_max:.1f} kPa)')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("poropressao_CU.png", dpi=150)
    print(f"  ✓ poropressao_CU.png")
    plt.close()
    
    return results, val


def exemplo_UU():
    """Ensaio UU - Unconsolidated Undrained"""
    print("\n" + "="*60)
    print("ENSAIO UU (Unconsolidated Undrained)")
    print("="*60)
    print_parametros()
    
    results = rodar_ensaio('UU')
    q_max = max(results['q'])
    epsv_max = max(abs(np.array(results['volumetric_strain'])))
    
    print(f"\n  Resultado: q_max={q_max:.2f} kPa, cu={q_max/2:.2f} kPa")
    print(f"  Verificação: |εv|_max = {epsv_max:.2e}")
    
    val = validar_ensaio_UU(results['volumetric_strain'])
    
    plot_resultado(results, f"Ensaio UU: σ₃={PARAMS['sigma3']} kPa", "ensaio_UU.png")
    return results, val


# ============================================================================
# EXECUÇÃO PRINCIPAL
# ============================================================================
if __name__ == "__main__":
    print("\n" + "█"*60)
    print("█  ENSAIOS TRIAXIAIS - MODELO DE MOHR-COULOMB".ljust(59) + "█")
    print("█"*60)
    
    # Executar ensaios
    _, val_CD = exemplo_CD()
    _, val_CU = exemplo_CU()
    _, val_UU = exemplo_UU()
    
    # Resumo
    print("\n" + "="*60)
    print("RESUMO DAS VALIDAÇÕES")
    print("="*60)
    for nome, val in [('CD', val_CD), ('CU', val_CU), ('UU', val_UU)]:
        status = '✓ OK' if val['aprovado'] else '✗ FALHA'
        print(f"  {nome}: {status}")
    print("="*60)
    
    print("\nPróximos passos:")
    print("  python exemplos_hardening.py  → Hardening e Softening")
    print("  python analise_mohr.py        → Círculos de Mohr")
    print("="*60 + "\n")
