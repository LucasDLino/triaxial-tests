"""
Exemplos Completos e Demonstrações Avançadas

Contém exemplos detalhados com hardening, dilatância, círculos de Mohr e análises comparativas
"""

import numpy as np
import matplotlib.pyplot as plt
from mohr_courlomb import MohrCoulombModel
from triaxial import TriaxialTest

# Configuração de estilo para os gráficos
plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10
plt.rcParams['lines.linewidth'] = 2


def plot_triaxial_results(tests_data, title="Resultados dos Ensaios Triaxiais", save_path=None):
    """
    Plota resultados de ensaios triaxiais
    
    Parameters:
    -----------
    tests_data : list of dict
        Lista com dicionários contendo:
        - 'label': Nome do teste
        - 'eps': Deformação axial
        - 'q': Tensão desviadora
        - 'p': Tensão média
        - 'epsv': Deformação volumétrica
        - 'color': Cor da linha (opcional)
    title : str
        Título geral dos gráficos
    save_path : str, optional
        Caminho para salvar a figura
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(title, fontsize=16, fontweight='bold')
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    for idx, data in enumerate(tests_data):
        label = data['label']
        eps = np.array(data['eps']) * 100  # Converter para %
        q = np.array(data['q'])
        p = np.array(data['p'])
        epsv = np.array(data['epsv']) * 100  # Converter para %
        color = data.get('color', colors[idx % len(colors)])
        
        # Gráfico 1: q vs εa
        axes[0, 0].plot(eps, q, label=label, color=color, marker='o', 
                        markevery=len(eps)//10, markersize=4)
        
        # Gráfico 2: εv vs εa
        axes[0, 1].plot(eps, epsv, label=label, color=color, marker='s',
                        markevery=len(eps)//10, markersize=4)
        
        # Gráfico 3: Trajetória de tensões (p-q)
        axes[1, 0].plot(p, q, label=label, color=color, marker='^',
                        markevery=len(eps)//10, markersize=4)
        
        # Gráfico 4: p vs εa
        axes[1, 1].plot(eps, p, label=label, color=color, marker='d',
                        markevery=len(eps)//10, markersize=4)
    
    # Configurar gráfico 1
    axes[0, 0].set_xlabel('Deformação Axial ε$_a$ (%)', fontweight='bold')
    axes[0, 0].set_ylabel('Tensão Desviadora q (kPa)', fontweight='bold')
    axes[0, 0].set_title('Curva Tensão-Deformação', fontweight='bold')
    axes[0, 0].legend(loc='best', framealpha=0.9)
    axes[0, 0].grid(True, alpha=0.3)
    
    # Configurar gráfico 2
    axes[0, 1].set_xlabel('Deformação Axial ε$_a$ (%)', fontweight='bold')
    axes[0, 1].set_ylabel('Deformação Volumétrica ε$_v$ (%)', fontweight='bold')
    axes[0, 1].set_title('Variação Volumétrica', fontweight='bold')
    axes[0, 1].legend(loc='best', framealpha=0.9)
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].axhline(y=0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
    
    # Configurar gráfico 3
    axes[1, 0].set_xlabel('Tensão Média p (kPa)', fontweight='bold')
    axes[1, 0].set_ylabel('Tensão Desviadora q (kPa)', fontweight='bold')
    axes[1, 0].set_title('Trajetória de Tensões', fontweight='bold')
    axes[1, 0].legend(loc='best', framealpha=0.9)
    axes[1, 0].grid(True, alpha=0.3)
    
    # Configurar gráfico 4
    axes[1, 1].set_xlabel('Deformação Axial ε$_a$ (%)', fontweight='bold')
    axes[1, 1].set_ylabel('Tensão Média p (kPa)', fontweight='bold')
    axes[1, 1].set_title('Evolução da Tensão Média', fontweight='bold')
    axes[1, 1].legend(loc='best', framealpha=0.9)
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"\n✓ Gráfico salvo: {save_path}")
    
    plt.show()


def plot_hardening_analysis(eps, q, p, epsv, model, save_path=None):
    """
    Plota análise detalhada incluindo evolução de hardening
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Análise Detalhada com Hardening/Softening', fontsize=16, fontweight='bold')
    
    eps_pct = np.array(eps) * 100
    
    # Gráfico 1: q vs εa com marcadores de pico
    axes[0, 0].plot(eps_pct, q, 'b-', linewidth=2, label='q')
    q_max_idx = np.argmax(q)
    axes[0, 0].plot(eps_pct[q_max_idx], q[q_max_idx], 'r*', markersize=15, 
                    label=f'Pico: q={q[q_max_idx]:.2f} kPa')
    axes[0, 0].set_xlabel('Deformação Axial ε$_a$ (%)', fontweight='bold')
    axes[0, 0].set_ylabel('Tensão Desviadora q (kPa)', fontweight='bold')
    axes[0, 0].set_title('Tensão Desviadora', fontweight='bold')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend()
    
    # Gráfico 2: Curva de hardening
    if hasattr(model, 'sampling_pairs') and len(model.sampling_pairs) > 2:
        eps_p_values = model.sampling_pairs[:-1, 0] * 100
        c_values = model.sampling_pairs[:-1, 1]
        axes[0, 1].plot(eps_p_values, c_values, 'ro-', linewidth=2, 
                        markersize=8, label='Curva de Hardening')
        axes[0, 1].plot(model.equivalent_plastic_strain * 100, 
                        model.plfun(model.equivalent_plastic_strain), 
                        'g*', markersize=15, label='Estado Final')
        axes[0, 1].set_xlabel('Deformação Plástica Equivalente ε$_p$ (%)', fontweight='bold')
        axes[0, 1].set_ylabel('Coesão c (kPa)', fontweight='bold')
        axes[0, 1].set_title('Evolução de Hardening/Softening', fontweight='bold')
        axes[0, 1].grid(True, alpha=0.3)
        axes[0, 1].legend()
    else:
        axes[0, 1].text(0.5, 0.5, 'Modelo Perfeitamente Plástico\n(sem hardening)', 
                        ha='center', va='center', fontsize=12, transform=axes[0, 1].transAxes)
        axes[0, 1].set_title('Evolução de Hardening/Softening', fontweight='bold')
    
    # Gráfico 3: Trajetória de tensões com envoltória
    axes[1, 0].plot(p, q, 'b-', linewidth=2, label='Trajetória')
    axes[1, 0].plot(p[0], q[0], 'go', markersize=10, label='Início')
    axes[1, 0].plot(p[-1], q[-1], 'rs', markersize=10, label='Fim')
    
    # Envoltória de Mohr-Coulomb
    p_range = np.linspace(0, max(p) * 1.5, 100)
    phi = model.phi
    c_final = model.plfun(model.equivalent_plastic_strain)
    M = 6 * np.sin(phi) / (3 - np.sin(phi))
    q_envelope = M * p_range + 2 * np.sqrt(6) * c_final * np.cos(phi) / (3 - np.sin(phi))
    axes[1, 0].plot(p_range, q_envelope, 'r--', linewidth=1.5, alpha=0.7, 
                    label=f'Envoltória MC (c={c_final:.1f} kPa)')
    
    axes[1, 0].set_xlabel('Tensão Média p (kPa)', fontweight='bold')
    axes[1, 0].set_ylabel('Tensão Desviadora q (kPa)', fontweight='bold')
    axes[1, 0].set_title('Trajetória de Tensões e Envoltória', fontweight='bold')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend()
    axes[1, 0].set_xlim(left=0)
    axes[1, 0].set_ylim(bottom=0)
    
    # Gráfico 4: Deformações
    axes[1, 1].plot(eps_pct, np.array(epsv)*100, 'r-', linewidth=2, label='ε$_v$ (volumétrica)')
    axes[1, 1].set_xlabel('Deformação Axial ε$_a$ (%)', fontweight='bold')
    axes[1, 1].set_ylabel('Deformação Volumétrica ε$_v$ (%)', fontweight='bold')
    axes[1, 1].set_title('Evolução das Deformações', fontweight='bold')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].axhline(y=0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
    axes[1, 1].legend()
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"\n✓ Gráfico salvo: {save_path}")
    
    plt.show()


def plot_mohr_circles(tests_data, title="Círculos de Mohr", save_path=None):
    """Plota círculos de Mohr para diferentes ensaios"""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    max_sigma = 0
    min_sigma = 0
    max_tau = 0
    
    for idx, data in enumerate(tests_data):
        label = data['label']
        sigma1 = data['sigma1']
        sigma3 = data['sigma3']
        color = data.get('color', colors[idx % len(colors)])
        
        center = (sigma1 + sigma3) / 2
        radius = (sigma1 - sigma3) / 2
        
        theta = np.linspace(0, np.pi, 100)
        sigma_n = center + radius * np.cos(theta)
        tau = radius * np.sin(theta)
        
        ax.plot(sigma_n, tau, linewidth=2.5, label=label, color=color)
        ax.plot([sigma3, sigma1], [0, 0], 'o', markersize=8, color=color)
        
        max_sigma = max(max_sigma, sigma1)
        min_sigma = min(min_sigma, sigma3)
        max_tau = max(max_tau, radius)
    
    # Adicionar envoltória de Mohr-Coulomb
    if 'phi' in tests_data[0] and 'cohesion' in tests_data[0]:
        phi = tests_data[0]['phi']
        c = tests_data[0]['cohesion']
        
        sigma_range = np.linspace(min_sigma * 0.1, max_sigma * 1.3, 100)
        tau_envelope = c + sigma_range * np.tan(phi)
        
        ax.plot(sigma_range, tau_envelope, 'r--', linewidth=2.5, 
                label=f'Envoltória MC (c={c:.1f} kPa, φ={np.degrees(phi):.1f}°)', alpha=0.7)
    
    ax.set_xlabel('Tensão Normal σ (kPa)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Tensão Cisalhante τ (kPa)', fontweight='bold', fontsize=12)
    ax.set_title(title, fontweight='bold', fontsize=14)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    
    x_margin = (max_sigma - min_sigma) * 0.1
    y_margin = max_tau * 0.1
    ax.set_xlim(min_sigma - x_margin, max_sigma + x_margin)
    ax.set_ylim(-y_margin, max_tau + y_margin)
    ax.set_aspect('equal', adjustable='datalim')
    
    ax.legend(loc='upper left', framealpha=0.95, fontsize=10)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"\n✓ Gráfico salvo: {save_path}")
    
    plt.show()


def exemplo_comparacao_tipos():
    """Exemplo: Comparação entre CD, CU e UU"""
    print("="*70)
    print("EXEMPLO 1: Comparação de Tipos de Ensaio (CD, CU, UU)")
    print("="*70)
    
    E = 50000  # kPa
    nu = 0.3
    phi_deg = 30
    cohesion = 20
    psi_deg = 0  # NÃO-ASSOCIADO: ψ = 0
    sigma3 = 100
    
    tests_data = []
    
    # CD
    model_CD = MohrCoulombModel(E=E, nu=nu, phi_deg=phi_deg, cohesion=cohesion, psi_deg=psi_deg)
    test_CD = TriaxialTest(model_CD, sigma3=sigma3, test_type="CD")
    results_CD = test_CD.run()
    
    # CU
    model_CU = MohrCoulombModel(E=E, nu=nu, phi_deg=phi_deg, cohesion=cohesion, psi_deg=psi_deg)
    test_CU = TriaxialTest(model_CU, sigma3=sigma3, test_type="CU")
    results_CU = test_CU.run()
    
    # UU
    model_UU = MohrCoulombModel(E=E, nu=nu, phi_deg=phi_deg, cohesion=cohesion, psi_deg=psi_deg)
    test_UU = TriaxialTest(model_UU, sigma3=sigma3, test_type="UU")
    results_UU = test_UU.run()
    
    tests_data = [
        {'label': 'CD', 'eps': results_CD['axial_strain'], 'q': results_CD['q'], 
         'p': results_CD['p'], 'epsv': results_CD['volumetric_strain'], 'color': '#1f77b4'},
        {'label': 'CU', 'eps': results_CU['axial_strain'], 'q': results_CU['q'], 
         'p': results_CU['p'], 'epsv': results_CU['volumetric_strain'], 'color': '#ff7f0e'},
        {'label': 'UU', 'eps': results_UU['axial_strain'], 'q': results_UU['q'], 
         'p': results_UU['p'], 'epsv': results_UU['volumetric_strain'], 'color': '#2ca02c'},
    ]
    
    plot_triaxial_results(tests_data, 
                          title="Comparação: CD vs CU vs UU",
                          save_path="comparacao_tipos_ensaio.png")


def exemplo_hardening():
    """Exemplo: Hardening Piecewise Linear"""
    print("\n" + "="*70)
    print("EXEMPLO 2: Modelo com Hardening Piecewise Linear")
    print("="*70)
    
    E = 50000
    nu = 0.3
    phi_deg = 30
    psi_deg = 0  # NÃO-ASSOCIADO: ψ = 0
    sigma3 = 100
    
    sampling_pairs = [
        [0.0, 10],
        [0.02, 25],
        [0.05, 30],
        [0.10, 25],
    ]
    
    model_CD = MohrCoulombModel(
        E=E, nu=nu, phi_deg=phi_deg, cohesion=10, psi_deg=psi_deg,
        sampling_pairs=sampling_pairs
    )
    test_CD = TriaxialTest(model_CD, sigma3=sigma3, test_type="CD")
    results = test_CD.run(eps_max=0.20, steps=300)
    
    plot_hardening_analysis(results['axial_strain'], results['q'], results['p'], 
                           results['volumetric_strain'], model_CD,
                           save_path="analise_hardening.png")


def exemplo_dilatancia():
    """Exemplo: Efeito da Dilatância"""
    print("\n" + "="*70)
    print("EXEMPLO 3: Efeito da Dilatância")
    print("="*70)
    
    E = 50000
    nu = 0.3
    phi_deg = 30
    cohesion = 20
    sigma3 = 100
    
    tests_data = []
    psi_angles = [0, 5, 10, 15]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    for idx, psi_deg in enumerate(psi_angles):
        model = MohrCoulombModel(E=E, nu=nu, phi_deg=phi_deg, cohesion=cohesion, psi_deg=psi_deg)
        test = TriaxialTest(model, sigma3=sigma3, test_type="CD")
        results = test.run()
        
        tests_data.append({
            'label': f'ψ = {psi_deg}°',
            'eps': results['axial_strain'],
            'q': results['q'],
            'p': results['p'],
            'epsv': results['volumetric_strain'],
            'color': colors[idx]
        })
    
    plot_triaxial_results(tests_data,
                          title="Efeito da Dilatância em Ensaios CD",
                          save_path="efeito_dilatancia.png")


if __name__ == "__main__":
    exemplo_comparacao_tipos()
    exemplo_hardening()
    exemplo_dilatancia()
    
    print("\n" + "="*70)
    print("✓ Todos os exemplos avançados concluídos!")
    print("="*70)
