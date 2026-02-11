"""
Exemplo demonstrando as melhorias no ensaio CU
Comparação com implementação Julia (Cam Clay)

Principais ajustes implementados:
1. Rastreamento de poropressão (u)
2. Distinção entre tensões totais e efetivas
3. Relações normalizadas (q/σ3, u/σ3)
4. Compatibilidade retroativa mantida
"""

import numpy as np
import matplotlib.pyplot as plt
from mohr_courlomb import MohrCoulombModel
from triaxial import TriaxialTest

# Configuração do modelo
E = 30000  # kPa
nu = 0.3
phi_deg = 30
cohesion = 10  # kPa
psi_deg = 0  # NÃO-ASSOCIADO: ψ = 0 (mais realista para solos)

model = MohrCoulombModel(E, nu, phi_deg, cohesion, psi_deg)

# Executar ensaio CU
sigma3_total = 100  # kPa (pressão de câmara - constante)
test_CU = TriaxialTest(model, sigma3=sigma3_total, test_type="CU")

# NOVA FORMA: Acesso por dicionário
results = test_CU.run(eps_max=0.15, steps=200)

# Acessar novos dados
eps_a = np.array(results['axial_strain']) * 100  # Converter para %
q = np.array(results['q'])
p_prime = np.array(results['p'])
sigma1_total = np.array(results['sigma1'])
sigma3_total_array = np.array(results['sigma3'])
u = np.array(results['pore_pressure'])
q_sigma3 = np.array(results['q_sigma3_ratio'])
u_sigma3 = np.array(results['u_sigma3_ratio'])

# OU FORMA ANTIGA: Desempacotamento como tupla (compatibilidade)
model2 = MohrCoulombModel(E, nu, phi_deg, cohesion, psi_deg)
test_CU2 = TriaxialTest(model2, sigma3=sigma3_total, test_type="CU")
eps, q_old, p_old, epsv = test_CU2.run(eps_max=0.15, steps=200)
print("✓ Compatibilidade retroativa mantida!")

# Plotar resultados (estilo Julia)
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('Ensaio CU - Implementação Ajustada\n'
             f'σ3_total = {sigma3_total} kPa, φ = {phi_deg}°', 
             fontsize=14, fontweight='bold')

# 1. q vs εa
axes[0, 0].plot(eps_a, q, 'b-', linewidth=2)
axes[0, 0].set_xlabel('Deformação Axial εa (%)')
axes[0, 0].set_ylabel('Tensão Desviadora q (kPa)')
axes[0, 0].grid(True, alpha=0.3)
axes[0, 0].set_title('q vs εa')

# 2. Trajetória de tensões (p' vs q)
axes[0, 1].plot(p_prime, q, 'r-', linewidth=2, label='Trajetória CU')
axes[0, 1].set_xlabel("p' (kPa)")
axes[0, 1].set_ylabel('q (kPa)')
axes[0, 1].grid(True, alpha=0.3)
axes[0, 1].set_title("Trajetória de Tensões (p' vs q)")
axes[0, 1].legend()

# 3. NOVO: Poropressão vs εa
axes[0, 2].plot(eps_a, u, 'g-', linewidth=2)
axes[0, 2].set_xlabel('Deformação Axial εa (%)')
axes[0, 2].set_ylabel('Poropressão u (kPa)')
axes[0, 2].grid(True, alpha=0.3)
axes[0, 2].set_title('Poropressão vs εa')

# 4. NOVO: σ1 e σ3 (totais)
axes[1, 0].plot(eps_a, sigma1_total, 'b-', linewidth=2, label='σ1 (total)')
axes[1, 0].plot(eps_a, sigma3_total_array, 'r--', linewidth=2, label='σ3 (total)')
axes[1, 0].set_xlabel('Deformação Axial εa (%)')
axes[1, 0].set_ylabel('Tensões Principais (kPa)')
axes[1, 0].grid(True, alpha=0.3)
axes[1, 0].set_title('Tensões Totais vs εa')
axes[1, 0].legend()

# 5. NOVO: q/σ3 normalizado (estilo laboratório)
axes[1, 1].plot(eps_a, q_sigma3, 'm-', linewidth=2)
axes[1, 1].set_xlabel('Deformação Axial εa (%)')
axes[1, 1].set_ylabel('q/σ3')
axes[1, 1].grid(True, alpha=0.3)
axes[1, 1].set_title('Relação q/σ3 (normalizado)')

# 6. NOVO: u/σ3 normalizado (parâmetro de Skempton)
axes[1, 2].plot(eps_a, u_sigma3, 'c-', linewidth=2)
axes[1, 2].set_xlabel('Deformação Axial εa (%)')
axes[1, 2].set_ylabel('u/σ3')
axes[1, 2].grid(True, alpha=0.3)
axes[1, 2].set_title('Relação u/σ3 (normalizado)')

plt.tight_layout()
plt.savefig('ensaio_CU_ajustado.png', dpi=300, bbox_inches='tight')
print("✓ Gráfico salvo: ensaio_CU_ajustado.png")
plt.show()

# Verificar conservação de volume (CU)
eps_v = np.array(results['volumetric_strain'])
print(f"\nVerificação de volume constante (CU):")
print(f"  εv máximo: {np.max(np.abs(eps_v)):.2e}")
print(f"  εv final: {eps_v[-1]:.2e}")
if np.max(np.abs(eps_v)) < 1e-10:
    print("  ✓ Volume constante verificado!")
else:
    print(f"  ⚠ Pequena variação volumétrica detectada (esperado para incrementos finitos)")

# Resumo final
print(f"\nResumo do ensaio:")
print(f"  σ3_total (constante): {sigma3_total} kPa")
print(f"  q_max: {np.max(q):.2f} kPa")
print(f"  u_max: {np.max(u):.2f} kPa")
print(f"  (q/σ3)_max: {np.max(q_sigma3):.3f}")
print(f"  (u/σ3)_max: {np.max(u_sigma3):.3f}")
