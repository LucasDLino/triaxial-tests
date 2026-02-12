# Ensaios Triaxiais - Modelo de Mohr-Coulomb

Simulacao numerica de ensaios triaxiais (CD, CU, UU) usando o modelo constitutivo de Mohr-Coulomb com return mapping completo e validacao analitica.

---

## Inicio Rapido

```bash
# Instale dependencias
pip install -r requirements.txt

# Execute validacao completa
python main.py
```

**Saida esperada:**
- Ensaios CD, CU, UU executados e validados
- Hardening e softening verificados
- Circulos de Mohr comparativos
- Graficos salvos: `validacao_completa.png`, `circulos_mohr_comparativo.png`

---

## Estrutura do Projeto

```
triaxial/
├── main.py                 # Ponto de entrada (executa validacao.py)
├── validacao.py            # Validacao completa (CD, CU, UU, hardening, softening)
├── triaxial.py             # Simulador de ensaios triaxiais
├── mohr_courlomb.py        # Modelo constitutivo Mohr-Coulomb
├── utils.py                # Funcoes auxiliares
├── ALGORITMO_TRIAXIAL.md   # Documentacao tecnica detalhada
├── README.md               # Este arquivo
└── requirements.txt        # Dependencias (numpy, matplotlib)
```

---

## Arquitetura

### Separacao de Responsabilidades

| Componente | Responsabilidade |
|------------|------------------|
| **MohrCoulombModel** | Modelo constitutivo: recebe deformacao (e), retorna tensao efetiva (s') |
| **TriaxialTest** | Simulador: impoe condicoes de contorno, calcula poropressao |
| **validacao.py** | Validacao: compara resultados numericos com solucao analitica |

### Por Que o Modelo Trabalha com Tensoes Efetivas?

O principio de Terzaghi: comportamento mecanico do solo e governado por tensoes efetivas.
- Criterio de Mohr-Coulomb: tau_f = c' + sigma' * tan(phi')
- A poropressao (u) e calculada externamente pelo simulador de ensaio

Ver [ALGORITMO_TRIAXIAL.md](ALGORITMO_TRIAXIAL.md) para explicacao completa.

---

## Tipos de Ensaio

| Ensaio | Descricao | Condicao | Drenagem |
|--------|-----------|----------|----------|
| **CD** | Consolidated Drained | sigma3' = constante | u = 0 |
| **CU** | Consolidated Undrained | eps_v = 0 | u calculado |
| **UU** | Unconsolidated Undrained | eps_v = 0 | u calculado |

---

## Exemplo de Uso

### Ensaio CD com Validacao

```python
from mohr_courlomb import MohrCoulombModel
from triaxial import TriaxialTest
import numpy as np

# Criar modelo
model = MohrCoulombModel(
    E=25000,        # kPa
    nu=0.3,
    phi_deg=25,     # graus
    cohesion=20,    # kPa
    psi_deg=0       # dilatancia (0 = nao-associado)
)

# Executar ensaio CD
test = TriaxialTest(model, sigma3=50, test_type="CD")
r = test.run(eps_max=0.10, steps=100)

# Acessar resultados
q_max = max(r['q'])
sigma1_max = max(r['sigma1'])
sigma3 = r['sigma3'][0]

# Validar contra solucao analitica
phi_rad = np.radians(25)
Nf = (1 + np.sin(phi_rad)) / (1 - np.sin(phi_rad))
q_teo = sigma3 * (Nf - 1) + 2 * 20 * np.sqrt(Nf)

print(f"q numerico: {q_max:.1f} kPa")
print(f"q teorico:  {q_teo:.1f} kPa")
print(f"Erro: {100*abs(q_max-q_teo)/q_teo:.1f}%")
```

### Ensaio CU com Poropressao

```python
# Executar ensaio CU
test = TriaxialTest(model, sigma3=100, test_type="CU")
r = test.run(eps_max=0.10, steps=100)

# Acessar tensoes efetivas e poropressao
sigma1_eff = r['sigma1']        # Tensao axial EFETIVA
sigma3_eff = r['sigma3']        # Tensao radial EFETIVA
u = r['pore_pressure']          # Poropressao

# Tensoes totais (se necessario)
sigma1_total = [s1 + ui for s1, ui in zip(sigma1_eff, u)]
sigma3_total = 100  # Constante (pressao de camara)
```

### Hardening (Endurecimento)

```python
# Coesao evolui com deformacao plastica
sampling_pairs = [
    [0.00, 20],   # eps_p = 0    -> c = 20 kPa
    [0.05, 50],   # eps_p = 5%   -> c = 50 kPa
]

model = MohrCoulombModel(
    E=25000, nu=0.3, phi_deg=25, cohesion=20, psi_deg=0,
    sampling_pairs=sampling_pairs
)
```

---

## Resultados da Validacao

Executando `python main.py`:

```
============================================================
VALIDACAO COMPLETA DO MODELO MOHR-COULOMB
============================================================

[OK] CD: Erro < 5% vs solucao analitica (erro = 0.0%)
[OK] CU: Geracao de poropressao verificada
[OK] UU: cu independente de sigma3
[OK] Hardening: Aumento progressivo de q
[OK] Softening: Reducao de resistencia pos-pico
```

### Graficos Gerados

1. **validacao_completa.png** - Curvas q vs eps_a, poropressao, hardening/softening
2. **circulos_mohr_comparativo.png** - Circulos de Mohr para CD, CU, UU com envoltorias

---

## API Reference

### TriaxialTest

```python
TriaxialTest(model, sigma3, test_type="CU")
```

**Parametros:**
- `model`: Instancia de MohrCoulombModel
- `sigma3`: Tensao de confinamento (kPa)
  - CD: tensao efetiva sigma3'
  - CU/UU: tensao total (pressao de camara)
- `test_type`: "CD", "CU" ou "UU"

**Metodo run():**

```python
run(eps_max=0.15, steps=200) -> dict
```

**Retorna dicionario com:**
- `axial_strain`: Deformacao axial
- `q`: Tensao desviadora (efetiva)
- `p`: Tensao media (efetiva)
- `volumetric_strain`: Deformacao volumetrica
- `sigma1`: Tensao principal maior (EFETIVA)
- `sigma3`: Tensao principal menor (EFETIVA)
- `pore_pressure`: Poropressao u

### MohrCoulombModel

```python
MohrCoulombModel(E, nu, phi_deg, cohesion, psi_deg=0, sampling_pairs=None)
```

**Parametros:**
- `E`: Modulo de Young (kPa)
- `nu`: Coeficiente de Poisson
- `phi_deg`: Angulo de atrito (graus)
- `cohesion`: Coesao (kPa)
- `psi_deg`: Angulo de dilatancia (graus)
- `sampling_pairs`: Lista [[eps_p, c], ...] para hardening/softening

---

## Documentacao Tecnica

Para detalhes sobre:
- Algoritmo de integracao elastoplastica
- Return mapping para psi = 0
- Calculo de poropressao em ensaios nao-drenados
- Distincao tensoes efetivas vs totais

Ver **[ALGORITMO_TRIAXIAL.md](ALGORITMO_TRIAXIAL.md)**

---

## Dependencias

- Python 3.7+
- numpy
- matplotlib

```bash
pip install -r requirements.txt
```

---

## Licenca

MIT License
