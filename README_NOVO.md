# 🏗️ Ensaios Triaxiais - Modelo de Mohr-Coulomb

Simulação numérica de ensaios triaxiais (CD, CU, UU) usando o modelo constitutivo de Mohr-Coulomb com return mapping completo e **validação analítica**.

[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org/)
[![NumPy](https://img.shields.io/badge/NumPy-required-green.svg)](https://numpy.org/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-required-orange.svg)](https://matplotlib.org/)

---

## 🚀 **Início Rápido**

```bash
# Clone o repositório
cd triaxial

# Instale dependências
pip install -r requirements.txt

# Execute exemplos básicos com validação analítica 
python main.py
```

**Saída esperada:**
- 3 ensaios triaxiais (CD, CU, UU) executados
- Validação analítica automática
- Gráficos salvos em PNG
- Aprovação/reprovação com % de erro

---

## 📦 **Instalação**

### Requisitos
- Python 3.7+
- numpy
- matplotlib

```bash
pip install -r requirements.txt
```

---

## 📁 **Estrutura do Projeto**

```
triaxial/
├── 🎯 main.py                      # ⭐ Exemplos básicos (COMECE AQUI!)
├── 🔬 validacao_analitica.py       # Funções de validação analítica
├── 🎨 exemplos_completos.py        # Análises avançadas
├── 🔧 triaxial.py                  # Classe TriaxialTest (♻️ ATUALIZADO)
├── 📐 mohr_courlomb.py             # Modelo constitutivo
├── 🛠️ utils.py                     # Funções auxiliares
├── 📖 ESTRUTURA_ARQUIVOS.md        # Documentação detalhada
├── 📝 CHANGELOG_TRIAXIAL.md        # Histórico de alterações
└── 📋 requirements.txt             # Dependências
```

### 🎯 Fluxo Recomendado

| Nível | Arquivo | Descrição |
|-------|---------|-----------|
| **Iniciante** 🌱 | `main.py` | 3 exemplos básicos com validação |
| **Intermediário** 🌿 | `exemplos_completos.py` | Hardening, dilatância, comparações |
| **Avançado** 🌳 | `validacao_analitica.py` | Criar suas próprias validações |

---

## 🎯 **Principais Recursos**

### ✅ Tipos de Ensaio
- **CD** (Consolidated Drained) - Drenado, volume varia
- **CU** (Consolidated Undrained) - Não drenado, volume constante, **com poropressão**
- **UU** (Unconsolidated Undrained) - Não consolidado, não drenado

### ✅ Modelo Constitutivo de Mohr-Coulomb
- ✓ Return mapping completo (1, 2 e multi-vector)
- ✓ Hardening/softening piecewise linear
- ✓ Dilatância não-associada configurável
- ✓ Integração elastoplástica robusta

### ✅ **NOVO: Validação Analítica** 🔬
- ✓ Comparação automática com solução analítica
- ✓ Verificação de conservação de volume (CU/UU)
- ✓ Cálculo de erro relativo (%)
- ✓ Aprovação/reprovação com tolerância configurável

### ✅ **NOVO: Rastreamento Completo (CU)** 📊
- ✓ **Poropressão** u = σ3_total - σ3'
- ✓ **Tensões totais** (σ1, σ3)
- ✓ **Relações normalizadas** (q/σ3, u/σ3)
- ✓ **100% compatível** com código antigo

---

## 📖 **Exemplos de Uso**

### 🔹 Exemplo 1: Ensaio CD com Validação Automática

```python
from mohr_courlomb import MohrCoulombModel
from triaxial import TriaxialTest
from validacao_analitica import validar_ensaio_CD

# Criar modelo
model = MohrCoulombModel(
    E=50000,        # kPa
    nu=0.3,
    phi_deg=30,     # graus
    cohesion=20,    # kPa
    psi_deg=0       # sem dilatância
)

# Executar ensaio CD
test = TriaxialTest(model, sigma3=100, test_type="CD")
results = test.run(eps_max=0.15, steps=200)

# Validar automaticamente
validacao = validar_ensaio_CD(
    sigma3=100,
    phi_deg=30,
    cohesion=20,
    q_numerico=results['q'],
    p_numerico=results['p']
)

# Resultado: 
# ✓ VALIDAÇÃO APROVADA
# Erro em q: 0.02%
```

**Saída:**
```
VALIDAÇÃO ANALÍTICA - ENSAIO CD
======================================================================
q (tensão desviadora no pico):
  Analítico:  173.21 kPa
  Numérico:   173.18 kPa
  Erro:       0.02% ✓

✓ VALIDAÇÃO APROVADA
```

---

### 🔹 Exemplo 2: Ensaio CU com Poropressão

```python
from mohr_courlomb import MohrCoulombModel
from triaxial import TriaxialTest
from validacao_analitica import validar_ensaio_CU

# Criar modelo
model = MohrCoulombModel(E=50000, nu=0.3, phi_deg=30, cohesion=20)

# Executar ensaio CU
test = TriaxialTest(model, sigma3=100, test_type="CU")
results = test.run()

# NOVO: Acessar poropressão e tensões totais
eps_a = results['axial_strain']
q = results['q']                    # Desviadora efetiva
u = results['pore_pressure']        # Poropressão (NOVO!)
sigma1_total = results['sigma1']    # Tensão axial total (NOVO!)
sigma3_total = results['sigma3']    # Tensão radial total (NOVO!)

# Validar conservação de volume
validacao = validar_ensaio_CU(results['volumetric_strain'])

# Resultado: 
# ✓ VALIDAÇÃO APROVADA - Volume conservado
# |εv|_max: 1.23e-12
```

---

### 🔹 Exemplo 3: Compatibilidade Retroativa

```python
# Código antigo ainda funciona! ✓
test = TriaxialTest(model, sigma3=100, test_type="CD")

# Desempacotamento como tupla (forma antiga)
eps, q, p, epsv = test.run()  # ✓ Funciona!

# OU acesso por dicionário (forma nova)
results = test.run()
u = results['pore_pressure']   # ✓ Novos dados disponíveis!
```

---

### 🔹 Exemplo 4: Hardening Piecewise Linear

```python
# Definir curva de hardening/softening
sampling_pairs = [
    [0.0, 10],    # εp = 0    → c = 10 kPa
    [0.02, 25],   # εp = 0.02 → c = 25 kPa (hardening)
    [0.05, 30],   # εp = 0.05 → c = 30 kPa (pico)
    [0.10, 25],   # εp = 0.10 → c = 25 kPa (softening)
]

model = MohrCoulombModel(
    E=50000, nu=0.3, phi_deg=30, cohesion=10, psi_deg=10,
    sampling_pairs=sampling_pairs  # Curva personalizada
)

test = TriaxialTest(model, sigma3=100, test_type="CD")
results = test.run(eps_max=0.20, steps=300)
```

---

## 📊 **Visualizações**

### Gráficos Gerados Automaticamente

Ao executar `main.py`:

1. **ensaio_CD.png** - Resultados do ensaio drenado
   - q vs εa (com pico marcado)
   - Trajetória de tensões (p' vs q)

2. **ensaio_CU.png** - Resultados do ensaio não drenado
   - q vs εa
   - Trajetória de tensões

3. **poropressao_CU.png** - Evolução da poropressão
   - u vs εa (importante para análise CU!)

4. **ensaio_UU.png** - Resultados não consolidado

### Exemplos Avançados (`exemplos_completos.py`)

- **comparacao_tipos_ensaio.png** - CD vs CU vs UU
- **analise_hardening.png** - Evolução de coesão + trajetórias
- **efeito_dilatancia.png** - Influência do ângulo ψ

---

## 🔬 **Validação Analítica**

### Solução Analítica de Mohr-Coulomb (Ensaio CD)

Para um ensaio drenado, a tensão principal maior no pico é:

$$
\sigma_1' = \sigma_3' \times N_\phi + 2c\sqrt{N_\phi}
$$

onde:

$$
N_\phi = \frac{1 + \sin\phi}{1 - \sin\phi}
$$

**Implementação:**

```python
from validacao_analitica import sigma1_pico_analitico, q_pico_analitico

# Calcular valores analíticos
sigma1_analitico = sigma1_pico_analitico(sigma3=100, phi_rad=0.524, cohesion=20)
q_analitico = q_pico_analitico(sigma3=100, phi_rad=0.524, cohesion=20)

print(f"σ1' analítico: {sigma1_analitico:.2f} kPa")
print(f"q analítico: {q_analitico:.2f} kPa")
```

### Conservação de Volume (Ensaio CU/UU)

Para ensaios não drenados:

$$
\varepsilon_v = \varepsilon_a + 2\varepsilon_r = 0
$$

**Validação:**

```python
from validacao_analitica import validar_ensaio_CU

# Verifica se |εv| < tolerância
validacao = validar_ensaio_CU(
    epsv_numerico=results['volumetric_strain'],
    tolerancia_epsv=1e-6
)

print(f"|εv|_max: {validacao['epsv_max']:.2e}")
# Saída: |εv|_max: 1.23e-12  ✓
```

---

## 🧪 **Suíte de Testes**

Execute todos os testes de validação:

```bash
python main.py
```

**Saída esperada:**

```
█████████████████████████████████████████████████████████████████████
█                                                                   █
█    ENSAIOS TRIAXIAIS - MODELO DE MOHR-COULOMB                    █
█                                                                   █
█████████████████████████████████████████████████████████████████████

EXEMPLO 1: Ensaio CD (Consolidated Drained)
======================================================================
✓ VALIDAÇÃO APROVADA - Erro: 0.02%

EXEMPLO 2: Ensaio CU (Consolidated Undrained)
======================================================================
✓ VALIDAÇÃO APROVADA - Volume conservado

EXEMPLO 3: Ensaio UU (Unconsolidated Undrained)
======================================================================
✓ VALIDAÇÃO APROVADA - Volume conservado

RESUMO DAS VALIDAÇÕES
======================================================================
  CD: ✓ APROVADO
  CU: ✓ APROVADO
  UU: ✓ APROVADO
======================================================================
```

---

## 📚 **Documentação Adicional**

- **[ESTRUTURA_ARQUIVOS.md](ESTRUTURA_ARQUIVOS.md)** - Guia detalhado de cada arquivo
- **[CHANGELOG_TRIAXIAL.md](CHANGELOG_TRIAXIAL.md)** - Histórico completo de alterações
- Comentários inline em todo o código

---

## 🔧 **API Reference**

### Classe `TriaxialTest`

```python
TriaxialTest(model, sigma3, test_type="CD")
```

**Parâmetros:**
- `model`: Instância de `MohrCoulombModel`
- `sigma3`: Tensão de confinamento (kPa)
  - CD/UU: tensão efetiva
  - CU: tensão total constante (pressão de câmara)
- `test_type`: "CD", "CU" ou "UU"

**Métodos:**

```python
run(eps_max=0.15, steps=200)
```

**Retorna:** `TriaxialResults` (dicionário + compatibilidade com tupla)

```python
{
    'axial_strain': [...],       # εa
    'q': [...],                  # Desviadora q (efetiva)
    'p': [...],                  # Média p' (efetiva)
    'volumetric_strain': [...],  # εv
    'sigma1': [...],             # σ1 total
    'sigma3': [...],             # σ3 total
    'pore_pressure': [...],      # u (poropressão)
    'q_sigma3_ratio': [...],     # q/σ3
    'u_sigma3_ratio': [...]      # u/σ3
}
```

---

### Classe `MohrCoulombModel`

```python
MohrCoulombModel(E, nu, phi_deg, cohesion, psi_deg=0, sampling_pairs=None)
```

**Parâmetros:**
- `E`: Módulo de Young (kPa)
- `nu`: Coeficiente de Poisson (adimensional)
- `phi_deg`: Ângulo de atrito interno (graus)
- `cohesion`: Coesão (kPa)
- `psi_deg`: Ângulo de dilatância (graus, padrão=0)
- `sampling_pairs`: Lista de `[εp, c]` para hardening (opcional)

---

## 🤝 **Contribuindo**

Melhorias bem-vindas! Áreas de interesse:

1. ✨ Modelos constitutivos adicionais (Cam Clay, Drucker-Prager)
2. 🔄 Integração implícita para ensaios CU (Newton-Raphson)
3. 📊 Mais validações analíticas
4. 🧪 Testes unitários automatizados
5. 📈 Interface gráfica interativa

---

## 📜 **Licença**

MIT License - livre para uso acadêmico e comercial

---

## ✨ **Novidades da Versão Atual**

### v2.0 (Fevereiro 2026)

#### 🆕 Adições
- ✅ **Validação analítica completa** para CD, CU, UU
- ✅ **Rastreamento de poropressão** em ensaios CU
- ✅ **Tensões totais vs efetivas** claramente distinguidas
- ✅ **Relações normalizadas** (q/σ3, u/σ3)
- ✅ Arquivo `validacao_analitica.py` com funções de validação
- ✅ `ESTRUTURA_ARQUIVOS.md` - documentação organizacional

#### 🔄 Melhorias
- ✅ `main.py` simplificado - 3 exemplos básicos
- ✅ `exemplos_completos.py` - análises avançadas separadas
- ✅ Compatibilidade 100% retroativa mantida
- ✅ Documentação expandida em todos os arquivos

#### 🐛 Correções
- ✅ Cálculo correto de poropressão em CU
- ✅ Distinção clara σ_total vs σ' no modelo
- ✅ Comentários detalhados sobre limitações

---

## 📞 **Contato**

Para questões técnicas ou sugestões, abra uma issue no repositório.

---

**⭐ Se este projeto foi útil, considere dar uma estrela!**
