# Algoritmo de Simulação de Ensaios Triaxiais

Este documento explica o algoritmo implementado em `triaxial.py` para simulação de ensaios triaxiais CD, CU e UU usando o modelo constitutivo de Mohr-Coulomb.

---

## 1. Conceito Fundamental: Tensões Efetivas vs Totais

### 1.1 Definições

| Símbolo | Nome | Descrição |
|---------|------|-----------|
| **σ** | Tensão Total | Tensão medida externamente (água + esqueleto) |
| **σ'** | Tensão Efetiva | Tensão transmitida pelo esqueleto sólido |
| **u** | Poropressão | Pressão da água nos poros |

**Relação fundamental (Terzaghi)**:
$$\sigma = \sigma' + u$$

Portanto:
$$\sigma' = \sigma - u$$

### 1.2 Regra de Ouro

> **O modelo de Mohr-Coulomb SEMPRE trabalha com TENSÕES EFETIVAS (σ').**

O critério de ruptura de Mohr-Coulomb é formulado em termos de tensões efetivas:
$$\tau_f = c' + \sigma'_n \cdot \tan(\phi')$$

Isso significa que:
- A **entrada** do modelo (deformações) produz **tensões efetivas** na saída
- A conversão para tensões totais é feita **externamente** quando necessário (em ensaios CU/UU)

---

## 2. Fluxo de Dados no Sistema

```
┌─────────────────────────────────────────────────────────────────────┐
│                        TRIAXIAL TEST                                 │
│  (classe TriaxialTest em triaxial.py)                               │
├─────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  ENTRADA DO USUÁRIO:                                                │
│  ┌─────────────────────────────────────────────────────────┐        │
│  │ sigma3 = 50 kPa  (tensão de confinamento)              │        │
│  │ test_type = "CU" ou "CD" ou "UU"                        │        │
│  └─────────────────────────────────────────────────────────┘        │
│                              │                                       │
│                              ▼                                       │
│  INTERPRETAÇÃO POR TIPO DE ENSAIO:                                  │
│  ┌─────────────────────────────────────────────────────────┐        │
│  │ CD: sigma3 = σ₃' (efetiva, pois u=0)                   │        │
│  │ CU: sigma3 = σ₃_total (câmara), σ₃'₀ = σ₃ após consol. │        │
│  │ UU: sigma3 = σ₃_total (câmara), σ₃' = fixo (1 kPa)     │        │
│  └─────────────────────────────────────────────────────────┘        │
│                              │                                       │
│                              ▼                                       │
│  SIMULAÇÃO DO CARREGAMENTO:                                         │
│  ┌─────────────────────────────────────────────────────────┐        │
│  │ Para cada passo:                                        │        │
│  │   1. Calcular incremento de deformação (ε)             │        │
│  │   2. Passar ε para o modelo Mohr-Coulomb               │        │
│  │                        │                                │        │
│  │                        ▼                                │        │
│  │   ┌────────────────────────────────────────────┐       │        │
│  │   │  MODELO MOHR-COULOMB (mohr_courlomb.py)    │       │        │
│  │   │  • Recebe: tensor de deformações ε         │       │        │
│  │   │  • Calcula: integração elastoplástica      │       │        │
│  │   │  • Retorna: tensor de TENSÕES EFETIVAS σ'  │       │        │
│  │   └────────────────────────────────────────────┘       │        │
│  │                        │                                │        │
│  │                        ▼                                │        │
│  │   3. Extrair σ₁' e σ₃' do tensor retornado             │        │
│  │   4. Calcular q = σ₁' - σ₃' (sempre efetivo)           │        │
│  │   5. Para CU/UU: calcular u e tensões totais           │        │
│  └─────────────────────────────────────────────────────────┘        │
│                                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 2.1 Arquitetura: Modelo Constitutivo vs Simulador de Ensaio

### Por Que o Modelo Trabalha Apenas com Tensões Efetivas?

O modelo constitutivo de Mohr-Coulomb implementa a **relação tensão-deformação do esqueleto sólido**. Esta é a base da mecânica dos solos saturados:

- O **princípio de Terzaghi** estabelece que o comportamento mecânico do solo (deformação, ruptura) é governado pelas **tensões efetivas**
- A poropressão (u) não causa deformação do esqueleto — apenas transmite pressão hidrostática
- O critério de ruptura de Mohr-Coulomb é formulado em tensões efetivas: τf = c' + σ'·tan(φ')

Portanto, a separação é:

| Componente | Responsabilidade |
|------------|------------------|
| **MohrCoulombModel** | Recebe ε, retorna σ' (relação constitutiva em tensões efetivas) |
| **TriaxialTest** | Impõe condições de contorno (σ₃' const, εᵥ=0), calcula u |

### Interface do Modelo Constitutivo

O modelo expõe um único método principal:

```python
stress = model.update(strain_total)
#   ▲                      ▲
#   │                      │
#   └── σ' (efetivo)       └── ε (tensor 3×3)
```

**Entradas e Saídas**:
- **Entrada**: Tensor de deformação total ε (3×3)
- **Saída**: Tensor de tensão efetiva σ' (3×3)

### Acoplamento com o Simulador de Ensaio

⚠️ **Limitação atual**: O `triaxial.py` acessa atributos internos do modelo:

```python
model.K           # Módulo volumétrico
model.E, model.nu # Propriedades elásticas  
model.strain      # Estado interno de deformação
model.is_plastic  # Flag de plastificação
model.equivalent_plastic_strain  # Para backup/restore
```

Isso significa que para **substituir** o MohrCoulombModel por outro modelo (ex: Cam-Clay, Drucker-Prager), o novo modelo precisaria:

1. Implementar `update(strain_total) → stress`
2. Expor os mesmos atributos: `K`, `E`, `nu`, `strain`, `is_plastic`
3. Trabalhar com a mesma convenção (compressão positiva)

**Para tornar o sistema mais modular**, seria necessário definir uma interface abstrata (ABC) com apenas os métodos essenciais.

---

## 3. Cálculo da Poropressão (u) em Ensaios Não-Drenados

### 3.1 De Onde Vem σ₃_total?

No ensaio CU ou UU, o usuário especifica a **pressão de câmara** (confinamento aplicado externamente). Esta é uma **tensão total** porque inclui a pressão da água.

```python
# Ao criar o ensaio CU:
test = TriaxialTest(model, sigma3=100, test_type='CU')
#                          ▲
#                          │
#                          └── Esta é σ₃_total = 100 kPa
#                              (pressão de câmara = água + esqueleto)
```

### 3.2 Interpretação por Tipo de Ensaio

| Ensaio | sigma3 fornecido | Interpretação |
|--------|------------------|---------------|
| **CD** | σ₃' = σ₃_total | São iguais (u = 0, drenado) |
| **CU** | σ₃_total (câmara) | Após consolidação: σ₃' inicial = σ₃_total |
| **UU** | σ₃_total (câmara) | σ₃' inicial = fixo (~1 kPa), u₀ = σ₃_total - σ₃' |

### 3.3 Fórmula da Poropressão

**Durante o cisalhamento não-drenado (CU/UU)**:

A pressão de câmara σ₃_total permanece **constante** (controlada externamente). O modelo Mohr-Coulomb retorna tensões efetivas σ'. A poropressão é calculada por:

$$u = \sigma_{3,total} - \sigma_3'$$

Onde:
- **σ₃_total**: Valor fornecido pelo usuário (pressão de câmara, constante)
- **σ₃'**: Calculado pelo modelo Mohr-Coulomb (varia durante o ensaio)
- **u**: Poropressão resultante (cresce durante carregamento não-drenado)

### 3.4 Exemplo Numérico (CU)

```
Entrada: sigma3 = 100 kPa (pressão de câmara)

Após consolidação (início do cisalhamento):
  σ₃_total = 100 kPa (constante)
  σ₃'      = 100 kPa (consolidado)
  u        = 100 - 100 = 0 kPa

Durante cisalhamento (εᵥ = 0 forçado):
  σ₃_total = 100 kPa (ainda constante!)
  σ₃'      = 70 kPa  (calculado pelo modelo, diminui)
  u        = 100 - 70 = 30 kPa (aumenta)
  
Na ruptura:
  σ₃_total = 100 kPa
  σ₃'      = 50 kPa
  u        = 100 - 50 = 50 kPa
```

---

## 4. Tipos de Ensaio: Detalhamento

### 4.1 Ensaio CD (Consolidated Drained)

**Características físicas**:
- Drenagem aberta durante cisalhamento
- Poropressão u = 0 sempre
- σ' = σ (tensões efetivas = totais)

**Restrição**:
- σ₃' deve permanecer constante

**Algoritmo**:
```
Para cada passo:
  1. Estimar ε_r (deformação radial)
  2. Aplicar incremento: d_ε = [d_ε_a, ε_r, ε_r]
  3. Calcular σ' via Mohr-Coulomb
  4. Se σ₃' ≠ σ₃_target:
     - Ajustar ε_r iterativamente
     - Repetir até σ₃' ≈ σ₃_target
```

**Saída**: σ₁', σ₃', q, p' (todas efetivas = totais)

### 4.2 Ensaio CU (Consolidated Undrained)

**Características físicas**:
- Consolidado até σ₃_total antes do cisalhamento
- Drenagem fechada durante cisalhamento
- Volume constante (εᵥ = 0)
- Poropressão se desenvolve (u ≠ 0)

**Restrição**:
- εᵥ = ε_a + 2·ε_r = 0 → ε_r = -ε_a/2

**Algoritmo**:
```
Para cada passo:
  1. Calcular ε_r = -d_ε_a / 2  (volume constante)
  2. Aplicar incremento: d_ε = [d_ε_a, ε_r, ε_r]
  3. Calcular σ' via Mohr-Coulomb
  4. Calcular u = σ₃_total - σ₃'
  5. Tensões totais: σ_total = σ' + u
```

**Saída**: 
- σ₁', σ₃', q, p' (efetivas - do modelo)
- u (calculada)
- σ₁_total, σ₃_total (totais - calculadas)

### 4.3 Ensaio UU (Unconsolidated Undrained)

**Características físicas**:
- SEM consolidação prévia
- Tensões efetivas iniciais = estado da amostra (fixo)
- Toda pressão de câmara vira poropressão inicial

**Diferença chave do CU**:
- No CU: σ₃' inicial = σ₃_total (consolida até igualar)
- No UU: σ₃' inicial = fixo (~1 kPa), independente de σ₃_total

**Resultado**:
- q_max é o mesmo para qualquer σ₃_total (cu constante)
- Envoltória total: φᵤ = 0 (horizontal)

---

## 5. Por Que o Controle Iterativo é Necessário no CD?

### 5.1 O Problema

O modelo Mohr-Coulomb é **strain-driven** (recebe deformações, retorna tensões).

Para manter σ₃' = constante, não sabemos ε_r a priori:
- No regime **elástico**: ε_r ≈ -ν·ε_a funciona bem
- No regime **plástico**: depende de ψ (dilatância), não de ν

### 5.2 A Solução

Usar iteração para encontrar ε_r que satisfaça σ₃' = σ₃_target:

```
1. Estimar ε_r inicial (aproximação elástica)
2. Aplicar deformação e calcular σ'
3. Verificar erro = σ₃' - σ₃_target
4. Se erro > tolerância:
   - Corrigir ε_r proporcionalmente ao erro
   - Voltar ao passo 2
```

### 5.3 Por Que CU/UU Não Precisam de Iteração?

Porque a restrição deles (εᵥ = 0) é sobre **deformações**, não tensões.

Podemos impor diretamente: ε_r = -ε_a/2

Não precisamos verificar tensão nenhuma — o modelo calcula as tensões resultantes.

### 5.4 O Mesmo Modelo Constitutivo é Usado em Todos os Ensaios

Independente do tipo de ensaio (CD, CU, UU), a chamada ao modelo constitutivo é **sempre a mesma**:

```python
stress = model.update(total_strain)
```

O que **muda entre os ensaios** é:
- **Como** o tensor de deformação `total_strain` é construído
- **O que** se faz com o tensor de tensão retornado

| Ensaio | Construção de ε | Pós-processamento de σ' |
|--------|-----------------|-------------------------|
| **CD** | ε_r iterado até σ₃' = constante | u = 0, σ = σ' |
| **CU** | ε_r = -ε_a/2 (volume constante) | u = σ₃_total - σ₃' |
| **UU** | ε_r = -ε_a/2 (volume constante) | u = σ₃_total - σ₃' |

Isso significa que o **modelo constitutivo é completamente substituível**. Se você implementar outro modelo (Cam-Clay, Drucker-Prager, hipoplástico) com a mesma interface `update(ε) → σ'`, ele funcionará automaticamente em todos os tipos de ensaio.

---

## 6. Correção do Return Mapping para ψ = 0

### 6.1 O Problema Original

A formulação clássica de return mapping assume ψ ≠ 0. Com ψ = 0, ocorrem inversões numéricas inválidas.

### 6.2 A Solução

Para ψ = 0, projeta-se σ₁ diretamente para a superfície de ruptura mantendo σ₃ fixo:

$$\sigma_1^{limite} = K_p \cdot \sigma_3 + 2c\sqrt{K_p}$$

$$K_p = \frac{1 + \sin(\phi')}{1 - \sin(\phi')}$$

---

## 7. Círculos de Mohr: Por Que a Envoltória Total é Calculada por Regressão?

### 7.1 Envoltória Efetiva (usando parâmetros de entrada)

A envoltória de tensões **efetivas** usa diretamente φ' e c' fornecidos:

$$\tau = c' + \sigma' \cdot \tan(\phi')$$

### 7.2 Envoltória Total (requer regressão)

A envoltória de tensões **totais** não pode usar os mesmos parâmetros porque:
- A poropressão u **varia** com σ₃ de forma não-linear
- A inclinação da envoltória total **difere** da efetiva

Por isso, calculamos (φ_cu, c_cu) por regressão linear sobre os círculos totais:
- Estes são os parâmetros **aparentes** que seriam medidos se apenas tensões totais fossem conhecidas

---

## 8. O Papel da Dilatância (ψ) nos Diferentes Ensaios

### 8.1 O Que é Dilatância?

A dilatância (ψ) é um **parâmetro do material** que descreve a tendência do solo de mudar de volume durante cisalhamento plástico:

| Valor de ψ | Comportamento | Solo típico |
|------------|---------------|-------------|
| ψ > 0 | Dilata (expande) | Areia densa, argila OC |
| ψ = 0 | Volume constante | Argila NC, areia fofa |
| ψ < 0 | Contrai | Solo muito fofo |

### 8.2 Dilatância em Cada Tipo de Ensaio

| Ensaio | Volume livre? | Efeito de ψ |
|--------|---------------|-------------|
| **CD** | ✅ Sim | ψ controla εᵥ plástico diretamente |
| **CU** | ❌ Não (εᵥ = 0) | ψ afeta redistribuição de tensões → poropressão |
| **UU** | ❌ Não (εᵥ = 0) | Mesmo que CU |

### 8.3 Por Que Usamos ψ = 0 na Validação?

1. **Simplicidade numérica**: ψ = 0 evita problemas no return mapping clássico
2. **Representativo para argilas**: Argilas saturadas normalmente consolidadas têm ψ ≈ 0
3. **Conservador**: ψ = 0 é a hipótese mais conservadora (sem dilatância)

### 8.4 Quando Usar ψ ≠ 0?

Para **areias densas** ou **argilas pré-adensadas**, seria necessário:
- ψ = φ/3 (aproximação comum para areias)
- ψ = φ (plasticidade associada - limite superior)

**No ensaio CD com ψ > 0**:
- O solo dilata durante cisalhamento
- εᵥ aumenta (expansão volumétrica)
- O controle iterativo de σ₃ ajusta ε_r automaticamente

**No ensaio CU com ψ > 0**:
- O solo "quer" dilatar mas não pode (εᵥ = 0 forçado)
- Isso gera **poropressão negativa** (sucção)
- O modelo calcula σ₃' maior → u = σ₃_total - σ₃' pode ser negativo

### 8.5 Resumo

A dilatância é uma propriedade do **material**, não do ensaio. A escolha de ψ depende do solo que está sendo modelado:

```
Argila NC saturada:     ψ ≈ 0
Argila OC / Areia densa: ψ = φ/3 a φ/2
Validação numérica:      ψ = 0 (simplifica implementação)
```

---

## 9. Fórmulas de Referência

### Critério de Mohr-Coulomb (tensões efetivas)
$$f = \sigma_1' - K_p \cdot \sigma_3' - 2c'\sqrt{K_p} = 0$$

### Tensão Desviadora de Ruptura (CD)
$$q_f = \sigma_3' \cdot (K_p - 1) + 2c'\sqrt{K_p}$$

### Poropressão (CU)
$$u = \sigma_{3,total} - \sigma_3'$$

### Resistência Não-Drenada
$$c_u = \frac{q_{max}}{2}$$
