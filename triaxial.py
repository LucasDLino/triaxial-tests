import numpy as np
from utils import deviatoric, q_invariant


class TriaxialResults(dict):
    """
    Classe para resultados de ensaio triaxial
    Suporta acesso por dicionário E desempacotamento como tupla (compatibilidade)
    """
    def __iter__(self):
        # COMPATIBILIDADE: Permite desempacotar como tupla
        # eps, q, p, epsv = test.run()
        return iter([self['axial_strain'], self['q'], self['p'], self['volumetric_strain']])


class TriaxialTest:
    """
    Classe para simulação de ensaios triaxiais
    
    Tipos de ensaio:
    - CU: Consolidated Undrained (adensado não drenado)
    - CD: Consolidated Drained (adensado drenado)
    - UU: Unconsolidated Undrained (não adensado não drenado)
    """
    
    def __init__(self, model, sigma3, test_type="CU", OCR=1.0):
        """
        Inicializa o ensaio triaxial
        
        Parameters:
        -----------
        model : MohrCoulombModel
            Modelo constitutivo
        sigma3 : float
            Tensão de confinamento (efetiva para CD/UU, total para CU)
        test_type : str
            Tipo de ensaio ("CU", "CD" ou "UU")
        OCR : float, optional
            Over-Consolidation Ratio (padrão: 1.0 = normalmente consolidado)
            OCR > 1: solo pré-adensado (sobreconsolidado)
            OCR = 1: solo normalmente consolidado
            Usado para definir histórico de tensões (p_c' = OCR × p'0)
        
        CONVENÇÃO DE SINAIS:
        -------------------
        Este código usa a convenção geotécnica padrão:
        - COMPRESSÃO = POSITIVA em p' (tensão média)
        - COMPRESSÃO = NEGATIVA no tensor de tensões (σ < 0)
        - q (desviadora) sempre positiva em compressão triaxial
        """
        self.model = model
        self.test_type = test_type
        self.OCR = OCR
        
        # =====================================================================
        # AJUSTE: Distinguir tensão total vs efetiva por tipo de ensaio
        # =====================================================================
        # CD: σ3 aplicado é tensão EFETIVA (drenado, u=0)
        # CU: σ3 aplicado é tensão TOTAL, consolidação prévia até σ3' = σ3
        # UU: σ3 aplicado é tensão TOTAL, mas NÃO CONSOLIDA
        #     Toda pressão σ3 vai para poropressão (Δu = Δσ3)
        #     Tensões efetivas permanecem constantes (estado inicial da amostra)
        # =====================================================================
        
        if test_type == "CU":
            self.sigma3_total = sigma3
            self.sigma3_effective = sigma3  # Consolidado até σ3
            self._u0 = 0.0  # Poropressão inicial após consolidação
            
        elif test_type == "UU":
            # UU: NÃO CONSOLIDA - tensões efetivas NÃO mudam com σ3
            # Usamos um valor fixo pequeno para σ3' inicial (estado da amostra)
            # σ3_effective_inicial representa o estado de tensões in-situ ou residual
            self.sigma3_total = sigma3
            self._sigma3_eff_inicial = 1.0  # kPa - estado inicial fixo da amostra
            self.sigma3_effective = self._sigma3_eff_inicial  # NÃO depende de σ3!
            self._u0 = sigma3 - self._sigma3_eff_inicial  # Poropressão inicial
            
        else:  # CD
            self.sigma3_total = sigma3
            self.sigma3_effective = sigma3  # Drenado, u = 0
            self._u0 = 0.0
        
        # Listas para armazenar resultados
        self.axial_strain = []
        self.q = []
        self.p = []  # p' efetivo
        self.volumetric_strain = []
        
        # AJUSTE: Adicionar rastreamento de tensões principais e poropressão
        self.sigma1 = []  # Tensão axial total
        self.sigma3 = []  # Tensão radial total
        self.pore_pressure = []  # Poropressão u
        self.q_sigma3_ratio = []  # q/σ3 (normalizado)
        self.u_sigma3_ratio = []  # u/σ3 (normalizado)
        
        # CORREÇÃO: Para ensaios não drenados, rastrear estado de ruptura
        self._undrained_failure_reached = False
        self._q_at_failure = 0.0
        self._stress_at_failure = None
        
        # =====================================================================
        # INICIALIZAÇÃO DO ESTADO DE TENSÕES E DEFORMAÇÕES
        # =====================================================================
        # Convenção: compressão POSITIVA (convenção geotécnica)
        # 
        # NOTA IMPORTANTE: A inicialização de deformação elástica é necessária
        # para manter CONSISTÊNCIA entre tensão e deformação no estado inicial.
        # O modelo constitutivo assume que:
        #   σ = f(ε_total) através de ε_elastic e ε_plastic
        # 
        # Se inicializássemos com σ ≠ 0 mas ε = 0, haveria inconsistência.
        # A deformação elástica inicial corresponde à aplicação da tensão
        # de confinamento σ3' de forma isotrópica.
        # =====================================================================
        
        # Estado isotrópico de confinamento (tensão efetiva)
        # CONVENÇÃO: Compressão POSITIVA (geotécnica)
        self.model.stress = self.sigma3_effective * np.eye(3)
        
        # Deformação elástica correspondente à tensão de confinamento
        # Para estado isotrópico: ε_iso = σ_iso / (3K) (positivo = contração)
        eps_iso = self.sigma3_effective / (3 * self.model.K)
        self.model.elastic_strain = eps_iso * np.eye(3)
        self.model.strain = self.model.elastic_strain.copy()
        
        # Deformação plástica equivalente (ajustado por OCR se necessário)
        # Para OCR > 1: solo teve consolidação prévia maior
        if hasattr(self.model, 'equivalent_plastic_strain'):
            # Para modelos com histórico de consolidação
            # OCR > 1 implica que o solo já foi carregado até p_c' > p'0
            # Aqui apenas documentamos; implementação completa requer
            # modelo com memória de pré-adensamento (ex: Cam Clay)
            self.model.equivalent_plastic_strain = 0.0
        
        # CORREÇÃO: Deformação total acumulada DEVE iniciar igual à deformação do modelo
        # Isso garante que o incremento (strain_total - self.strain) seja correto no 1º passo
        # IMPORTANTE: O modelo calcula internamente strain_increment
        # como (strain_total - self.strain) a cada passo
        self.total_strain = self.model.strain.copy()
    
    def get_radial_strain(self, d_eps_a):
        """
        Calcula a deformação radial baseada no tipo de ensaio
        
        Parameters:
        -----------
        d_eps_a : float
            Incremento de deformação axial
            
        Returns:
        --------
        float
            Deformação radial
        """
        if self.test_type == "CU":
            # Não drenado (consolidado): volume constante (εv = 0)
            # Para deformações infinitesimais: εv = εa + 2εr = 0 → εr = -εa/2
            # Esta aproximação é adequada para incrementos pequenos
            eps_r = -d_eps_a / 2.0
            
        elif self.test_type == "CD":
            # Drenado: deformação lateral baseada em Poisson ELÁSTICO
            # 
            # ⚠️ LIMITAÇÃO IMPORTANTE:
            # Esta é uma APROXIMAÇÃO válida apenas no regime ELÁSTICO.
            # Após plastificação, a deformação lateral real depende do
            # ângulo de dilatância (ψ) do modelo, NÃO de ν.
            # 
            # Para ensaio CD rigoroso, seria necessário:
            # 1. Manter σ3' constante (não εr constante)
            # 2. Calcular εr iterativamente para satisfazer σ3' = const
            # 3. O modelo então determinaria εr baseado em ψ no regime plástico
            # 
            # A aproximação atual (εr = -ν εa) é aceitável para:
            # - Análises qualitativas
            # - Pequenas deformações
            # - Quando comportamento é predominantemente elástico
            eps_r = -self.model.nu * d_eps_a
            
        elif self.test_type == "UU":
            # Não consolidado, não drenado: volume constante (εv = 0)
            # 
            # CORREÇÃO: UU deve usar a MESMA condição de volume constante que CU!
            # A diferença entre CU e UU está no ESTADO INICIAL:
            # - CU: consolidado isotropicamente antes do cisalhamento
            # - UU: aplicação rápida de σ3 sem drenagem/consolidação
            # 
            # Ambos têm caminho de deformação não drenado: εv = 0
            eps_r = -d_eps_a / 2.0  # Volume constante (corrigido)
        else:
            raise ValueError(f"Tipo de ensaio inválido: {self.test_type}")
        
        return eps_r
    
    def run(self, eps_max=0.15, steps=200):
        """
        Executa o ensaio triaxial
        
        Parameters:
        -----------
        eps_max : float
            Deformação axial máxima (incremental a partir do confinamento)
        steps : int
            Número de passos
            
        Returns:
        --------
        dict
            Dicionário com resultados:
            - 'axial_strain': Deformação axial εa
            - 'q': Tensão desviadora q (efetiva)
            - 'p': Tensão média p' (efetiva)
            - 'volumetric_strain': Deformação volumétrica εv
            - 'sigma1': Tensão axial total σ1
            - 'sigma3': Tensão radial total σ3
            - 'pore_pressure': Poropressão u (apenas CU)
            - 'q_sigma3_ratio': q/σ3_total
            - 'u_sigma3_ratio': u/σ3_total (apenas CU)
        """
        d_eps = eps_max / steps
        
        # Armazenar deformação volumétrica inicial (do confinamento)
        eps_v_inicial = np.trace(self.model.strain)
        
        for i in range(steps):
            if self.test_type == "CD":
                # =====================================================================
                # ENSAIO CD: Controle iterativo para manter σ3' = constante
                # =====================================================================
                # Usar MUITOS subpassos com correção proporcional limitada
                # =====================================================================
                sigma3_target = self.sigma3_effective
                
                # Muitos subpassos para estabilidade com hardening/softening
                n_sub = 20
                d_eps_sub = d_eps / n_sub
                
                # Manter valor de eps_r entre subpassos (continuidade)
                if not hasattr(self, '_eps_r_history'):
                    self._eps_r_history = -self.model.nu * d_eps_sub
                eps_r_base = self._eps_r_history
                
                for sub_idx in range(n_sub):
                    # Salvar estado para possível rollback
                    strain_backup = self.total_strain.copy()
                    model_backup = {
                        'strain': self.model.strain.copy(),
                        'elastic_strain': self.model.elastic_strain.copy(),
                        'plastic_strain': self.model.plastic_strain.copy(),
                        'stress': self.model.stress.copy(),
                        'eps_p': self.model.equivalent_plastic_strain
                    }
                    
                    # Aplicar incremento inicial
                    d_strain = np.array([
                        [d_eps_sub, 0, 0],
                        [0, eps_r_base, 0],
                        [0, 0, eps_r_base]
                    ])
                    self.total_strain = strain_backup + d_strain
                    stress = self.model.update(self.total_strain)
                    
                    # Verificar σ3
                    sigma3_atual = np.min(np.linalg.eigvalsh(stress))
                    erro = sigma3_atual - sigma3_target
                    
                    # Correção iterativa (máximo 10 iterações)
                    for iter_count in range(10):
                        if abs(erro) < 0.1:  # Tolerância 0.1 kPa
                            break
                        
                        # Restaurar estado
                        self.model.strain = model_backup['strain'].copy()
                        self.model.elastic_strain = model_backup['elastic_strain'].copy()
                        self.model.plastic_strain = model_backup['plastic_strain'].copy()
                        self.model.stress = model_backup['stress'].copy()
                        self.model.equivalent_plastic_strain = model_backup['eps_p']
                        
                        # Correção de εr: se σ3 > target, precisamos de menos contração (εr mais positivo)
                        # Usar sensibilidade elástica: dσ3/dεr ≈ 2G (para regime elástico)
                        # No regime plástico, sensibilidade é MUITO menor para evitar oscilações
                        G = self.model.E / (2 * (1 + self.model.nu))
                        
                        # Fator de relaxação progressivo para garantir convergência
                        relaxation = 0.5 / (1 + iter_count * 0.2)  # Diminui a cada iteração
                        
                        if self.model.is_plastic:
                            # Plástico: usar sensibilidade muito baixa
                            sensitivity = 2 * G * 0.1 * relaxation
                        else:
                            sensitivity = 2 * G * relaxation
                            
                        delta_eps_r = -erro / sensitivity
                        
                        # Limitar correção para estabilidade - mais conservador
                        max_corr = 0.3 * abs(d_eps_sub)
                        delta_eps_r = np.clip(delta_eps_r, -max_corr, max_corr)
                        eps_r_base = eps_r_base + delta_eps_r
                        
                        # Reaplicar com novo εr
                        d_strain = np.array([
                            [d_eps_sub, 0, 0],
                            [0, eps_r_base, 0],
                            [0, 0, eps_r_base]
                        ])
                        self.total_strain = strain_backup + d_strain
                        stress = self.model.update(self.total_strain)
                        
                        sigma3_atual = np.min(np.linalg.eigvalsh(stress))
                        erro = sigma3_atual - sigma3_target
                
                # Guardar eps_r para próximo passo (continuidade)
                self._eps_r_history = eps_r_base
            else:
                # CU/UU: usar método original (volume constante)
                eps_r = self.get_radial_strain(d_eps)
                
                # Incremento de deformação
                d_strain = np.array([
                    [d_eps, 0, 0],
                    [0, eps_r, 0],
                    [0, 0, eps_r]
                ])
                
                # Atualizar deformação total acumulada
                self.total_strain += d_strain
                
                # =====================================================================
                # CORREÇÃO PARA CU/UU: Detectar e manter estado de ruptura
                # =====================================================================
                # Em ensaios não drenados perfeitamente plásticos, uma vez que o
                # material atinge a superfície de cedência, deveria permanecer nela
                # com q = constante (resistência não drenada cu = q_ruptura/2).
                # =====================================================================
                
                if self._undrained_failure_reached:
                    # Após ruptura: manter tensões constantes (perfeitamente plástico)
                    stress = self._stress_at_failure
                else:
                    # Atualizar modelo normalmente
                    stress = self.model.update(self.total_strain)
                    
                    # Verificar se atingiu ruptura
                    eigenvalues_check = np.linalg.eigvalsh(stress)
                    sigma3_check = np.min(eigenvalues_check)
                    q_current = np.max(eigenvalues_check) - sigma3_check
                    
                    # Detectar ruptura: quando σ3' fica muito baixo ou o modelo plastifica
                    # e q está próximo do máximo
                    if sigma3_check < 1.0 or (self.model.is_plastic and len(self.q) > 0 and q_current < self.q[-1] * 0.9):
                        # Encontrar o estado de pico (q máximo até agora)
                        if len(self.q) > 0:
                            q_max_idx = np.argmax(self.q)
                            self._q_at_failure = self.q[q_max_idx]
                            # Reconstruir tensão de ruptura
                            # Em CU/UU, a ruptura ocorre quando q = q_max
                            # σ1' = q + σ3', σ3' ajustado para q_max
                            sigma1_failure = self.sigma1[q_max_idx] - self.pore_pressure[q_max_idx]
                            sigma3_failure = sigma1_failure - self._q_at_failure
                            self._stress_at_failure = np.diag([sigma1_failure, sigma3_failure, sigma3_failure])
                            self._undrained_failure_reached = True
                            stress = self._stress_at_failure
            
            # Calcular invariantes (tensões EFETIVAS)
            # Convenção: compressão POSITIVA
            eigenvalues = np.linalg.eigvalsh(stress)
            sigma1_prime = np.max(eigenvalues)  # Maior tensão principal (axial)
            sigma3_prime = np.min(eigenvalues)  # Menor tensão principal (radial)
            
            q = sigma1_prime - sigma3_prime
            p_prime = (sigma1_prime + 2 * sigma3_prime) / 3.0
            
            # AJUSTE: Calcular poropressão e tensões totais
            if self.test_type == "CU" or self.test_type == "UU":
                # Em CU/UU: σ3_total = constante (pressão de câmara)
                # Poropressão: u = σ3_total - σ3'
                u = self.sigma3_total - sigma3_prime
                sigma3_total = self.sigma3_total
                sigma1_total = sigma1_prime + u  # σ1_total = σ1' + u
                
                # Proteção contra divisão por zero nas normalizações
                # Usar max() para garantir denominador mínimo
                sigma3_min = max(abs(self.sigma3_total), 1e-10)
                u_ratio = u / sigma3_min
            else:
                # Em CD: u = 0 (drenado)
                u = 0.0
                sigma3_total = sigma3_prime
                sigma1_total = sigma1_prime
                u_ratio = 0.0
            
            # Deformação volumétrica INCREMENTAL (descontando a inicial)
            eps_v_total = np.trace(self.total_strain)
            eps_v_incremental = eps_v_total - eps_v_inicial
            
            # AJUSTE: Armazenar resultados - tensões EFETIVAS
            # O modelo constitutivo trabalha com tensões efetivas
            # As tensões totais podem ser calculadas: σ_total = σ' + u
            self.axial_strain.append((i + 1) * d_eps)
            self.q.append(q)
            self.p.append(p_prime)
            self.volumetric_strain.append(eps_v_incremental)
            self.sigma1.append(sigma1_prime)  # Tensão EFETIVA
            self.sigma3.append(sigma3_prime)  # Tensão EFETIVA
            self.pore_pressure.append(u)
            
            # Proteção contra divisão por zero nas normalizações
            sigma3_min = max(abs(self.sigma3_total), 1e-10)
            self.q_sigma3_ratio.append(q / sigma3_min)
            self.u_sigma3_ratio.append(u_ratio)
        
        # AJUSTE: Retornar dicionário completo ao invés de tupla
        # MANTÉM: Compatibilidade com código antigo via TriaxialResults
        return TriaxialResults({
            'axial_strain': self.axial_strain,
            'q': self.q,
            'p': self.p,
            'volumetric_strain': self.volumetric_strain,
            'sigma1': self.sigma1,
            'sigma3': self.sigma3,
            'pore_pressure': self.pore_pressure,
            'q_sigma3_ratio': self.q_sigma3_ratio,
            'u_sigma3_ratio': self.u_sigma3_ratio
        })
