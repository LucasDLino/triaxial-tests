import numpy as np


class MohrCoulombModel:
    """
    Modelo de Mohr-Coulomb com return mapping completo
    Baseado no algoritmo de integração elastoplástica com:
    - Return to main plane (1-vector)
    - Return to edge (2-vector)
    - Return to apex (multi-vector)
    - Hardening piecewise linear
    
    NOTA IMPORTANTE SOBRE ENSAIOS TRIAXIAIS:
    ----------------------------------------
    Este modelo trabalha com TENSÕES EFETIVAS (σ').
    
    - Em ensaios DRENADOS (CD): σ' = σ_total (poropressão u = 0)
    - Em ensaios NÃO DRENADOS (CU): σ' ≠ σ_total
      * O modelo calcula σ' a partir das deformações
      * A classe TriaxialTest calcula u = σ3_total - σ3'
      * Tensões totais: σ_total = σ' + u
    
    A imposição de condições de contorno (ex: σ3_total constante em CU)
    é feita pela classe TriaxialTest, NÃO por este modelo constitutivo.
    """
    
    def __init__(self, E, nu, phi_deg, cohesion, psi_deg=0, sampling_pairs=None, H=None, 
                 compression_positive=True):
        """
        Inicializa o modelo de Mohr-Coulomb
        
        Parameters:
        -----------
        E : float
            Módulo de Young
        nu : float
            Coeficiente de Poisson
        phi_deg : float
            Ângulo de atrito interno (graus)
        cohesion : float
            Coesão inicial
        psi_deg : float
            Ângulo de dilatância (graus)
        sampling_pairs : array-like, optional
            Pares [eps_plastica, coesao] para hardening piecewise linear
            Se None, usa modelo perfeitamente plástico
        H : array-like, optional
            Módulos de hardening para cada segmento
        compression_positive : bool
            Se True, usa convenção geotécnica (compressão positiva)
            Se False, usa convenção de mecânica dos sólidos (tração positiva)
        """
        self.E = E
        self.nu = nu
        self.phi = np.radians(phi_deg)
        self.psi = np.radians(psi_deg)
        self.compression_positive = compression_positive
        
        # Fator de sinal para convenção de tensões
        # Para compressão positiva: inverte o sinal de sin(phi) no critério
        self.sign_conv = -1.0 if compression_positive else 1.0
        
        # Propriedades elásticas
        self.G = E / (2 * (1 + nu))
        self.K = E / (3 * (1 - 2 * nu))
        
        # Configuração de hardening
        if sampling_pairs is None:
            # Perfeitamente plástico
            self.sampling_pairs = np.array([[0.0, cohesion], [1e100, cohesion]])
            self.H = np.array([0.0, 0.0])
        else:
            self.sampling_pairs = np.array(sampling_pairs)
            # Adiciona ponto no infinito para facilitar computação
            self.sampling_pairs = np.vstack([self.sampling_pairs, [1e100, self.sampling_pairs[-1, 1]]])
            
            if H is None:
                # Calcula H automaticamente a partir dos pares
                self.H = np.zeros(len(self.sampling_pairs))
                for i in range(len(self.sampling_pairs) - 1):
                    delta_c = self.sampling_pairs[i+1, 1] - self.sampling_pairs[i, 1]
                    delta_eps = self.sampling_pairs[i+1, 0] - self.sampling_pairs[i, 0]
                    if delta_eps > 0:
                        self.H[i] = delta_c / delta_eps
                    else:
                        self.H[i] = 0.0
            else:
                self.H = np.array(list(H) + [0.0])
        
        self.n_hard = len(self.sampling_pairs)
        
        # Variáveis de estado
        self.stress = np.zeros((3, 3))
        self.strain = np.zeros((3, 3))
        self.elastic_strain = np.zeros((3, 3))
        self.plastic_strain = np.zeros((3, 3))
        self.equivalent_plastic_strain = 0.0
        self.effective_stress = 0.0
        self.effective_strain = 0.0
        
        # Flags
        self.is_plastic = False
        self.is_fail = False
        
        # Tendência de dilatação (para cálculo de poropressão em ensaios não-drenados)
        # Valor positivo = tendência de contração, negativo = tendência de dilatação
        self.plastic_vol_tendency = 0.0
    
    def plfun(self, eps_p):
        """
        Piecewise linear function para cohesão em função da deformação plástica
        
        Parameters:
        -----------
        eps_p : float
            Deformação plástica equivalente
            
        Returns:
        --------
        float
            Cohesão correspondente
        """
        if self.n_hard == 1:
            return self.sampling_pairs[0, 1]
        
        if eps_p < self.sampling_pairs[0, 0]:
            return self.sampling_pairs[0, 1]
        
        for i in range(1, self.n_hard):
            if eps_p < self.sampling_pairs[i, 0]:
                # Interpolação linear
                eps0 = self.sampling_pairs[i-1, 0]
                eps1 = self.sampling_pairs[i, 0]
                c0 = self.sampling_pairs[i-1, 1]
                c1 = self.sampling_pairs[i, 1]
                c = c0 + ((eps_p - eps0) * (c1 - c0)) / (eps1 - eps0)
                return c
        
        return self.sampling_pairs[-1, 1]
    
    def update(self, strain_total):
        """
        Algoritmo de integração para material elastoplástico com Mohr-Coulomb
        
        IMPORTANTE: Este método calcula e retorna TENSÕES EFETIVAS (σ')
        
        Parameters:
        -----------
        strain_total : ndarray (3x3)
            Deformação total no passo atual
            
        Returns:
        --------
        ndarray (3x3)
            Tensor de TENSÕES EFETIVAS atualizado (σ')
            Para obter tensões totais em ensaios não drenados:
            σ_total = σ' + u*I (onde u é a poropressão)
            
        NOTA SOBRE INCREMENTOS:
        -----------------------
        O método recebe strain_total mas calcula internamente o incremento:
            strain_increment = strain_total - self.strain
        
        Isso permite que a classe TriaxialTest acumule deformações totais
        enquanto o modelo gerencia internamente os incrementos necessários
        para a integração elastoplástica.
        """
        # Flags
        self.is_plastic = False
        self.is_fail = False
        self.plastic_vol_tendency = 0.0  # Reset a cada passo
        
        # Incremento de deformação (calculado internamente)
        # IMPORTANTE: O modelo gerencia os incrementos, não a classe TriaxialTest
        strain_increment = strain_total - self.strain
        
        # PREDITOR ELÁSTICO
        ee_trial = self.elastic_strain + strain_increment
        eps_trial = self.equivalent_plastic_strain
        eps = eps_trial
        
        # Componentes hidrostática e desviadora
        eev_trial = np.trace(ee_trial)
        eed_trial = ee_trial - eev_trial * np.eye(3) / 3.0
        
        p_trial = self.K * eev_trial
        s_trial = 2.0 * self.G * eed_trial
        
        stress_trial = s_trial + p_trial * np.eye(3)
        
        # Decomposição espectral
        eigenvalues, eigenvectors = np.linalg.eigh(stress_trial)
        # Ordenar em ordem decrescente
        ind = np.argsort(eigenvalues)[::-1]
        pstrs_trial = eigenvalues[ind]
        pdirs = eigenvectors[:, ind]
        pstrs = pstrs_trial.copy()
        
        # Verificar admissibilidade plástica: sigma1 >= sigma2 >= sigma3
        cohesion_trial = self.plfun(eps_trial)
        # Critério MC: Φ = (σ1 - σ3) + sign*(σ1 + σ3)*sin(φ) - 2c*cos(φ)
        # sign = -1 para compressão positiva, +1 para tração positiva
        Phi_trial = (pstrs_trial[0] - pstrs_trial[2] + 
                     self.sign_conv * (pstrs_trial[0] + pstrs_trial[2]) * np.sin(self.phi) - 
                     2 * cohesion_trial * np.cos(self.phi))
        
        if Phi_trial > 0:
            # PASSO PLÁSTICO: Return mapping
            self.is_plastic = True
            
            # =====================================================================
            # =============================================================
            # Return mapping de Borst/Crisfield (formulação unificada)
            # =============================================================
            # NOTA: Para ψ = 0 em geometria triaxial (σ₂ = σ₃), o 1-vector
            # return naturalmente viola σ₁ ≥ σ₂ ≥ σ₃ (produz σ₃ > σ₂),
            # então cai automaticamente no 2-vector return (edge).
            # Isso é correto e esperado — não requer tratamento especial.
            #
            # CONVENÇÃO: As fórmulas de de Borst assumem tração positiva.
            # Para compressão positiva, aplicamos sign_conv a sinφ e sinψ:
            #   sc_sinφ = sign_conv * sinφ  (-sinφ para comp+)
            #   sc_sinψ = sign_conv * sinψ  (-sinψ para comp+)
            # O denominador 'a' é invariante (sc² = 1).
            # =============================================================
            sc = self.sign_conv
            sin_phi = np.sin(self.phi)
            sin_psi = np.sin(self.psi)
            cos_phi = np.cos(self.phi)
            sc_sin_phi = sc * sin_phi
            sc_sin_psi = sc * sin_psi
            
            dgamma = 0
            # 'a' é invariante à convenção de sinais (sc² = 1)
            a = (4.0 * self.G * (1.0 + 1.0/3.0 * sin_phi * sin_psi) + 
                 4.0 * self.K * sin_phi * sin_psi)
            
            # 1-VECTOR RETURN (Main Plane)
            self.is_fail = True
            for i in range(self.n_hard):
                c_i = (self.sampling_pairs[i, 1] + 
                       (eps_trial - self.sampling_pairs[i, 0]) * self.H[i])
                dgamma = ((sc_sin_phi + 1) * pstrs_trial[0] + 
                         (sc_sin_phi - 1) * pstrs_trial[2] - 
                         2 * cos_phi * c_i) / \
                        (a + 4 * self.H[i] * cos_phi * cos_phi)
                
                eps = eps_trial + 2 * cos_phi * dgamma
                
                if eps >= self.sampling_pairs[i, 0] and eps <= self.sampling_pairs[i + 1, 0]:
                    self.is_fail = False
                    break
            
            pstrs[0] = pstrs_trial[0] - (2.0 * self.G * (1.0 + sc_sin_psi/3.0) + 
                                          2.0 * self.K * sc_sin_psi) * dgamma  # Sigma1
            pstrs[1] = pstrs_trial[1] + (4.0/3.0 * self.G - 2.0 * self.K) * sc_sin_psi * dgamma  # Sigma2
            pstrs[2] = pstrs_trial[2] + (2.0 * self.G * (1.0 - sc_sin_psi/3.0) - 
                                          2.0 * self.K * sc_sin_psi) * dgamma  # Sigma3
            
            # Tendência de dilatação volumétrica: dεᵥᵖ = -sin(ψ) * dgamma
            # Negativo = tendência a DILATAR (expandir)
            # Positivo = tendência a CONTRAIR
            self.plastic_vol_tendency = -sin_psi * dgamma
            TOL = max(np.abs(pstrs)) * 1e-6
            
            # Verificar validade do 1-vector return (apenas para ψ ≠ 0)
            if not ((pstrs[0] + TOL) >= pstrs[1] and (pstrs[1] + TOL) >= pstrs[2]):
                # 2-VECTOR RETURN (Edge)
                dgammaA = 0
                dgammaB = 0
                
                # Identificar lado da edge
                edgeSide = (pstrs_trial[0] * (1 - np.sin(self.psi)) - 
                           2 * pstrs_trial[1] + 
                           pstrs_trial[2] * (1 + np.sin(self.psi)))
                
                sigmaA = (pstrs_trial[0] - pstrs_trial[2] + 
                         (pstrs_trial[0] + pstrs_trial[2]) * sc_sin_phi)
                
                if edgeSide > 0:  # RIGHT side
                    b = (2.0 * self.G * (1.0 + sc_sin_phi + sc_sin_psi - 
                         1.0/3.0 * sin_phi * sin_psi) + 
                         4.0 * self.K * sin_phi * sin_psi)
                    sigmaB = (pstrs_trial[0] - pstrs_trial[1] + 
                             (pstrs_trial[0] + pstrs_trial[1]) * sc_sin_phi)
                else:  # LEFT side
                    b = (2.0 * self.G * (1.0 - sc_sin_phi - sc_sin_psi - 
                         1.0/3.0 * sin_phi * sin_psi) + 
                         4.0 * self.K * sin_phi * sin_psi)
                    sigmaB = (pstrs_trial[1] - pstrs_trial[2] + 
                             (pstrs_trial[1] + pstrs_trial[2]) * sc_sin_phi)
                
                # Esquema explícito para hardening piecewise linear
                self.is_fail = True
                for i in range(self.n_hard):
                    term1 = 4.0 * self.H[i] * (sigmaA - sigmaB) * np.cos(self.phi) * np.cos(self.phi)
                    term2 = 2.0 * (a - b) * np.cos(self.phi) * (self.sampling_pairs[i, 1] + 
                            (eps_trial - self.sampling_pairs[i, 0]) * self.H[i])
                    denom = (a - b) * (8.0 * self.H[i] * np.cos(self.phi) * np.cos(self.phi) + a + b)
                    
                    dgammaA = (term1 - term2 + a * sigmaA - b * sigmaB) / denom
                    dgammaB = (-term1 - term2 + a * sigmaB - b * sigmaA) / denom
                    
                    eps = eps_trial + 2 * np.cos(self.phi) * (dgammaA + dgammaB)
                    
                    if eps >= self.sampling_pairs[i, 0] and eps <= self.sampling_pairs[i + 1, 0]:
                        self.is_fail = False
                        break
                
                if edgeSide > 0:  # RIGHT side
                    pstrs[0] = pstrs_trial[0] - (2.0 * self.G * (1.0 + sc_sin_psi/3.0) + 
                                                  2.0 * self.K * sc_sin_psi) * (dgammaA + dgammaB)  # Sigma1
                    pstrs[1] = (pstrs_trial[1] + (4.0/3.0 * self.G - 2.0 * self.K) * sc_sin_psi * dgammaA + 
                               (2.0 * self.G * (1.0 - sc_sin_psi/3.0) - 2.0 * self.K * sc_sin_psi) * dgammaB)  # Sigma2
                    pstrs[2] = (pstrs_trial[2] + (2.0 * self.G * (1.0 - sc_sin_psi/3.0) - 2.0 * self.K * sc_sin_psi) * dgammaA + 
                               (4.0/3.0 * self.G - 2.0 * self.K) * sc_sin_psi * dgammaB)  # Sigma3
                else:  # LEFT side
                    pstrs[0] = (pstrs_trial[0] - (2.0 * self.G * (1.0 + sc_sin_psi/3.0) + 
                                2.0 * self.K * sc_sin_psi) * dgammaA + 
                               (4.0/3.0 * self.G - 2.0 * self.K) * sc_sin_psi * dgammaB)  # Sigma1
                    pstrs[1] = (pstrs_trial[1] + (4.0/3.0 * self.G - 2.0 * self.K) * sc_sin_psi * dgammaA - 
                               (2.0 * self.G * (1.0 + sc_sin_psi/3.0) + 2.0 * self.K * sc_sin_psi) * dgammaB)  # Sigma2
                    pstrs[2] = (pstrs_trial[2] + (2.0 * self.G * (1.0 - sc_sin_psi/3.0) - 
                                2.0 * self.K * sc_sin_psi) * (dgammaA + dgammaB))  # Sigma3
            
                
                # Atualizar tendência de dilatação para 2-vector return
                self.plastic_vol_tendency = -sin_psi * (dgammaA + dgammaB)
            TOL = max(np.abs(pstrs)) * 1e-6
            
            # Verificar validade do 2-vector return
            if not ((pstrs[0] + TOL) >= pstrs[1] and (pstrs[1] + TOL) >= pstrs[2]):
                # MULTI-VECTOR RETURN (Apex)
                dEpv = 0
                alpha = np.cos(self.phi) / np.sin(self.psi) if np.sin(self.psi) != 0 else 0
                
                self.is_fail = True
                for i in range(self.n_hard):
                    cot_phi = 1.0 / np.tan(self.phi) if np.tan(self.phi) != 0 else 1e10
                    dEpv = ((p_trial - cot_phi * (self.sampling_pairs[i, 1] + 
                            (eps_trial - self.sampling_pairs[i, 0]) * self.H[i])) / 
                            (cot_phi * alpha * self.H[i] + self.K))
                    
                    eps = eps_trial + alpha * dEpv
                    
                    if eps >= self.sampling_pairs[i, 0] and eps <= self.sampling_pairs[i + 1, 0]:
                        self.is_fail = False
                        break
                
                pstrs[:] = p_trial - self.K * dEpv
            
            # Atualizar tensão
            self.stress = pdirs @ np.diag(pstrs) @ pdirs.T
            
            p = np.trace(self.stress) / 3.0
            s = self.stress - p * np.eye(3)
            
            # Atualizar estado para passo plástico
            # CORREÇÃO: Usar .copy() para evitar referência compartilhada
            self.strain = strain_total.copy()
            self.equivalent_plastic_strain = eps
            self.effective_stress = np.sqrt(3.0 * 0.5 * np.trace(s @ s))
            
            # Calcular deformação elástica a partir da tensão
            self.elastic_strain = ((9.0 * self.K * self.stress + np.trace(self.stress) * np.eye(3) * 
                                   (2.0 * self.G - 3.0 * self.K)) / (18.0 * self.G * self.K))
            
            self.plastic_strain = self.strain - self.elastic_strain
            
            ed = self.strain - (np.trace(self.strain) / 3.0) * np.eye(3)
            self.effective_strain = np.sqrt((2.0 * np.trace(ed @ ed)) / 3.0)
            
        else:
            # PASSO ELÁSTICO
            # CORREÇÃO: Usar .copy() para evitar referência compartilhada
            self.strain = strain_total.copy()
            s = s_trial
            self.stress = s + p_trial * np.eye(3)
            self.effective_stress = np.sqrt(3.0 * 0.5 * np.trace(s @ s))
            self.elastic_strain = ee_trial.copy()
            self.plastic_strain = self.strain - self.elastic_strain
            self.equivalent_plastic_strain = eps_trial
            
            ed = self.strain - (np.trace(self.strain) / 3.0) * np.eye(3)
            self.effective_strain = np.sqrt((2.0 * np.trace(ed @ ed)) / 3.0)
        
        return self.stress
