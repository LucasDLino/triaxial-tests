import numpy as np


def deviatoric(stress):
    """
    Calcula tensor desviador e tensão média
    
    Convenção: COMPRESSÃO POSITIVA (padrão geotécnico)
    p = tr(σ)/3 é positivo em compressão
    """
    p = np.trace(stress) / 3.0
    return stress - p * np.eye(3), p

def q_invariant(stress):
    """
    Calcula invariante q (tensão desviadora de von Mises)
    
    q = √(3J₂) = √(3/2 s:s)
    
    NOTA: q é sempre positivo (norma do desviador)
    """
    s, _ = deviatoric(stress)
    return np.sqrt(1.5 * np.sum(s*s))
