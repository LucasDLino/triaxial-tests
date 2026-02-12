"""
Main - Validação do Modelo de Mohr-Coulomb para Ensaios Triaxiais

Este módulo executa a validação completa do modelo, incluindo:
- Ensaios CD, CU, UU com validação analítica
- Hardening e Softening
- Círculos de Mohr comparativos

Para executar:
    python main.py

Isso chamará automaticamente o módulo validacao.py que contém
todos os testes unificados.
"""

from validacao import main as executar_validacao


if __name__ == "__main__":
    # Executa a validação completa
    executar_validacao()
