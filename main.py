"""
Modelo de Mohr-Coulomb — Simulação de Ensaios Triaxiais
========================================================

Módulos:
  mohr_courlomb.py  — Modelo constitutivo (return mapping de Borst/Crisfield)
  triaxial.py       — Simulação de ensaios (CD, CU, UU)
  validacao.py      — Verificação numérica (PASS/FAIL) + gráficos
  utils.py          — Funções auxiliares

Uso:
  python main.py                   → validação completa (testes + gráficos)
  python main.py --check           → apenas testes PASS/FAIL
  python main.py --plot            → apenas gráficos (salva PNG)
  python main.py --plot --show     → gráficos salvos + exibe na tela

Documentação do algoritmo: ALGORITMO_TRIAXIAL.md
"""

from validacao import main as executar_validacao


if __name__ == "__main__":
    executar_validacao()
