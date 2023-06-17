from main import *


nome_arquivo, dados_txt = ler_arquivo("teste.txt")
tipo_calculo = escolhe_calculo(2)
tem_solo = escolhe_solo(0)

# depuração de leitura do arquivo toml formatado e imprime na tela
new_toml_string = toml.dumps(dados_txt)
print(new_toml_string)

nome_arquivo, NomePastaResultados, copia_msh = gerencia_pastas(
    nome_arquivo, tipo_calculo)

grava_geo(nome_arquivo, dados_txt, tem_solo)
if not copia_msh:
    executa_gmsh(nome_arquivo + '.geo', dados_txt)
    pass
grava_fbd(nome_arquivo, dados_txt, tem_solo)
grava_solver(nome_arquivo, dados_txt, tipo_calculo, tem_solo)
if not copia_msh:
    executa_cgx(nome_arquivo + '.fbd')
    pass

executa_ccx(nome_arquivo)
converte_resultados(nome_arquivo)

gerencia_arquivos(nome_arquivo, NomePastaResultados)

grava_resultados(nome_arquivo, NomePastaResultados, tipo_calculo)
