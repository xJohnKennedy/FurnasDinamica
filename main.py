from numpy import integer
import toml
import os
import math
import ccx2paraview
import numpy


def adiciona_cone(arquivo, tag: integer,  x0, y0, z0, dx, dy, dz, r0, r1):
    # definicao do tronco de cone:
    arquivo.append('Cone(%i) = {%f, %f, %f, %f, %f, %f, %f, %f};'
                   % (tag,  x0, y0, z0, dx, dy, dz, r0, r1))
    pass


def cria_pontos_cargas(arquivo,  x0, y0, z0, r0, num_points):
    # distribuição dos pontos de aplicacao das cargas sobre um circulo:
    num_div = 2 * numpy.pi / num_points
    for i in range(0, num_points):
        dy = r0 * numpy.cos(i*num_div)
        dz = r0 * numpy.sin(i*num_div)
        arquivo.append('Point(%i) = {%f, %f, %f};   // ponto %d aplicacao carga'
                       % (i + 1, x0, dy + y0, dz + z0, i + 1))
        pass
    pass


def adiciona_cilindro(arquivo, tag: integer,  x0, y0, z0, dx, dy, dz, r0, coment=""):
    # definicao do tronco de cone:
    arquivo.append('Cylinder(%i) = {%f, %f, %f, %f, %f, %f, %f};%s'
                   % (tag,  x0, y0, z0, dx, dy, dz, r0, coment))
    pass


def lateral_cilindro_tag(tag: integer):
    return (tag*3 - 2)


def ponta_cilindro_tag(tag: integer):
    return (tag*3 - 1)


def grava_geo(nome_arquivo, dados_txt):

    # header .geo
    arquivo = [
        '// ================== RADIER EM TRONCO DE CONE COM {:n} ESTACAS ========================= '.format(
            dados_txt['estacas']["num_est"]),
        "//+\nSetFactory(\"OpenCASCADE\");\n//+"]

    # definicao das variaveis de geometria:
    h_base = dados_txt['geometria']['h_base'] / 100
    h_cone = dados_txt['geometria']['h_cone'] / 100
    D_base = dados_txt['geometria']['diam_base'] / 100
    D_topo = dados_txt['geometria']['diam_topo'] / 100
    h_topo = dados_txt['geometria']['h_topo'] / 100

    num_cargas = dados_txt['cargas']['num_cargas']
    diam_circ_carga = dados_txt['cargas']['diam_carga']/100

    # definicao dos pontos de aplicação do carregamento:
    cria_pontos_cargas(arquivo, h_base + h_cone + h_topo,
                       0, 0, diam_circ_carga/2, num_cargas)

    # definicao do tronco de cone:
    adiciona_cone(arquivo, 1, h_base, 0, 0, h_cone, 0, 0, D_base/2, D_topo/2)

    # definicao dos cilindros de base e topo
    adiciona_cilindro(arquivo, 2,  0, 0, 0, h_base, 0, 0, D_base/2)

    adiciona_cilindro(arquivo, 3,  h_base+h_cone, 0, 0, h_topo, 0, 0, D_topo/2)

    # adiciona pontos de cargas na geometria
    arquivo.append('Point{%i : %i} In Surface {%i};' % (1, num_cargas, 8))
    arquivo.append('Point{%i : %i} In Volume {%i};' % (1, num_cargas, 3))

    arquivo.append(
        'Physical Point ("nos_carga") = {%i : %i};' % (1, num_cargas))
    arquivo.append('Physical Volume ("bloco01") = {1,2,3};')

    # definicao das estacas
    # definicao das variaveis de geometria das estacas:
    D_estq = dados_txt['estacas']['diam_est'] / 100
    H_estq = dados_txt['estacas']['h_est'] / 100
    n_estq = dados_txt['estacas']['num_est']
    cob_est = dados_txt['estacas']['cob_est']/100
    Proj_estq = (D_base - 2*cob_est)

    y0 = Proj_estq/2
    z0 = 0.0

    adiciona_cilindro(arquivo, 4,  0, y0, z0, -H_estq,
                      0, 0, D_estq, "  // estaca %i" % (1))
    grupo_estacas = "Physical Volume (\"estacas\") = { 4"
    grupo_Lateral_estacas = "Physical Surface (\"LatEstacas\") = { %i" % (
        lateral_cilindro_tag(4))

    grupo_ponta_estacas = "Physical Surface (\"PontEstacas\") = { %i" % (
        ponta_cilindro_tag(4))

    cont_estq = 5
    Dalpha = 360/n_estq
    alpha = Dalpha

    for i in range(1, n_estq):
        y1 = (Proj_estq/2)*math.cos(math.pi*alpha/180)
        z1 = (Proj_estq/2)*math.sin(math.pi*alpha/180)

        adiciona_cilindro(arquivo, cont_estq,  0, y1, z1, -H_estq,
                          0, 0, D_estq, "  // estaca %i" % (i+1))
        grupo_estacas += ", %i" % (cont_estq)
        grupo_Lateral_estacas += ", %i" % (lateral_cilindro_tag(cont_estq))
        grupo_ponta_estacas += ", %i" % (ponta_cilindro_tag(cont_estq))
        alpha = alpha + Dalpha
        cont_estq = cont_estq + 1

    grupo_estacas += "};"
    grupo_Lateral_estacas += "};"
    grupo_ponta_estacas += "};"

    arquivo.append(grupo_estacas)
    arquivo.append(grupo_Lateral_estacas)
    arquivo.append(grupo_ponta_estacas)

    arquivo.append(
        """//make all interfaces conformal
        b() = BooleanFragments{ Physical Volume {1}; Delete; }{ Physical Volume {2}; Delete; };
        Physical Volume (\"bloco\") = {b()};
        Mesh.MeshSizeFactor = %f;
        Mesh 3;
        Save \"%s_mesh.inp\";""" % (dados_txt['geral']['MeshSizeFactor'], nome_arquivo))

    with open(nome_arquivo + '.geo', 'w') as file_out:
        file_out.writelines('\n'.join(arquivo))
        file_out.close()
    pass


def executa_gmsh(nome_arquivo):
    os.system("gmsh %s -" % (nome_arquivo))
    pass


def executa_cgx(nome_arquivo):
    if os.name == 'nt':
        comando = "cmd /c \"C:\\Program Files (x86)\\bConverged\\common\\site\\cmdStartup_mod.bat\" cgx -bg % s" % (
            nome_arquivo)
        os.system(comando)
        pass
    elif os.name == 'posix':
        comando = "bash -i -c cgx -bg % s" % (
            nome_arquivo)
        os.system(comando)

        comando = "bash -i -c ccx %s_solve" % (
            nome_arquivo)
        os.system(comando)

        pass
    pass


def grava_fbd(nome_arquivo: str, dados_txt):
    # header _cgx.geo
    arquivo = [
        """#
# // ================== RADIER EM TRONCO DE CONE COM %i ESTACAS =========================
# // Arquivo .fbd geracao da malha, condicoes de contorno e carregamentos entendidos pelo Calculix
#""" % (dados_txt['estacas']["num_est"])]

    arquivo.append("""#
read %s_mesh.inp
ulin gmsh without virtual topology
# remove elementos de superficie
zap +CPS3
# salva definicoes de malha e condicoes de contorno
send all abq
send nos_carga abq
send LatEstacas abq

""" % (nome_arquivo))

    if os.name == 'nt':
        arquivo.append("""
# solve
sys ccx %s_solve

""" % (nome_arquivo))

    with open(nome_arquivo + '.fbd', 'w') as file_out:
        file_out.writelines('\n'.join(arquivo))
        file_out.close()
    pass


def grava_solver(nome_arquivo: str, dados_txt):
    # header _cgx.geo
    arquivo = [
        """// ================== RADIER EM TRONCO DE CONE COM %i ESTACAS =========================
*include,input=%s.msh
*include,input=%s.msh
""" % (dados_txt['estacas']["num_est"], 'all', 'LatEstacas')]

    arquivo.append("""** condicoes de contorno
*boundary
NLatEstacas,1 , 3, 0
""" % ())

    arquivo.append("""** definicao do material
*material, name=mat_radier
*elastic
%e,%e,0
*density
%e
""" % (dados_txt['geometria']['mod_E'], dados_txt['geometria']['poi'], dados_txt['geometria']['den']))

    arquivo.append("""** aplicacao do material
*solid section, elset=Eall, material=mat_radier
""" % ())

    arquivo.append("""** frequencia natural
*step
*frequency,solver=arpack,storage=yes
%d""" % (dados_txt['freq_natural']['num_modos']))

    menor_freq = dados_txt['freq_natural']['menor_freq']
    maior_freq = dados_txt['freq_natural']['maior_freq']

    if (menor_freq == 0) and (maior_freq == 0):
        arquivo.append(arquivo.pop() + """\n""" % ())
    elif ((maior_freq == 0)):
        arquivo.append(arquivo.pop() + (""",%.2f\n""" % (menor_freq)))
    else:
        arquivo.append(arquivo.pop() + (""",%.2f,%.2f\n""" %
                       (menor_freq, maior_freq)))

    arquivo.append("""** gravacao dos resultados
*node file
U
*end step
""" % ())

    arquivo.append("""** carregamento peso proprio
*step
*static
*DLOAD
Eall,GRAV,%.3f,-1.,0.,0.
""" % (dados_txt['geral']['gravidade']))

    arquivo.append("""** gravacao dos resultados
*el file
S
*end step
""" % ())

    arquivo.append("""** pontos de aplicacao do carregamento dinamico
*include,input=nos_carga.msh
**
** pontos de leitura do carregamento dinamico
*NSET,NSET=NO_DESLOC_CENTRO
8591
""")
    arquivo.append("""*Amplitude, name=%s
**
** LER UM ARQUIVO EXTERNO (t, f(t))
**
*include,input=%s
**

""" % (dados_txt['cargas']['car_arq_0'], dados_txt['cargas']['car_arq_0']))

    arquivo.append("""** calculo dinamico - vibracao forcada
*STEP, INC=%d
*MODAL DYNAMIC,SOLVER=ITERATIVE SCALING
%e,%f
*MODAL DAMPING
1,%d, 0.02
*CLOAD,AMPLITUDE=%s
Nnos_carga, 1, -1.
*NODE FILE,NSET=NO_DESLOC_CENTRO
U
*EL FILE,ELSET=Eall,TOTALS=ONLY
ELSE,ELKE,EVOL
*END STEP
""" % (dados_txt['cargas']['tempo_final'] / dados_txt['cargas']['incremento_tempo'] + 1,
       dados_txt['cargas']['incremento_tempo'], dados_txt['cargas']['tempo_final'],
        dados_txt['freq_natural']['num_modos'],
        dados_txt['cargas']['car_arq_0']))

    with open(nome_arquivo + '_solve.inp', 'w') as file_out:
        file_out.writelines('\n'.join(arquivo))
        file_out.close()
    pass


def converte_resultados(nome_arquivo: str):
    arquivo = nome_arquivo + "_solve.frd"
    c = ccx2paraview.Converter(arquivo, ['vtu'])
    print("Convertendo arquivos para Paraview")
    c.run()
    print("Conversao realizada!")

    pass


def main_func():
    nome_arquivo = input(
        "Nome do arquivo a ser lido (incluir extensão, se houver): ")
    if nome_arquivo != "":
        dados_txt = toml.load(nome_arquivo)
    else:
        nome_arquivo = None

    # depuração de leitura do arquivo toml formatado
    new_toml_string = toml.dumps(dados_txt)
    print(new_toml_string)

    # identifica extensao do nome do arquivo e cria arquivo .geo com o mesmo nome do .txt
    index = nome_arquivo.rfind(".")
    nome_arquivo = nome_arquivo[:index]

    grava_geo(nome_arquivo, dados_txt)
    executa_gmsh(nome_arquivo + '.geo')
    grava_fbd(nome_arquivo, dados_txt)
    grava_solver(nome_arquivo, dados_txt)
    executa_cgx(nome_arquivo + '.fbd')
    converte_resultados(nome_arquivo)

    pass


if __name__ == "__main__":
    main_func()
