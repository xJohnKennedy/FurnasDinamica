from numpy import integer
import toml
import os
import shutil
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
    arquivo.append('//Point{%i : %i} In Surface {%i};' % (1, num_cargas, 8))
    arquivo.append('//Point{%i : %i} In Volume {%i};' % (1, num_cargas, 3))

    arquivo.append(
        'Physical Point ("nos_carga") = {%d : %d};' % (1, num_cargas))
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
        b() = BooleanFragments{ Volume {1:3}; Delete; }{ Volume {%d:%d}; Delete; };
        Physical Volume (\"bloco\") = {b()};
        Point{%i : %i} In Surface {%i};
        MeshSize {:} = %f;
        Mesh.SaveAll=1;
        Mesh.SaveGroupsOfElements = 1;
        Mesh.SaveGroupsOfNodes = 1;
        Mesh 3;
        Save \"%s_mesh.inp\";""" % (4, cont_estq-1, 1, num_cargas, n_estq * 3 + 16,
                                    dados_txt['geral']['TamanhoMalha'], nome_arquivo))

    with open(nome_arquivo + '.geo', 'w') as file_out:
        file_out.writelines('\n'.join(arquivo))
        file_out.close()
    pass


def executa_gmsh(nome_arquivo, dados_txt):
    os.system("gmsh %s -clmin %f -clmax %f -clscale %f -" % (nome_arquivo,
                                                             dados_txt['geral']['MalhaMin'],
                                                             dados_txt['geral']['MalhaMax'],
                                                             dados_txt['geral']['FatorDeMalha']))
    pass


def executa_cgx(nome_arquivo):
    if os.name == 'nt':
        comando = "cmd /c \"C:\\Program Files (x86)\\bConverged\\common\\site\\cmdStartup_mod.bat\" cgx -bg % s" % (
            nome_arquivo)
        os.system(comando)
        pass
    elif os.name == 'posix':
        comando = "cgx_2.20.1 -bg % s" % (
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
# """ % (dados_txt['estacas']["num_est"])]

    arquivo.append("""#
read %s_mesh.inp
ulin gmsh without virtual topology
# remove elementos de barra
zap +T3D2
# remove elementos de superficie
zap +CPS3
# salva definicoes de malha e condicoes de contorno
send all abq
send nos_carga abq
send nos_carga lst
send LatEstacas abq

""" % (nome_arquivo))

    if os.name == 'nt':
        arquivo.append("""
# solve
sys ccx %s_solve

""" % (nome_arquivo))

    if os.name == 'posix':
        arquivo.append("""
# solve
sys ccx_2.19_MT %s_solve

""" % (nome_arquivo))

    with open(nome_arquivo + '.fbd', 'w') as file_out:
        file_out.writelines('\n'.join(arquivo))
        file_out.close()
    pass


def grava_solver(nome_arquivo: str, dados_txt, tipo_calculo: str):
    #####################################################
    # header _cgx.geo
    arquivo = [
        """// ================== RADIER EM TRONCO DE CONE COM %i ESTACAS =========================
*include,input=%s.msh
*include,input=%s.msh
""" % (dados_txt['estacas']["num_est"], 'all', 'LatEstacas')]

    #####################################################
    # pontos de aplicacao do carregamento estatico e dinamico
    arquivo.append("""** pontos de aplicacao do carregamento dinamico
*include,input=nos_carga.msh
**
** pontos de leitura do carregamento dinamico
*NSET,NSET=NO_DESLOC_CENTRO
8591
""")

    #####################################################
    # condicoes de contorno

    arquivo.append("""** condicoes de contorno
*boundary
NLatEstacas,1 , 3, 0
""" % ())

    #####################################################
    # definicao dos materiais

    arquivo.append("""** definicao do material
*material, name=mat_radier
*elastic
%e,%e,0
*density
%e
""" % (dados_txt['geometria']['mod_E'], dados_txt['geometria']['poi'], dados_txt['geometria']['den']))

    #####################################################
    # aplicacao dos materiais

    arquivo.append("""** aplicacao do material
*solid section, elset=Eall, material=mat_radier
""" % ())

    #####################################################
    # carregamento estatico do peso proprio

    if tipo_calculo == "estatico":
        arquivo.append("""** carregamento peso proprio
*step
*static
*DLOAD
Eall,GRAV,%.3f,-1.,0.,0.
""" % (dados_txt['geral']['gravidade']))

    ##########
    # se existir carregamento estatico aplicado
        arquivo.append("""** carregamento estatico""")

        for i in range(0, dados_txt['cargas']['num_cargas'] + 1):
            try:
                carga = dados_txt['cargas']['car_est_%i' % (i)]
                pass
            except:
                carga = 0
                pass
            if carga != 0 and i == 0:
                arquivo.append("""*CLOAD
%s, 1, %f
""" % ('Nnos_carga', carga * (-1)))
                pass
            elif carga != 0 and i != 0:
                arquivo.append("""*CLOAD
%i, 1, %f
""" % (i, carga * (-1)))
                pass
            pass

        arquivo.append("""** gravacao dos resultados
*node file
U, RF
*el file
S
*end step
""" % ())
        pass

    #####################################################
    # calculo da frequencia natural

    if tipo_calculo == "modal" or tipo_calculo == "dinamico":
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

        if tipo_calculo == "modal":
            arquivo.append("""** gravacao dos resultados
*node file
U""")
            pass

        arquivo.append("""*end step
    """ % ())
        pass

    #####################################################
    # inclusao das forcas dinamicas

    if tipo_calculo == "dinamico":

        for i in range(0, dados_txt['cargas']['num_cargas'] + 1):
            try:
                carga = dados_txt['cargas']['car_est_%i' % (i)]
                carga_din = dados_txt['cargas']['car_din_%i' % (i)]
                pass
            except:
                carga = 0
                carga_din = 0
                pass
            if carga != 0 and carga_din != 0:
                arquivo.append("""*Amplitude, name=%s
**
** LER UM ARQUIVO EXTERNO (t, f(t))
**
*include,input=%s
**

""" % ('car_din_%d' % (i), dados_txt['cargas']['car_din_%d' % (i)]))
                pass
            pass
        pass

    #####################################################
    # calculo da vibracao forcada
        num_modos_ini = dados_txt['vib_forcada']['num_modos_ini']
        if num_modos_ini == 0:
            num_modos_ini = 1
            pass

        num_modos_ult = dados_txt['vib_forcada']['num_modos_ult']
        if num_modos_ult == 0:
            num_modos_ult = dados_txt['freq_natural']['num_modos']
            pass

        arquivo.append("""** calculo dinamico - vibracao forcada
*STEP, INC=%d
*MODAL DYNAMIC,SOLVER=ITERATIVE SCALING
%e,%f
*MODAL DAMPING
%d,%d, %f
""" % (dados_txt['cargas']['tempo_final'] / dados_txt['cargas']['incremento_tempo'] + 1,
            dados_txt['cargas']['incremento_tempo'], dados_txt['cargas']['tempo_final'],
            num_modos_ini, num_modos_ult, dados_txt['vib_forcada']['amortecimento']))

    ##########
    # se existir carregamento estatico aplicado e configuração de carga dinamica
        arquivo.append("""** carregamento dinamico""")

        for i in range(0, dados_txt['cargas']['num_cargas'] + 1):
            try:
                carga = dados_txt['cargas']['car_est_%i' % (i)]
                carga_din = dados_txt['cargas']['car_din_%d' % (i)]
                pass
            except:
                carga = 0
                carga_din = 0
                pass
            if carga != 0 and carga_din != 0 and i == 0:
                arquivo.append("""*CLOAD, AMPLITUDE = %s
%s, 1, %f
""" % ('car_din_%d' % (i), 'Nnos_carga', carga * (-1)))
                pass
            elif carga != 0 and carga_din != 0 and i != 0:
                arquivo.append("""*CLOAD, AMPLITUDE = %s
%i, 1, %f
""" % ('car_din_%d' % (i), i, carga * (-1)))
                pass
            pass

        arquivo.append("""*NODE FILE, NSET = Nall
U, RF
*EL FILE, ELSET = Eall, TOTALS = ONLY
ELSE, ELKE, EVOL
*END STEP
""" % ())
        pass

    #####################################################
    # salva inp

    with open(nome_arquivo + '_solve.inp', 'w') as file_out:
        file_out.writelines('\n'.join(arquivo))
        file_out.close()
    pass


def converte_resultados(nome_arquivo: str, NomePastaResultados: str):
    arquivo = nome_arquivo + "_solve.frd"
    c = ccx2paraview.Converter(arquivo, ['vtu'])
    print("Convertendo arquivos para Paraview")
    c.run()
    print("Conversao realizada!")

    pass


def grava_resultados(nome_arquivo: str, NomePastaResultados: str):
    pass


def main_func():
    nome_arquivo = input(
        "Nome do arquivo a ser lido (incluir extensão, se houver): ")
    if nome_arquivo != "":
        dados_txt = toml.load(nome_arquivo)
        pass
    else:
        nome_arquivo = None
        pass

    print("\n\nQual tipo de calculo deseja executar?\n" +
          "[1] = calculo estatico\n" +
          "[2] = calculo modal\n" +
          "[3] = calculo dinamico\n")

    tipo_calculo = input(
        "\nDigite opcao:  ")

    tipo_calculo = int(tipo_calculo)

    if tipo_calculo == 1 or tipo_calculo == 2 or tipo_calculo == 3:
        pass
    else:
        print("%s nao e uma opcao valida." % (tipo_calculo))
        exit()
        pass

    tipo_calculo = ["estatico", "modal", "dinamico"][tipo_calculo - 1]

    # depuração de leitura do arquivo toml formatado
    new_toml_string = toml.dumps(dados_txt)
    print(new_toml_string)

    # identifica extensao do nome do arquivo e cria arquivo .geo com o mesmo nome do .txt
    index = nome_arquivo.rfind(".")
    nome_arquivo = nome_arquivo[:index]

    nome_arquivo = nome_arquivo + "_" + tipo_calculo
    NomePastaResultados = "resultados_" + nome_arquivo

    if os.path.exists(NomePastaResultados):

        apagar_pasta = input(
            "Deseja apagar a pasta %s [S/n]: " % (NomePastaResultados))

        if apagar_pasta.lower() == 'n':
            print("Altere o nome da pasta %s." % (NomePastaResultados))
            exit()
            pass

    shutil.rmtree(NomePastaResultados, ignore_errors=True)

    os.mkdir(NomePastaResultados)

    print(os.getcwd())

    grava_geo(nome_arquivo, dados_txt)
    executa_gmsh(nome_arquivo + '.geo', dados_txt)
    grava_fbd(nome_arquivo, dados_txt)
    grava_solver(nome_arquivo, dados_txt, tipo_calculo)
    executa_cgx(nome_arquivo + '.fbd')
    converte_resultados(nome_arquivo, NomePastaResultados)

    if os.name == 'nt':
        os.system("move /Y %s*.* %s" % (nome_arquivo, NomePastaResultados))
        os.system("move /Y *.msh %s" % (NomePastaResultados))
        pass
    elif os.name == 'posix':
        os.system("mv -f %s*.* %s" % (nome_arquivo, NomePastaResultados))
        os.system("mv -f *.msh %s" % (NomePastaResultados))
        pass
    pass

    pass


if __name__ == "__main__":
    main_func()
