from numpy import integer
import toml
import os
import shutil
import math
import ccx2paraview
import numpy


def adiciona_cone(arquivo, tag: integer,  x0, y0, z0, dx, dy, dz, r0, r1):
    """Adiociona Cone() para o volume do cone do radier

    Args:
        arquivo (array): array que contem as strings para descrever o arquivo .geo 
        tag (integer): tag 3D do cone
        x0 (float): coordenada X do centro do circulo que está no primeiro plano do tronco de cone
        y0 (float): coordenada Y do centro do circulo que está no primeiro plano do tronco de cone
        z0 (float): coordenada Z do centro do circulo que está no primeiro plano do tronco de cone
        dx (float): valor na direcao de X do vetor que descreve o centro do tronco de cone que está no segundo plano
        dy (float): valor na direcao de Y do vetor que descreve o centro do tronco de cone que está no segundo plano
        dz (float): valor na direcao de Z do vetor que descreve o centro do tronco de cone que está no segundo plano
        r0 (float): valor do raio do circulo no primeiro plano
        r1 (float): valor do raio do circulo no segundo plano
    """
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
    """Adiociona Cylinder() para o volume do cilindro

    Args:
        arquivo (array): array que contem as strings para descrever o arquivo .geo 
        tag (integer): tag 3D do cilindro
        x0 (float): coordenada X do centro do circulo que está no primeiro plano do cilindro
        y0 (float): coordenada Y do centro do circulo que está no primeiro plano do cilindro
        z0 (float): coordenada Z do centro do circulo que está no primeiro plano do cilindro
        dx (float): valor na direcao de X do vetor que descreve o centro do circulo que está no segundo plano
        dy (float): valor na direcao de Y do vetor que descreve o centro do circulo que está no segundo plano
        dz (float): valor na direcao de Z do vetor que descreve o centro do circulo que está no segundo plano
        r0 (float): valor do raio do cilindro
        coment (str, optional): comentario a ser inserido na linha do .geo. Defaults to "".
    """
    arquivo.append('Cylinder(%i) = {%f, %f, %f, %f, %f, %f, %f};%s'
                   % (tag,  x0, y0, z0, dx, dy, dz, r0, coment))
    pass


def adiciona_bloco_solo(arquivo, tag: integer,  x0, y0, z0, dx, dy, dz, coment=""):
    """Adiociona Box() para o volume do solo envolvente

    Args:
        arquivo (array): array que contem as strings para descrever o arquivo .geo 
        tag (integer): tag 3D do solo
        x0 (double): coordenada X do vertice inicial do paralelepipedo
        y0 (_type_): coordenada Y do vertice inicial do paralelepipedo
        z0 (_type_): coordenada Z do vertice inicial do paralelepipedo
        dx (_type_): valor na direcao de X da aresta do paralelepipedo
        dy (_type_): valor na direcao de Y da aresta do paralelepipedo
        dz (_type_): valor na direcao de Z da aresta do paralelepipedo
        coment (str, optional): comentario a ser inserido na linha do .geo. Defaults to "".
    """

    arquivo.append('Box(%i) = {%f, %f, %f, %f, %f, %f}; %s'
                   % (tag,  x0, y0, z0, dx, dy, dz, coment))
    pass


def lateral_cilindro_tag(tag: integer):
    """Baseado na tag do cilindro, retorna a tag do plano que representa a lateral do cilindro

    Args:
        tag (integer): numero do cilindro

    Returns:
        int: tag do plano da lateral da estaca
    """
    return (tag*3 - 2)


def ponta_cilindro_tag(tag: integer):
    """Baseado na tag do cilindro, retorna a tag do plano que representa a ponta de menor cota X,
    assumindo que o vetor da direção do cilindro está na direção de -X
    ** Usar essa funcao apenas para obter o plano da ponta da estaca em contato com o solo

    Args:
        tag (integer): numero do cilindro

    Returns:
        int: tag do plano da ponta da estaca
    """
    return (tag*3 - 1)


def grava_geo(nome_arquivo: str, dados_txt, tem_solo: int = 0):
    """Formata e grava o arquivo .geo a ser interpretado pelo GMSH para geração
    da geometria e malha do radier

    Args:
        nome_arquivo (str): nome do arquivo de configuracao lido, sem a extensao do formato
        dados_txt (dictionary): dicionario {key:value} com o conteudo do arquivo toml lido
        coment (int, optional): comentario a ser inserido na linha do .geo. Defaults to 0.
    """

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

    arquivo.append(
        'Physical Point ("nos_carga") = {%d : %d};' % (1, num_cargas))
    arquivo.append('Physical Volume ("bloco01") = {1,2,3};')

    # definicao das estacas
    # definicao das variaveis de geometria das estacas:
    D_estq = dados_txt['estacas']['diam_est'] / 100 / 2
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

    tag_ultimo_vol = cont_estq-1
    tag_ultimo_plano = n_estq * 3 + 16

    if tem_solo:
        # definicao das variaveis de geometria d solo:
        dx_solo = dados_txt['solo']['h_solo'] / 100
        dy_solo = dados_txt['solo']['base_y'] / 100
        dz_solo = dados_txt['solo']['base_z'] / 100
        y0_solo = -dy_solo/2
        z0_solo = -dz_solo/2

        adiciona_bloco_solo(arquivo, cont_estq,  0, y0_solo, z0_solo, -dx_solo, dy_solo, dz_solo,
                            "// Cubo solo envolvente")

        plano_fundo_solo = 3 * n_estq + 10
        arquivo.append(
            """
Physical Surface ("LatSolo") = {%d:%d};
//retira volume das estacas do volume de solo
b() = BooleanDifference { Volume {%d}; Delete; }{ Volume {1:%d};};
Physical Volume ("solo") = {b()};"""
            % (plano_fundo_solo, plano_fundo_solo + 4, cont_estq, cont_estq - 1))

        tag_ultimo_vol = cont_estq
        tag_ultimo_plano = n_estq * 3 + 22
        pass

    arquivo.append("""
//make all interfaces conformal
b() = BooleanFragments{ Volume {1:3}; Delete; }{ Volume {%d:%d}; Delete; };
Physical Volume (\"bloco\") = {b()};
Point{%i : %i} In Surface {%i};""" % (4, tag_ultimo_vol, 1, num_cargas, tag_ultimo_plano))

    arquivo.append("""
         MeshSize {:} = %f;
         Mesh.SaveAll=1;
         Mesh.SaveGroupsOfElements = 1;
         Mesh.SaveGroupsOfNodes = 1;
         Mesh 3;
         Save \"%s_mesh.inp\";""" % (
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
    """Executa CGX

    Args:
        nome_arquivo (str): nome do arquivo de configuracao lido, sem a extensao do formato
    """
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


def executa_ccx(nome_arquivo: str):
    """Grava arquivo e executa CCX

    Args:
        nome_arquivo (str): nome do arquivo de configuracao lido, sem a extensao do formato
    """

    # header _cgx.geo
    arquivo = []

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

    with open(nome_arquivo + '_sys_solve.fbd', 'w') as file_out:
        file_out.writelines('\n'.join(arquivo))
        file_out.close()
        pass

    executa_cgx(nome_arquivo + '_sys_solve.fbd')


def grava_fbd(nome_arquivo: str, dados_txt, tem_solo: int = 0):
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

send nos_carga abq
send LatEstacas abq
send bloco01 abq
send estacas abq
send PontEstacas abq""" % (nome_arquivo))

    if tem_solo:
        arquivo.append("""# nos e elementos do solo
send LatSolo abq
send solo abq

""" % ())
        pass

    if os.name == 'nt':
        arquivo.append("""
# solve
# sys ccx %s_solve

""" % (nome_arquivo))

    if os.name == 'posix':
        arquivo.append("""
# solve
# sys ccx_2.19_MT %s_solve

""" % (nome_arquivo))

    with open(nome_arquivo + '.fbd', 'w') as file_out:
        file_out.writelines('\n'.join(arquivo))
        file_out.close()
    pass


def grava_solver(nome_arquivo: str, dados_txt, tipo_calculo: str, tem_solo: int = 0, copia_msh: bool = False):
    #####################################################
    # header _cgx.geo
    arquivo = [
        """// ================== RADIER EM TRONCO DE CONE COM %i ESTACAS =========================
*include,input=%s.msh
*include,input=%s.msh
*include,input=%s.msh
*include,input=%s.msh
""" % (dados_txt['estacas']["num_est"], 'bloco01', 'estacas', 'LatEstacas', 'PontEstacas')]

    if tem_solo:
        #####################################################
        # nos e elementos do solo envolvente
        arquivo.append(""" ** nos e elementos do solo envolvente
*include,input=%s.msh
*include,input=%s.msh
""" % ('solo', 'LatSolo'))
        pass

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
    if tem_solo:
        arquivo.append("""** condicoes de contorno
    *boundary
    NLatSolo,1 , 3, 0
    """ % ())
        pass
    else:
        arquivo.append("""** condicoes de contorno
    *boundary
    NLatEstacas,1 , 3, 0
    """ % ())
        pass

    #####################################################
    # definicao dos materiais do bloco e estacas

    arquivo.append("""** definicao do material do bloco
*material, name=mat_radier
*elastic
%e,%e,0
*density
%e
""" % (dados_txt['geometria']['mod_E'], dados_txt['geometria']['poi'], dados_txt['geometria']['den']))

    if tem_solo:
        #####################################################
        # definicao dos materiais do solo envolvente

        arquivo.append("""** definicao do material do solo
*material, name=mat_solo
*elastic
%e,%e,0
*density
%e
""" % (dados_txt['solo']['mod_E'], dados_txt['solo']['poi'], dados_txt['solo']['den']))
        pass

    #####################################################
    # aplicacao dos materiais

    arquivo.append("""** aplicacao do material no bloco
*solid section, elset=Ebloco01, material=mat_radier
""" % ())

    arquivo.append("""** aplicacao do material nas estacas
*solid section, elset=Eestacas, material=mat_radier
""" % ())

    if tem_solo:
        arquivo.append("""** aplicacao do material no solo
*solid section, elset=Esolo, material=mat_solo
""" % ())
        pass

    #####################################################
    # carregamento estatico do peso proprio

    if tipo_calculo == "estatico":
        gravidade = dados_txt['geral']['gravidade']
        arquivo.append("""** carregamento peso proprio
*step
*static
*DLOAD
Ebloco01,GRAV,%.3f,-1.,0.,0.
Eestacas,GRAV,%.3f,-1.,0.,0.""" % (gravidade, gravidade))

        if tem_solo:
            arquivo.append("""Esolo,GRAV,%.3f,-1.,0.,0.
""" % (gravidade))
            pass
        pass

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

    if tipo_calculo == "modal" or (tipo_calculo == "dinamico" and not copia_msh):
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
U
""")
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

        arquivo.append("""*NODE FILE
U, V, RF
*EL FILE, TOTALS = ONLY
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


def converte_resultados(nome_arquivo: str):
    """Converte os resultados do arquido FRD do Calculix para o formato VTU do Paraview

    Args:
        nome_arquivo (str): nome do arquivo de configuracao lido, sem a extensao do formato
    """

    arquivo = nome_arquivo + "_solve.frd"
    c = ccx2paraview.Converter(arquivo, ['vtu'])
    print("Convertendo arquivos para Paraview")
    c.run()
    print("Conversao realizada!")

    pass


def grava_resultados(nome_arquivo: str, NomePastaResultados: str, tipo_calculo: str):
    """Abre o Paraview Server em formato de background, formata e salva o script a ser utilizado
    na geracao dos resultados em formato txt, deleta o script

    Args:
        nome_arquivo (str): nome do arquivo de configuracao lido, sem a extensao do formato
        NomePastaResultados (str): nome da pasta onde estao salvos os arquivos vtu a sererm lidos e processads pelo Paraview
        tipo_calculo (str): tipo de calculo do processamento atual
    """

    os.chdir(NomePastaResultados)
    print(os.getcwd())

    # define caracter para acessar diretorios
    if os.name == 'nt':
        barraNomeArquivo = '\%s' % (nome_arquivo)
        pass
    elif os.name == 'posix':
        barraNomeArquivo = '/%s' % (nome_arquivo)
        pass

    if tipo_calculo == "estatico":
        arquivo = [
            """from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

teste_estatico_solvevtu = XMLUnstructuredGridReader(registrationName='%s.vtu', FileName=[
                                                    %r])
teste_estatico_solvevtu.PointArrayStatus = [
    'RF', 'S', 'S_Mises', 'S_Principal', 'U']

teste_estatico_solvevtu.TimeArray = 'None'

renderView1 = GetActiveViewOrCreate('RenderView')

teste_estatico_solvevtuDisplay = Show(
    teste_estatico_solvevtu, renderView1, 'UnstructuredGridRepresentation')

teste_estatico_solvevtuDisplay.Representation = 'Surface'
teste_estatico_solvevtuDisplay.ColorArrayName = [None, '']
teste_estatico_solvevtuDisplay.SelectTCoordArray = 'None'
teste_estatico_solvevtuDisplay.SelectNormalArray = 'None'
teste_estatico_solvevtuDisplay.SelectTangentArray = 'None'
teste_estatico_solvevtuDisplay.OSPRayScaleArray = 'S'
teste_estatico_solvevtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
teste_estatico_solvevtuDisplay.SelectOrientationVectors = 'None'
teste_estatico_solvevtuDisplay.ScaleFactor = 2.2200000762939456
teste_estatico_solvevtuDisplay.SelectScaleArray = 'None'
teste_estatico_solvevtuDisplay.GlyphType = 'Arrow'
teste_estatico_solvevtuDisplay.GlyphTableIndexArray = 'None'
teste_estatico_solvevtuDisplay.GaussianRadius = 0.11100000381469727
teste_estatico_solvevtuDisplay.SetScaleArray = ['POINTS', 'S']
teste_estatico_solvevtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
teste_estatico_solvevtuDisplay.OpacityArray = ['POINTS', 'S']
teste_estatico_solvevtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
teste_estatico_solvevtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
teste_estatico_solvevtuDisplay.PolarAxes = 'PolarAxesRepresentation'
teste_estatico_solvevtuDisplay.ScalarOpacityUnitDistance = 0.7304860160655032
teste_estatico_solvevtuDisplay.OpacityArrayName = ['POINTS', 'S']

teste_estatico_solvevtuDisplay.ScaleTransferFunction.Points = [
    -6000690.0, 0.0, 0.5, 0.0, 1318560.0, 1.0, 0.5, 0.0]

teste_estatico_solvevtuDisplay.OpacityTransferFunction.Points = [
    -6000690.0, 0.0, 0.5, 0.0, 1318560.0, 1.0, 0.5, 0.0]

renderView1.ResetCamera(False)

materialLibrary1 = GetMaterialLibrary()

renderView1.Update()

generateGlobalIds1 = GenerateGlobalIds(
    registrationName='GenerateGlobalIds1', Input=teste_estatico_solvevtu)

generateGlobalIds1Display = Show(
    generateGlobalIds1, renderView1, 'UnstructuredGridRepresentation')

generateGlobalIds1Display.Representation = 'Surface'
generateGlobalIds1Display.ColorArrayName = [None, '']
generateGlobalIds1Display.SelectTCoordArray = 'None'
generateGlobalIds1Display.SelectNormalArray = 'None'
generateGlobalIds1Display.SelectTangentArray = 'None'
generateGlobalIds1Display.OSPRayScaleArray = 'GlobalPointIds'
generateGlobalIds1Display.OSPRayScaleFunction = 'PiecewiseFunction'
generateGlobalIds1Display.SelectOrientationVectors = 'None'
generateGlobalIds1Display.ScaleFactor = 2.2200000762939456
generateGlobalIds1Display.SelectScaleArray = 'GlobalPointIds'
generateGlobalIds1Display.GlyphType = 'Arrow'
generateGlobalIds1Display.GlyphTableIndexArray = 'GlobalPointIds'
generateGlobalIds1Display.GaussianRadius = 0.11100000381469727
generateGlobalIds1Display.SetScaleArray = ['POINTS', 'GlobalPointIds']
generateGlobalIds1Display.ScaleTransferFunction = 'PiecewiseFunction'
generateGlobalIds1Display.OpacityArray = ['POINTS', 'GlobalPointIds']
generateGlobalIds1Display.OpacityTransferFunction = 'PiecewiseFunction'
generateGlobalIds1Display.DataAxesGrid = 'GridAxesRepresentation'
generateGlobalIds1Display.PolarAxes = 'PolarAxesRepresentation'
generateGlobalIds1Display.ScalarOpacityUnitDistance = 0.7304860160655032
generateGlobalIds1Display.OpacityArrayName = ['POINTS', 'GlobalPointIds']

generateGlobalIds1Display.ScaleTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 21532.0, 1.0, 0.5, 0.0]

generateGlobalIds1Display.OpacityTransferFunction.Points = [
    0.0, 0.0, 0.5, 0.0, 21532.0, 1.0, 0.5, 0.0]

Hide(teste_estatico_solvevtu, renderView1)

renderView1.Update()

SaveData(%r, proxy=generateGlobalIds1, PointDataArrays=['GlobalPointIds', 'RF', 'S', 'S_Mises', 'S_Principal', 'U', 'vtkGhostType'],
         CellDataArrays=['GlobalCellIds'])

layout1 = GetLayout()
layout1.SetSize(989, 494)

renderView1.CameraPosition = [28.888712451925134,
                              6.1795658897651125, 56.597073547304795]
renderView1.CameraFocalPoint = [-1.434999942779541, 0.0, 0.0]
renderView1.CameraViewUp = [0.8822946977651789, -
                            0.024212433113632765, -0.4700742753844009]
renderView1.CameraParallelScale = 16.695219357272986
""" % (nome_arquivo + '_solve', os.getcwd() + '%s_solve.vtu' % (barraNomeArquivo),
                os.getcwd() + '%s_resultados.txt' % (barraNomeArquivo))]
        pass

    if tipo_calculo == "modal":
        arquivo = [
            """from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

teste_modal_solvepvd = PVDReader(registrationName='%s.pvd', FileName=%r)
teste_modal_solvepvd.PointArrays = ['U']

animationScene1 = GetAnimationScene()

animationScene1.UpdateAnimationUsingDataTimeSteps()


renderView1 = GetActiveViewOrCreate('RenderView')

teste_modal_solvepvdDisplay = Show(teste_modal_solvepvd, renderView1, 'UnstructuredGridRepresentation')

teste_modal_solvepvdDisplay.Representation = 'Surface'
teste_modal_solvepvdDisplay.ColorArrayName = [None, '']
teste_modal_solvepvdDisplay.SelectTCoordArray = 'None'
teste_modal_solvepvdDisplay.SelectNormalArray = 'None'
teste_modal_solvepvdDisplay.SelectTangentArray = 'None'
teste_modal_solvepvdDisplay.OSPRayScaleArray = 'U'
teste_modal_solvepvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
teste_modal_solvepvdDisplay.SelectOrientationVectors = 'None'
teste_modal_solvepvdDisplay.ScaleFactor = 2.2200000762939456
teste_modal_solvepvdDisplay.SelectScaleArray = 'None'
teste_modal_solvepvdDisplay.GlyphType = 'Arrow'
teste_modal_solvepvdDisplay.GlyphTableIndexArray = 'None'
teste_modal_solvepvdDisplay.GaussianRadius = 0.11100000381469727
teste_modal_solvepvdDisplay.SetScaleArray = ['POINTS', 'U']
teste_modal_solvepvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
teste_modal_solvepvdDisplay.OpacityArray = ['POINTS', 'U']
teste_modal_solvepvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
teste_modal_solvepvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
teste_modal_solvepvdDisplay.PolarAxes = 'PolarAxesRepresentation'
teste_modal_solvepvdDisplay.ScalarOpacityUnitDistance = 0.7304860160655032
teste_modal_solvepvdDisplay.OpacityArrayName = ['POINTS', 'U']


teste_modal_solvepvdDisplay.ScaleTransferFunction.Points = [-0.000611832, 0.0, 0.5, 0.0, 3.8625099999999957e-05, 1.0, 0.5, 0.0]

teste_modal_solvepvdDisplay.OpacityTransferFunction.Points = [-0.000611832, 0.0, 0.5, 0.0, 3.8625099999999957e-05, 1.0, 0.5, 0.0]

renderView1.ResetCamera(False)

materialLibrary1 = GetMaterialLibrary()

renderView1.Update()

generateGlobalIds1 = GenerateGlobalIds(registrationName='GenerateGlobalIds1', Input=teste_modal_solvepvd)

generateGlobalIds1Display = Show(generateGlobalIds1, renderView1, 'UnstructuredGridRepresentation')

generateGlobalIds1Display.Representation = 'Surface'
generateGlobalIds1Display.ColorArrayName = [None, '']
generateGlobalIds1Display.SelectTCoordArray = 'None'
generateGlobalIds1Display.SelectNormalArray = 'None'
generateGlobalIds1Display.SelectTangentArray = 'None'
generateGlobalIds1Display.OSPRayScaleArray = 'GlobalPointIds'
generateGlobalIds1Display.OSPRayScaleFunction = 'PiecewiseFunction'
generateGlobalIds1Display.SelectOrientationVectors = 'None'
generateGlobalIds1Display.ScaleFactor = 2.2200000762939456
generateGlobalIds1Display.SelectScaleArray = 'GlobalPointIds'
generateGlobalIds1Display.GlyphType = 'Arrow'
generateGlobalIds1Display.GlyphTableIndexArray = 'GlobalPointIds'
generateGlobalIds1Display.GaussianRadius = 0.11100000381469727
generateGlobalIds1Display.SetScaleArray = ['POINTS', 'GlobalPointIds']
generateGlobalIds1Display.ScaleTransferFunction = 'PiecewiseFunction'
generateGlobalIds1Display.OpacityArray = ['POINTS', 'GlobalPointIds']
generateGlobalIds1Display.OpacityTransferFunction = 'PiecewiseFunction'
generateGlobalIds1Display.DataAxesGrid = 'GridAxesRepresentation'
generateGlobalIds1Display.PolarAxes = 'PolarAxesRepresentation'
generateGlobalIds1Display.ScalarOpacityUnitDistance = 0.7304860160655032
generateGlobalIds1Display.OpacityArrayName = ['POINTS', 'GlobalPointIds']

generateGlobalIds1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 21532.0, 1.0, 0.5, 0.0]

generateGlobalIds1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 21532.0, 1.0, 0.5, 0.0]

Hide(teste_modal_solvepvd, renderView1)

renderView1.Update()

SaveData(%r, proxy=generateGlobalIds1, WriteTimeSteps=1,
    WriteTimeStepsSeparately=1, ChooseArraysToWrite=1, 
    PointDataArrays=['GlobalPointIds', 'U', 'vtkGhostType'],
    CellDataArrays=['GlobalCellIds'])

layout1 = GetLayout()
layout1.SetSize(989, 494)

renderView1.CameraPosition = [-1.434999942779541, 0.0, 64.50537421100464]
renderView1.CameraFocalPoint = [-1.434999942779541, 0.0, 0.0]
renderView1.CameraParallelScale = 16.695219357272986
""" % (nome_arquivo + '_solve', os.getcwd() + '%s_solve.pvd' % (barraNomeArquivo),
                os.getcwd() + '%s_resultados.txt' % (barraNomeArquivo))]
        pass

    with open('script_paraview.py', 'w') as file_out:
        file_out.writelines('\n'.join(arquivo))
        file_out.close()
        pass

    # executa pvbatch para tratar exportar resultados para txt

    if os.name == 'nt':
        os.system('"C:\\Program Files\\ParaView 5.11.0\\bin\\pvbatch.exe" %s'
                  % ('script_paraview.py'))
        os.system('del /Q script_paraview.py')
        pass
    elif os.name == 'posix':
        os.system('"pvbatch" %s' % ('script_paraview.py'))
        os.system('rm -f script_paraview.py')
        pass

    return


def ler_arquivo(nome_arquivo: str = None):
    """Gerencia o nome do arquivo de configuracao no formato TOML a ser lido para o processamento

    Args:
        nome_arquivo (str, optional): se estiver definido, nao pergunta o nome do arquivo de configuracao.
        Defaults to None.

    Returns:
        str: nome do arquivo de configuracao, com a extensao do formato lido
        {key:value}: dicionario com o conteudo do arquivo toml lido

    """
    if nome_arquivo == None:
        nome_arquivo = input(
            "Nome do arquivo a ser lido (incluir extensão, se houver): ")
    if nome_arquivo != "" or nome_arquivo != None:
        dados_txt = toml.load(nome_arquivo)
        pass
    else:
        exit()
        pass
    return nome_arquivo, dados_txt


def escolhe_calculo(tipo_calculo: int = None):
    """Gerencia escolha do tipo de calculo a ser executado no processamento atual

    Args:
        tipo_calculo (int, optional): se definido algum inteiro valido, nao pergunta ao usuario 
        qual o tipo de calculo a ser executado. Defaults to None.

    Returns:
        tipo_calculo(str): retorna a string que contem o nome do tipo de calculo
    """

    if tipo_calculo == None:
        print("\n\nQual tipo de calculo deseja executar?\n" +
              "[1] = calculo estatico\n" +
              "[2] = calculo modal\n" +
              "[3] = calculo dinamico\n")

        tipo_calculo = input(
            "\nDigite opcao:  ")
        pass

    tipo_calculo = int(tipo_calculo)

    if tipo_calculo == 1 or tipo_calculo == 2 or tipo_calculo == 3:
        pass
    else:
        print("%s nao e uma opcao valida." % (tipo_calculo))
        exit()
        pass

    tipo_calculo = ["estatico", "modal", "dinamico"][tipo_calculo - 1]

    return tipo_calculo


def escolhe_solo(tem_solo: int = None):
    """Gerenciamento se possui solo envolvente no modelo do radier

    Args:
        tem_solo (int, optional):
        0: se nao possuir solo;
        1: se possuir solo.
        Defaults to None.

    Returns:
        int : retorna 0/1 se foi escolhido gerar o modelo com, ou sem, solo envolvente.
    """

    if tem_solo == None:
        print("\n\nConsiderar solo envolvente?\n" +
              "[0] = Nao\n" +
              "[1] = Sim\n")

        tem_solo = input(
            "\nDigite opcao:  ")
        pass

    tem_solo = int(tem_solo)

    if tem_solo == 1:
        return 1
        pass
    else:
        return 0
        pass


def gerencia_pastas(nome_arquivo: str, tipo_calculo: str, copia_msh: int = None):

    # identifica extensao do nome do arquivo e cria arquivo .geo com o mesmo nome do .txt
    index = nome_arquivo.rfind(".")
    nome_arquivo = nome_arquivo[:index]

    nome_arquivo_bruto = nome_arquivo
    nome_arquivo = nome_arquivo + "_" + tipo_calculo
    NomePastaResultados = "resultados_" + nome_arquivo

    if copia_msh == None:
        copia_msh = input("""Deseja copiar os dados de malha de algum processamento?: 
 [0] = Nao\n [1] = Estatico\n [2] = Modal\n [3] = Dinamico
\n Digite a opção: """)
        copia_msh = int(copia_msh)
        pass

    if copia_msh in [1, 2, 3]:

        if os.name == 'nt':
            os.system("robocopy %s . *.msh" %
                      ("resultados_" + nome_arquivo_bruto + "_" + escolhe_calculo(copia_msh)))
            if (tipo_calculo in ['dinamico']) and (copia_msh in [2, 3]):
                os.system("robocopy %s . *.eig" %
                          ("resultados_" + nome_arquivo_bruto + "_" + escolhe_calculo(copia_msh)))
                if tipo_calculo == 'dinamico' and copia_msh == 2:
                    os.system("ren *.eig %s.eig" %
                              (nome_arquivo_bruto + "_dinamico_solve"))
                    pass
                pass
            pass
        elif os.name == 'posix':
            os.system("cp -f %s/*.msh ." %
                      ("resultados_" + nome_arquivo_bruto + "_" + escolhe_calculo(copia_msh)))
            if (tipo_calculo in ['dinamico']) and (copia_msh in [2, 3]):
                os.system("cp -f %s/*.eig ." %
                          ("resultados_" + nome_arquivo_bruto + "_" + escolhe_calculo(copia_msh)))
                if tipo_calculo == 'dinamico' and copia_msh == 2:
                    os.system("mv *.eig %s.eig" %
                              (nome_arquivo_bruto + "_dinamico_solve"))
                    pass
                pass
            pass
        copia_msh = True
        pass
    else:
        copia_msh = False
        pass

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

    return nome_arquivo, NomePastaResultados, copia_msh


def gerencia_arquivos(nome_arquivo: str, NomePastaResultados: str):
    """Gerencia manipulacao de arquivos apos o processamento terminar

    Args:
        nome_arquivo (str): nome do arquivo de configuracao TOML que gerou o processamento
        NomePastaResultados (str): nome da pasta para onde os arquivos serao copiados
    """

    if os.name == 'nt':
        os.system("move /Y %s*.* %s" % (nome_arquivo, NomePastaResultados))
        os.system("move /Y *.msh %s" % (NomePastaResultados))
        pass
    elif os.name == 'posix':
        os.system("mv -f %s*.* %s" % (nome_arquivo, NomePastaResultados))
        os.system("mv -f *.msh %s" % (NomePastaResultados))
        pass
    return


def main_func():
    """Gerenciamento geral geral do processamento.
    """

    nome_arquivo, dados_txt = ler_arquivo()
    tipo_calculo = escolhe_calculo()
    tem_solo = escolhe_solo()

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
    grava_solver(nome_arquivo, dados_txt, tipo_calculo, tem_solo, copia_msh)
    if not copia_msh:
        executa_cgx(nome_arquivo + '.fbd')
        pass

    executa_ccx(nome_arquivo)

    converte_resultados(nome_arquivo)

    gerencia_arquivos(nome_arquivo, NomePastaResultados)

    grava_resultados(nome_arquivo, NomePastaResultados, tipo_calculo)

    pass


if __name__ == "__main__":
    main_func()
