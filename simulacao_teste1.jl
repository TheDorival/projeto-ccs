using JutulDarcy, Jutul, Statistics, CSV, DataFrames

# ═══════════════════════════════════════════════
# CONFIGURAÇÃO DE CAMINHOS
# ═══════════════════════════════════════════════

const PASTA_SCRIPT = @__DIR__
const CAMINHO_CSV_PADRAO = joinpath(PASTA_SCRIPT, "malha_1PIR_enriquecida.csv")
const CAMINHO_SAIDA_PADRAO = joinpath(PASTA_SCRIPT, "sensibilidade_campo_pilar.csv")

if !isfile(CAMINHO_CSV_PADRAO)
    error("Arquivo não encontrado: $CAMINHO_CSV_PADRAO")
end

# Gradiente de pressão calibrado no ponto original do Campo Pilar
# (21,62 MPa a 2207,7 m ≈ gradiente hidrostático normal de água)
const GRADIENTE_PRESSAO_PA_POR_M = 21.62e6 / 2207.7   # ≈ 9793.7 Pa/m

# ═══════════════════════════════════════════════
# 1. ENCONTRAR A MELHOR JANELA DE RESERVATÓRIO (restrita à faixa típica de CCS)
# ═══════════════════════════════════════════════

function encontrar_melhor_janela_reservatorio(caminho_csv, n_camadas, dz;
                                                profundidade_min = 800.0,
                                                profundidade_max = 2500.0)
    df = CSV.read(caminho_csv, DataFrame)
    sort!(df, :DEPT_TOPO)
    df_ok = filter(row -> row.COBERTURA_REAL, df)
    df_ok = filter(row -> profundidade_min <= row.DEPT_TOPO <= profundidade_max, df_ok)

    n = nrow(df_ok)
    if n < n_camadas
        error("Poucos blocos válidos ($n) no intervalo $profundidade_min–$profundidade_max m. " *
              "Considere ampliar a faixa de busca.")
    end

    melhor_kx_medio = -Inf
    melhor_indice_inicio = 1

    for inicio in 1:(n - n_camadas + 1)
        janela = df_ok[inicio:(inicio + n_camadas - 1), :]
        if maximum(janela.DEPT_TOPO) - minimum(janela.DEPT_TOPO) > (n_camadas * dz * 1.5)
            continue
        end
        kx_medio = mean(janela.KX)
        if kx_medio > melhor_kx_medio
            melhor_kx_medio = kx_medio
            melhor_indice_inicio = inicio
        end
    end

    melhor_janela = df_ok[melhor_indice_inicio:(melhor_indice_inicio + n_camadas - 1), :]
    prof_injecao_otima = (minimum(melhor_janela.DEPT_TOPO) + maximum(melhor_janela.DEPT_TOPO)) / 2 + dz/2

    println("[INFO] Melhor janela de reservatório (faixa $profundidade_min–$profundidade_max m): " *
            "$(minimum(melhor_janela.DEPT_TOPO))–$(maximum(melhor_janela.DEPT_TOPO))m | " *
            "prof_injecao=$prof_injecao_otima m | KX_médio=$(round(melhor_kx_medio,digits=1))mD | " *
            "arenito=$(sum(melhor_janela.FACIES.==1))/$n_camadas")

    return prof_injecao_otima
end

const PROF_INJECAO_OTIMA = encontrar_melhor_janela_reservatorio(CAMINHO_CSV_PADRAO, 12, 10.0)

# ═══════════════════════════════════════════════
# 2. CARREGAR DADOS REAIS DO POÇO 1PIR
# ═══════════════════════════════════════════════

function carregar_camadas_reservatorio(caminho_csv, prof_injecao, n_camadas, dz, cenario_geomec)
    df = CSV.read(caminho_csv, DataFrame)
    sort!(df, :DEPT_TOPO)

    df_ok = filter(row -> row.COBERTURA_REAL, df)

    meia_janela = (n_camadas / 2) * dz
    topo_janela = prof_injecao - meia_janela
    base_janela = prof_injecao + meia_janela

    df_zona = filter(row -> topo_janela <= row.DEPT_TOPO < base_janela, df_ok)

    if nrow(df_zona) < n_camadas
        error("Cobertura insuficiente: $(nrow(df_zona)) blocos válidos de $n_camadas " *
              "necessários entre $(topo_janela)m e $(base_janela)m.")
    end

    sort!(df_zona, :DEPT_TOPO)
    df_zona = df_zona[1:n_camadas, :]

    coesao_col = Symbol("COESAO_", cenario_geomec)
    atrito_col = Symbol("ATRITO_", cenario_geomec)
    pc_col     = Symbol("PC_ENTRADA_PA_", cenario_geomec)
    swirr_col  = Symbol("SWIRR_", cenario_geomec)
    sgr_col    = Symbol("SGR_RESIDUAL_", cenario_geomec)

    for col in [coesao_col, atrito_col, pc_col, swirr_col, sgr_col]
        if !(String(col) in names(df_zona))
            error("Coluna $col não encontrada no CSV.")
        end
    end

    return (
        kx     = df_zona.KX,
        kz     = df_zona.KZ,
        phi    = df_zona.PHIE,
        coesao = df_zona[!, coesao_col],
        atrito = df_zona[!, atrito_col],
        pc_entrada = df_zona[!, pc_col],
        swirr      = df_zona[!, swirr_col],
        sgr        = df_zona[!, sgr_col],
        facies = df_zona.FACIES,
        dept   = df_zona.DEPT_TOPO,
    )
end

# ═══════════════════════════════════════════════
# 3. CONSTRUÇÃO DO DOMÍNIO
# ═══════════════════════════════════════════════

function build_domain(nx, ny, nz, Lx, Lz, prof_injecao,
                      kx_reservatorio, phi_reservatorio,
                      k_selante_mD, k_FSASP_mD,
                      k_falha_mD, phi_selante,
                      phi_FSASP, phi_falha, i_falha, w_falha)

    n_camadas_res = length(kx_reservatorio)
    g  = CartesianMesh((nx, ny, nz), (Lx, 1.0, Lz))
    nc = nx * ny * nz
    perm = zeros(nc)
    poro = zeros(nc)

    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                idx = (k-1)*(nx*ny) + (j-1)*nx + i
                if k <= n_camadas_res
                    perm[idx] = kx_reservatorio[k] * 1e-3 * si_unit(:darcy)
                    poro[idx] = phi_reservatorio[k]
                elseif k <= 15
                    perm[idx] = k_selante_mD * 1e-3 * si_unit(:darcy)
                    poro[idx] = phi_selante
                else
                    perm[idx] = k_FSASP_mD   * 1e-3 * si_unit(:darcy)
                    poro[idx] = phi_FSASP
                end
                if abs(i - i_falha) <= w_falha
                    perm[idx] = k_falha_mD * 1e-3 * si_unit(:darcy)
                    poro[idx] = phi_falha
                end
            end
        end
    end

    return g, reservoir_domain(g,
        permeability = perm,
        porosity     = poro,
        depth        = prof_injecao)
end

# ═══════════════════════════════════════════════
# 4. FUNÇÃO: RODAR UM CENÁRIO
# ═══════════════════════════════════════════════

function rodar_cenario(;
        i_falha         = 20,
        k_falha_mD      = 100.0,
        P_injecao_fator = 1.10,
        anos            = 20,
        caminho_csv     = CAMINHO_CSV_PADRAO,
        prof_injecao    = PROF_INJECAO_OTIMA,
        cenario_geomec  = "CENTRAL",
        usar_capilaridade = true,
        debug_wd        = false)

    nx, ny, nz   = 30, 1, 20
    Lx, Lz       = 3000.0, 200.0
    dz           = Lz / nz

    n_camadas_res = 12
    camadas = carregar_camadas_reservatorio(
        caminho_csv, prof_injecao, n_camadas_res, dz, cenario_geomec)

    # CORREÇÃO: P_ini agora escala com a profundidade real (gradiente hidrostático
    # calibrado no ponto original do campo), em vez de ficar fixo em 21,62 MPa
    # independente de qual profundidade estamos simulando.
    P_ini        = GRADIENTE_PRESSAO_PA_POR_M * prof_injecao
    P_fratura    = 0.0181 * prof_injecao * 1e6

    coesao_Pa   = mean(camadas.coesao)
    atrito_grau = mean(camadas.atrito)
    P_reativ    = P_ini + coesao_Pa / tand(atrito_grau)

    pc_entrada_media = mean(camadas.pc_entrada)
    swirr_media      = mean(camadas.swirr)
    sgr_media        = mean(camadas.sgr)

    rho_agua, rho_CO2 = 1020.0, 650.0
    mu_agua,  mu_CO2  = 0.5e-3, 0.04e-3
    i_injetor    = 3

    k_selante_mD = 0.0001
    k_FSASP_mD   = 100.0
    phi_selante  = 0.02
    phi_FSASP    = 0.20
    phi_falha    = 0.08

    _, domain = build_domain(nx, ny, nz, Lx, Lz, prof_injecao,
        camadas.kx, camadas.phi,
        k_selante_mD, k_FSASP_mD,
        k_falha_mD, phi_selante,
        phi_FSASP, phi_falha, i_falha, 1)

    phases = (LiquidPhase(), VaporPhase())
    sys    = ImmiscibleSystem(phases,
                reference_densities=[rho_agua, rho_CO2])

    Injetor  = setup_vertical_well(domain, i_injetor, 1,
        name=:Injetor,  heel=1, toe=n_camadas_res)
    Produtor = setup_vertical_well(domain, nx-2, 1,
        name=:Produtor, heel=1, toe=n_camadas_res)

    model, parameters = setup_reservoir_model(domain, sys,
        wells=[Injetor, Produtor])

    c = [1e-6, 1e-4] / si_unit(:bar)

    relperm = BrooksCoreyRelativePermeabilities(
        sys, [2.5, 3.5], [swirr_media, sgr_media])

    replace_variables!(model,
        PhaseMassDensities = ConstantCompressibilityDensities(
            p_ref           = P_ini,
            density_ref     = [rho_agua, rho_CO2],
            compressibility = c),
        PhaseRelativePermeability = relperm,
        PhaseViscosities = [mu_agua, mu_CO2]
    )

    # ── PASSO 1: Pressão capilar — CORREÇÃO: pc precisa ser um vetor (uma
    # entrada por par de fases/região), não uma função escalar sozinha.
    pc_ativada = false
    if usar_capilaridade
        try
            pc_fn = s_water -> pc_entrada_media * (max(s_water - swirr_media, 1e-3) /
                                                     max(1.0 - swirr_media, 1e-3))^(-1.0/2.0)
            replace_variables!(model,
                CapillaryPressure = JutulDarcy.SimpleCapillaryPressure([pc_fn]))
            pc_ativada = true
        catch e1
            try
                # fallback: talvez precise do kwarg regions mesmo assim
                replace_variables!(model,
                    CapillaryPressure = JutulDarcy.SimpleCapillaryPressure([pc_fn]; regions = nothing))
                pc_ativada = true
            catch e2
                @warn "Pc não configurada. Tentativa 1: $e1 | Tentativa 2: $e2"
            end
        end
    end

    state0 = setup_reservoir_state(model,
        Pressure    = P_ini,
        Saturations = [1.0, 0.0])

    I_ctrl = InjectorControl(
        BottomHolePressureTarget(P_ini * P_injecao_fator),
        [0.0, 1.0], density=rho_CO2)
    P_ctrl = ProducerControl(
        BottomHolePressureTarget(P_ini * 0.95))

    forces = setup_reservoir_forces(model,
        control=Dict(:Injetor => I_ctrl, :Produtor => P_ctrl))

    dt = fill(365.0 * si_unit(:day), anos)

    wd, states, t = simulate_reservoir(state0, model, dt,
        parameters=parameters, forces=forces,
        info_level=-1)

    if debug_wd
        println("\n[DEBUG wd] wd.wells[:Injetor] chaves: ", keys(wd.wells[:Injetor]))
        println("[DEBUG wd] mass_rate, primeiros 3: ", first(wd.wells[:Injetor][:mass_rate], 3))
    end

    idx_falha_res = [(k-1)*(nx*ny) + i_falha for k in 1:n_camadas_res]
    idx_selante   = (nx*n_camadas_res+1):(nx*15)
    idx_FSASP     = (nx*15+1):(nx*nz)
    idx_falha_todas_camadas = [(k-1)*(nx*ny) + i_falha for k in 1:nz]

    P_falha_historico   = Float64[]
    S_selante_historico = Float64[]
    S_FSASP_historico    = Float64[]
    S_falha_historico    = Float64[]

    for s in states
        Sf_t = s[:Saturations][2, :]
        P_t  = s[:Pressure]
        push!(P_falha_historico,   maximum(P_t[idx_falha_res]))
        push!(S_selante_historico, maximum(Sf_t[idx_selante]))
        push!(S_FSASP_historico,   maximum(Sf_t[idx_FSASP]))
        push!(S_falha_historico,   maximum(Sf_t[idx_falha_todas_camadas]))
    end

    P_falha   = maximum(P_falha_historico)
    S_selante = maximum(S_selante_historico)
    S_FSASP   = maximum(S_FSASP_historico)
    S_falha   = maximum(S_falha_historico)
    ano_do_pico_pressao = argmax(P_falha_historico)

    margem_fratura  = (P_fratura - P_falha) / 1e6
    margem_reativ   = (P_reativ  - P_falha) / 1e6
    dist_m          = (i_falha - i_injetor) * (Lx/nx)

    reativacao = P_falha >= P_reativ

    risco = if S_FSASP > 1e-4
        "ALTO"
    elseif reativacao || S_falha > 0.01 || margem_fratura < 5.0
        "MEDIO"
    else
        "BAIXO"
    end

    massa_CO2_kg = NaN
    vazao_media_CO2_kg_dia = NaN
    try
        vazoes_massa = wd.wells[:Injetor][:mass_rate]
        massa_CO2_kg = sum(vazoes_massa .* dt)
        vazao_media_CO2_kg_dia = mean(vazoes_massa) * 86400.0
    catch e
        @warn "Erro ao extrair massa de CO2: $e"
    end

    return (
        dist_falha_m         = dist_m,
        k_falha_mD           = k_falha_mD,
        fator_pressao        = P_injecao_fator,
        prof_injecao_m       = prof_injecao,
        P_ini_MPa            = round(P_ini/1e6, digits=3),
        pc_ativada           = pc_ativada,
        cenario_geomec       = cenario_geomec,
        coesao_usada_MPa     = round(coesao_Pa/1e6, digits=3),
        atrito_usado_graus   = round(atrito_grau, digits=1),
        pc_entrada_usada_kPa = round(pc_entrada_media/1e3, digits=2),
        swirr_usado          = round(swirr_media, digits=3),
        sgr_usado            = round(sgr_media, digits=3),
        P_falha_MPa          = round(P_falha/1e6,     digits=3),
        ano_pico_pressao     = ano_do_pico_pressao,
        P_reativ_MPa         = round(P_reativ/1e6,    digits=3),
        margem_fratura       = round(margem_fratura,  digits=3),
        margem_reativ        = round(margem_reativ,   digits=3),
        reativacao           = reativacao,
        S_selante            = round(S_selante,       digits=6),
        S_FSASP              = round(S_FSASP,         digits=6),
        S_falha              = round(S_falha,         digits=4),
        massa_CO2_kg         = massa_CO2_kg,
        vazao_media_CO2_kg_dia = vazao_media_CO2_kg_dia,
        risco                = risco
    )
end

# ═══════════════════════════════════════════════
# TESTE ISOLADO
# ═══════════════════════════════════════════════

println("\n" * "═"^60)
println(" TESTE ISOLADO: profundidade típica de CCS + P_ini corrigido + Pc")
println("═"^60)

teste = rodar_cenario(
    i_falha=20, k_falha_mD=100.0, P_injecao_fator=1.10,
    anos=20, cenario_geomec="CENTRAL",
    usar_capilaridade=true, debug_wd=true)

println("\nResultado:")
println("  prof_injecao usada    : ", teste.prof_injecao_m, " m")
println("  P_ini calculado (MPa) : ", teste.P_ini_MPa)
println("  Pc ativada            : ", teste.pc_ativada)
println("  Massa CO2 (kt)        : ", round(teste.massa_CO2_kg/1e6, digits=2))
println("═"^60 * "\n")

# ═══════════════════════════════════════════════
# 5. ANÁLISE DE SENSIBILIDADE
# ═══════════════════════════════════════════════

dist_falhas     = [8, 14, 20]
k_falhas        = [10.0, 100.0, 1000.0, 10000.0]
fat_pressoes    = [1.05, 1.10, 1.20]
cenarios_geomec = ["PESSIMISTA", "CENTRAL", "OTIMISTA"]

total = length(dist_falhas) * length(k_falhas) * length(fat_pressoes) * length(cenarios_geomec)
println("Rodando $total cenários...\n")

resultados = []
cnt = 0

for i_f in dist_falhas
    for k_f in k_falhas
        for fp in fat_pressoes
            for cg in cenarios_geomec
                global cnt += 1
                print("Cenário $cnt/$total  ")
                r = rodar_cenario(
                    i_falha         = i_f,
                    k_falha_mD      = k_f,
                    P_injecao_fator = fp,
                    anos            = 20,
                    cenario_geomec  = cg,
                    caminho_csv     = CAMINHO_CSV_PADRAO,
                    usar_capilaridade = true,
                    debug_wd        = false)
                push!(resultados, r)
                reativ_str = r.reativacao ? " ⚠️ REATIVAÇÃO (ano $(r.ano_pico_pressao))" : ""
                println("dist=$(r.dist_falha_m)m  k=$(k_f)mD  fp=$(fp)  geomec=$cg  " *
                        "pc=$(r.pc_ativada)  massa=$(round(r.massa_CO2_kg/1e6, digits=1))kt  → $(r.risco)$reativ_str")
            end
        end
    end
end

println("\n" * "═"^70)
n_alto  = count(r -> r.risco == "ALTO",  resultados)
n_medio = count(r -> r.risco == "MEDIO", resultados)
n_baixo = count(r -> r.risco == "BAIXO", resultados)
n_reat  = count(r -> r.reativacao,       resultados)
println("RESUMO: ALTO=$n_alto  MÉDIO=$n_medio  BAIXO=$n_baixo  REATIVAÇÕES=$n_reat de $total")
println("═"^70)

df_saida = DataFrame(
    dist_falha_m           = [r.dist_falha_m           for r in resultados],
    k_falha_mD             = [r.k_falha_mD             for r in resultados],
    fator_pressao          = [r.fator_pressao          for r in resultados],
    prof_injecao_m         = [r.prof_injecao_m         for r in resultados],
    P_ini_MPa              = [r.P_ini_MPa              for r in resultados],
    pc_ativada             = [r.pc_ativada             for r in resultados],
    cenario_geomec         = [r.cenario_geomec         for r in resultados],
    coesao_usada_MPa       = [r.coesao_usada_MPa       for r in resultados],
    atrito_usado_graus     = [r.atrito_usado_graus     for r in resultados],
    pc_entrada_usada_kPa   = [r.pc_entrada_usada_kPa   for r in resultados],
    swirr_usado            = [r.swirr_usado            for r in resultados],
    sgr_usado              = [r.sgr_usado              for r in resultados],
    P_falha_MPa            = [r.P_falha_MPa            for r in resultados],
    ano_pico_pressao       = [r.ano_pico_pressao       for r in resultados],
    P_reativ_MPa           = [r.P_reativ_MPa           for r in resultados],
    margem_fratura         = [r.margem_fratura         for r in resultados],
    margem_reativ          = [r.margem_reativ          for r in resultados],
    reativacao             = [r.reativacao             for r in resultados],
    S_selante              = [r.S_selante              for r in resultados],
    S_FSASP                = [r.S_FSASP                for r in resultados],
    S_falha                = [r.S_falha                for r in resultados],
    massa_CO2_kg           = [r.massa_CO2_kg           for r in resultados],
    vazao_media_CO2_kg_dia = [r.vazao_media_CO2_kg_dia for r in resultados],
    risco                  = [r.risco                  for r in resultados]
)

CSV.write(CAMINHO_SAIDA_PADRAO, df_saida)
println("\nTabela salva em: $CAMINHO_SAIDA_PADRAO")