# ═══════════════════════════════════════════════════════════════════
# CAMPO PILAR — FUNÇÕES COMUNS (usado por simulacao_teste1.jl,
# simulacao_agss.jl, simulacao_co2brine.jl e validacao_convergencia.jl)
# ═══════════════════════════════════════════════════════════════════

using JutulDarcy, Jutul, Statistics, CSV, DataFrames

const PASTA_SCRIPT = @__DIR__
const CAMINHO_CSV_PADRAO = joinpath(PASTA_SCRIPT, "malha_1PIR_enriquecida.csv")

if !isfile(CAMINHO_CSV_PADRAO)
    error("Arquivo não encontrado: $CAMINHO_CSV_PADRAO")
end

# Gradiente de pressão calibrado no ponto original do Campo Pilar
const GRADIENTE_PRESSAO_PA_POR_M = 21.62e6 / 2207.7   # ≈ 9793.7 Pa/m

# Geometria vertical fixa do modelo (independente da resolução da malha)
const ESPESSURA_RESERVATORIO_M = 120.0   # 12 blocos de dados de 10 m
const ESPESSURA_SELANTE_M      = 30.0

# ═══════════════════════════════════════════════
# 1. MELHOR JANELA DE RESERVATÓRIO
#    Faixa cobre toda a Fm. Penedo (2250–3980 m no 1PIR). Mínimo de 800 m
#    garante CO2 supercrítico; zona ótima: 3200–3450 m (KX ~80 mD, PHIE ~17%).
# ═══════════════════════════════════════════════

function encontrar_melhor_janela_reservatorio(caminho_csv, n_camadas, dz;
                                                profundidade_min = 800.0,
                                                profundidade_max = 4000.0)
    df = CSV.read(caminho_csv, DataFrame)
    sort!(df, :DEPT_TOPO)
    df_ok = filter(row -> row.COBERTURA_REAL, df)
    df_ok = filter(row -> profundidade_min <= row.DEPT_TOPO <= profundidade_max, df_ok)

    n = nrow(df_ok)
    if n < n_camadas
        error("Poucos blocos válidos ($n) no intervalo $(profundidade_min)–$(profundidade_max) m.")
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

    println("[INFO] Melhor janela de reservatório (faixa $(profundidade_min)–$(profundidade_max) m): " *
            "$(minimum(melhor_janela.DEPT_TOPO))–$(maximum(melhor_janela.DEPT_TOPO))m | " *
            "prof_injecao=$prof_injecao_otima m | KX_médio=$(round(melhor_kx_medio,digits=1))mD | " *
            "arenito=$(sum(melhor_janela.FACIES.==1))/$n_camadas")

    return prof_injecao_otima
end

const PROF_INJECAO_OTIMA = encontrar_melhor_janela_reservatorio(CAMINHO_CSV_PADRAO, 12, 10.0)

# ═══════════════════════════════════════════════
# 2. CARREGAR DADOS REAIS DO POÇO 1PIR (inclui KZ e temperatura)
# ═══════════════════════════════════════════════

function carregar_camadas_reservatorio(caminho_csv, prof_injecao, n_camadas, dz_dados, cenario_geomec)
    df = CSV.read(caminho_csv, DataFrame)
    sort!(df, :DEPT_TOPO)
    df_ok = filter(row -> row.COBERTURA_REAL, df)

    meia_janela = (n_camadas / 2) * dz_dados
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
        kz     = df_zona.KZ,                    # permeabilidade vertical real
        phi    = df_zona.PHIE,
        coesao = df_zona[!, coesao_col],
        atrito = df_zona[!, atrito_col],
        pc_entrada = df_zona[!, pc_col],
        swirr      = df_zona[!, swirr_col],
        sgr        = df_zona[!, sgr_col],
        facies = df_zona.FACIES,
        dept   = df_zona.DEPT_TOPO,
        temperatura_C = df_zona.TEMPERATURA_C,  # temperatura real do poço
    )
end

# ═══════════════════════════════════════════════
# 3. DOMÍNIO — tensor de permeabilidade [kx, kx, kz], regiões de Pc,
#    geometria vertical definida por PROFUNDIDADE (permite refinar a malha
#    sem mudar a espessura física de reservatório/selante/FSASP).
# ═══════════════════════════════════════════════

function build_domain(nx, ny, nz, Lx, Ly, Lz, prof_injecao,
                      kx_res, kz_res, phi_res,
                      k_selante_mD, k_FSASP_mD, k_falha_mD,
                      phi_selante, phi_FSASP, phi_falha,
                      i_falha, w_falha;
                      kvkh_selante = 0.1,
                      kvkh_FSASP   = 0.1,
                      temperatura_C = missing)

    n_dados = length(kx_res)
    dz_dados = ESPESSURA_RESERVATORIO_M / n_dados
    g  = CartesianMesh((nx, ny, nz), (Lx, Ly, Lz))
    nc = nx * ny * nz
    dz = Lz / nz

    perm = zeros(3, nc)          # linhas: kx, ky, kz
    poro = zeros(nc)
    regiao_pc = ones(Int, nc)    # 1 = reservatório/falha/FSASP, 2 = selante

    mD = 1e-3 * si_unit(:darcy)

    for k in 1:nz, j in 1:ny, i in 1:nx
        idx = (k-1)*(nx*ny) + (j-1)*nx + i
        z_c = (k - 0.5) * dz     # profundidade do centro da célula no domínio
        if z_c <= ESPESSURA_RESERVATORIO_M
            b = clamp(ceil(Int, z_c / dz_dados), 1, n_dados)
            kh = kx_res[b] * mD
            kv = kz_res[b] * mD
            poro[idx] = phi_res[b]
        elseif z_c <= ESPESSURA_RESERVATORIO_M + ESPESSURA_SELANTE_M
            kh = k_selante_mD * mD
            kv = kvkh_selante * k_selante_mD * mD
            poro[idx] = phi_selante
            regiao_pc[idx] = 2
        else
            kh = k_FSASP_mD * mD
            kv = kvkh_FSASP * k_FSASP_mD * mD
            poro[idx] = phi_FSASP
        end
        if abs(i - i_falha) <= w_falha
            kh = k_falha_mD * mD
            kv = k_falha_mD * mD     # falha = conduíte vertical
            poro[idx] = phi_falha
            regiao_pc[idx] = 1
        end
        perm[1, idx] = kh
        perm[2, idx] = kh
        perm[3, idx] = kv
    end

    kwargs = Dict{Symbol, Any}(
        :permeability => perm,
        :porosity     => poro,
        :depth        => prof_injecao)
    if !ismissing(temperatura_C)
        kwargs[:temperature] = fill(convert_to_si(temperatura_C, :Celsius), nc)
    end

    return g, reservoir_domain(g; kwargs...), regiao_pc
end

# ═══════════════════════════════════════════════
# 4. PRESSÃO CAPILAR POR REGIÃO (reservatório × selante)
#    Brooks-Corey em função da saturação de gás; o selante recebe pressão
#    de entrada própria (ordem de MPa) — mecanismo real de contenção.
# ═══════════════════════════════════════════════

function tabela_pc(pc_entrada, swirr; fator_max = 50.0)
    sg_tab = collect(range(0.0, 1.0, length = 50))
    pc_tab = map(sg_tab) do sg
        sw = 1.0 - sg
        se = clamp((sw - swirr) / max(1.0 - swirr, 1e-3), 1e-3, 1.0)
        min(pc_entrada * se^(-0.5), fator_max * pc_entrada)
    end
    return Jutul.get_1d_interpolator(sg_tab, pc_tab)
end

function configurar_pc!(model, pc_entrada_res, swirr_res, pc_entrada_selante, regiao_pc)
    t_res = tabela_pc(pc_entrada_res, swirr_res)
    t_sel = tabela_pc(pc_entrada_selante, 0.05; fator_max = 5.0)
    rmodel = reservoir_model(model)
    rmodel.secondary_variables[:CapillaryPressure] =
        JutulDarcy.SimpleCapillaryPressure(([t_res, t_sel],); regions = regiao_pc)
    return true
end

# ═══════════════════════════════════════════════
# 5. VISCOSIDADES CONSTANTES (parâmetro na API atual do JutulDarcy)
# ═══════════════════════════════════════════════

function ajustar_viscosidades!(parameters, mu_liq, mu_gas)
    for (mkey, prm) in parameters
        if haskey(prm, :PhaseViscosities)
            mu = prm[:PhaseViscosities]
            mu[1, :] .= mu_liq
            mu[2, :] .= mu_gas
        end
    end
end

# ═══════════════════════════════════════════════
# 6. GEOMECÂNICA — SLIP TENDENCY (tensor de tensões simplificado)
#    σv integrado da densidade bulk dos dados do poço; σh = K0·σv_eff;
#    resolução de tensões no plano da falha com mergulho `dip`.
#    Retorna também a pressão crítica de reativação no plano da falha.
# ═══════════════════════════════════════════════

function analise_slip_tendency(caminho_csv, prof;
                                dip_graus = 60.0,
                                K0 = 0.7,
                                coesao_Pa = 5.0e6,
                                atrito_grau = 30.0)
    df = CSV.read(caminho_csv, DataFrame)
    sort!(df, :DEPT_TOPO)
    grav = 9.81

    # σv: integra densidade bulk em passos de 10 m usando o bloco de dados
    # mais próximo (lacunas preenchidas com o vizinho mais próximo)
    sigma_v = 0.0
    zs = 5.0:10.0:prof
    for z in zs
        i_prox = argmin(abs.(df.DEPT_TOPO .- z))
        phi = clamp(df.PHIE[i_prox], 0.0, 0.35)
        rho_bulk = (1.0 - phi) * 2650.0 + phi * 1030.0
        sigma_v += rho_bulk * grav * 10.0
    end

    p0 = GRADIENTE_PRESSAO_PA_POR_M * prof
    sv_eff = sigma_v - p0
    sh_eff = K0 * sv_eff

    dip = deg2rad(dip_graus)
    # σ1 = σv (falha normal); plano com mergulho dip:
    sn_eff = 0.5*(sv_eff + sh_eff) + 0.5*(sv_eff - sh_eff)*cos(2*dip)
    tau    = 0.5*(sv_eff - sh_eff)*sin(2*dip)

    ST = tau / max(sn_eff, 1.0)              # slip tendency
    mu_f = tand(atrito_grau)
    # Reativação quando τ ≥ c + μ(σn_eff − ΔP)  →  ΔP_crit:
    dP_crit = sn_eff - (tau - coesao_Pa) / mu_f
    P_reativ_slip = p0 + max(dP_crit, 0.0)

    return (
        sigma_v_MPa       = sigma_v / 1e6,
        sn_eff_MPa        = sn_eff / 1e6,
        tau_MPa           = tau / 1e6,
        slip_tendency     = ST,
        P_reativ_slip_MPa = P_reativ_slip / 1e6,
    )
end

# ═══════════════════════════════════════════════
# 7. CENÁRIO CCS IMISCÍVEL (usado por simulacao_teste1.jl e
#    validacao_convergencia.jl)
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
        # ── malha (o refinamento preserva a geometria física) ──
        nx = 30, ny = 1, nz = 20,
        Lx = 3000.0, Lz = 200.0,
        i_injetor = max(3, round(Int, nx/10)),
        # ── propriedades de selante/FSASP (sensibilidade da atividade 8) ──
        k_selante_mD   = 0.0001,
        pc_selante_MPa = 2.0,
        k_FSASP_mD     = 100.0,
        # ── geomecânica slip tendency ──
        dip_falha_graus = 60.0,
        K0_tensoes      = 0.7,
        debug_wd = false)

    Ly = ny == 1 ? 1.0 : 1000.0
    dz = Lz / nz
    n_res_grid = round(Int, ESPESSURA_RESERVATORIO_M / dz)   # camadas de malha no reservatório

    camadas = carregar_camadas_reservatorio(
        caminho_csv, prof_injecao, 12, 10.0, cenario_geomec)

    P_ini     = GRADIENTE_PRESSAO_PA_POR_M * prof_injecao
    P_fratura = 0.0181 * prof_injecao * 1e6

    coesao_Pa   = mean(camadas.coesao)
    atrito_grau = mean(camadas.atrito)
    P_reativ    = P_ini + coesao_Pa / tand(atrito_grau)   # critério simplificado (mantido)

    geomec = analise_slip_tendency(caminho_csv, prof_injecao;
        dip_graus = dip_falha_graus, K0 = K0_tensoes,
        coesao_Pa = coesao_Pa, atrito_grau = atrito_grau)

    pc_entrada_media = mean(camadas.pc_entrada)
    swirr_media      = mean(camadas.swirr)
    sgr_media        = mean(camadas.sgr)

    rho_agua, rho_CO2 = 1020.0, 650.0
    mu_agua,  mu_CO2  = 0.5e-3, 0.04e-3

    phi_selante  = 0.02
    phi_FSASP    = 0.20
    phi_falha    = 0.08

    _, domain, regiao_pc = build_domain(nx, ny, nz, Lx, Ly, Lz, prof_injecao,
        camadas.kx, camadas.kz, camadas.phi,
        k_selante_mD, k_FSASP_mD, k_falha_mD,
        phi_selante, phi_FSASP, phi_falha, i_falha, 1)

    phases = (LiquidPhase(), VaporPhase())
    sys    = ImmiscibleSystem(phases, reference_densities=[rho_agua, rho_CO2])

    Injetor  = setup_vertical_well(domain, i_injetor, 1,
        name=:Injetor,  heel=1, toe=n_res_grid)
    Produtor = setup_vertical_well(domain, nx-2, 1,
        name=:Produtor, heel=1, toe=n_res_grid)

    model = setup_reservoir_model(domain, sys, wells=[Injetor, Produtor])

    c = [1e-6, 1e-4] / si_unit(:bar)
    relperm = BrooksCoreyRelativePermeabilities(
        sys, [2.5, 3.5], [swirr_media, sgr_media])

    replace_variables!(model,
        PhaseMassDensities = ConstantCompressibilityDensities(
            p_ref           = P_ini,
            density_ref     = [rho_agua, rho_CO2],
            compressibility = c),
        RelativePermeabilities = relperm
    )

    pc_ativada = false
    if usar_capilaridade
        try
            configurar_pc!(model, pc_entrada_media, swirr_media,
                           pc_selante_MPa * 1e6, regiao_pc)
            pc_ativada = true
        catch e1
            @warn "Pc não configurada: $e1"
        end
    end

    parameters = setup_parameters(model)
    ajustar_viscosidades!(parameters, mu_agua, mu_CO2)

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
        parameters=parameters, forces=forces, info_level=-1)

    if debug_wd
        println("\n[DEBUG wd] wd.wells[:Injetor] chaves: ", keys(wd.wells[:Injetor]))
    end

    # índices por profundidade física (robustos ao refinamento da malha)
    z_c(k) = (k - 0.5) * dz
    ks_res    = [k for k in 1:nz if z_c(k) <= ESPESSURA_RESERVATORIO_M]
    ks_sel    = [k for k in 1:nz if ESPESSURA_RESERVATORIO_M < z_c(k) <= ESPESSURA_RESERVATORIO_M + ESPESSURA_SELANTE_M]
    ks_fsasp  = [k for k in 1:nz if z_c(k) > ESPESSURA_RESERVATORIO_M + ESPESSURA_SELANTE_M]
    idx_falha_res  = [(k-1)*(nx*ny) + i_falha for k in ks_res]
    # selante INTACTO: exclui as colunas da falha (a falha que corta o selo
    # é vazamento pela falha, não falha de contenção do capeador)
    eh_falha(i) = abs(i - i_falha) <= 1
    idx_selante    = [(k-1)*(nx*ny) + (j-1)*nx + i for k in ks_sel   for j in 1:ny for i in 1:nx if !eh_falha(i)]
    idx_falha_selo = [(k-1)*(nx*ny) + (j-1)*nx + i for k in ks_sel   for j in 1:ny for i in 1:nx if  eh_falha(i)]
    idx_FSASP      = vcat([collect(((k-1)*nx*ny + 1):(k*nx*ny)) for k in ks_fsasp]...)
    idx_falha_toda = [(k-1)*(nx*ny) + i_falha for k in 1:nz]

    P_falha_hist     = Float64[]
    S_selante_hist   = Float64[]
    S_falha_selo_hist = Float64[]
    S_FSASP_hist     = Float64[]
    S_falha_hist     = Float64[]
    for s in states
        Sf_t = s[:Saturations][2, :]
        P_t  = s[:Pressure]
        push!(P_falha_hist,    maximum(P_t[idx_falha_res]))
        push!(S_selante_hist,  maximum(Sf_t[idx_selante]))
        push!(S_falha_selo_hist, maximum(Sf_t[idx_falha_selo]))
        push!(S_FSASP_hist,    maximum(Sf_t[idx_FSASP]))
        push!(S_falha_hist,    maximum(Sf_t[idx_falha_toda]))
    end

    P_falha      = maximum(P_falha_hist)
    S_selante    = maximum(S_selante_hist)
    S_falha_selo = maximum(S_falha_selo_hist)
    S_FSASP      = maximum(S_FSASP_hist)
    S_falha      = maximum(S_falha_hist)
    ano_do_pico_pressao = argmax(P_falha_hist)

    margem_fratura = (P_fratura - P_falha) / 1e6
    margem_reativ  = (P_reativ  - P_falha) / 1e6
    dist_m         = (i_falha - i_injetor) * (Lx/nx)

    reativacao      = P_falha >= P_reativ
    reativacao_slip = P_falha >= geomec.P_reativ_slip_MPa * 1e6

    risco = if S_FSASP > 1e-4 || S_selante > 1e-3   # migração acima do selo OU capeador intacto invadido
        "ALTO"
    elseif reativacao || reativacao_slip || S_falha > 0.01 || margem_fratura < 5.0
        "MEDIO"
    else
        "BAIXO"
    end

    massa_CO2_kg = NaN
    vazao_media_CO2_kg_dia = NaN
    injetividade_kg_dia_MPa = NaN
    try
        vazoes_massa = wd.wells[:Injetor][:mass_rate]
        massa_CO2_kg = sum(vazoes_massa .* dt)
        vazao_media_CO2_kg_dia = mean(vazoes_massa) * 86400.0
        bhp = wd.wells[:Injetor][:bhp]
        dP_medio = mean(bhp) - P_ini
        if dP_medio > 0
            injetividade_kg_dia_MPa = vazao_media_CO2_kg_dia / (dP_medio / 1e6)
        end
    catch e
        @warn "Erro ao extrair dados do poço: $e"
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
        k_selante_mD         = k_selante_mD,
        pc_selante_MPa       = pc_selante_MPa,
        k_FSASP_mD           = k_FSASP_mD,
        nx = nx, ny = ny, nz = nz,
        P_falha_MPa          = round(P_falha/1e6,     digits=3),
        ano_pico_pressao     = ano_do_pico_pressao,
        P_reativ_MPa         = round(P_reativ/1e6,    digits=3),
        margem_fratura       = round(margem_fratura,  digits=3),
        margem_reativ        = round(margem_reativ,   digits=3),
        reativacao           = reativacao,
        slip_tendency        = round(geomec.slip_tendency, digits=3),
        P_reativ_slip_MPa    = round(geomec.P_reativ_slip_MPa, digits=3),
        reativacao_slip      = reativacao_slip,
        S_selante            = round(S_selante,       digits=6),
        S_falha_selo         = round(S_falha_selo,    digits=4),
        S_FSASP              = round(S_FSASP,         digits=6),
        S_falha              = round(S_falha,         digits=4),
        massa_CO2_kg         = massa_CO2_kg,
        vazao_media_CO2_kg_dia  = vazao_media_CO2_kg_dia,
        injetividade_kg_dia_MPa = injetividade_kg_dia_MPa,
        risco                = risco
    )
end

println("[INFO] campo_pilar_common.jl carregado.")
