# ═══════════════════════════════════════════════════════════════════
# CAMPO PILAR — AGSS: ARMAZENAMENTO GEOLÓGICO DE GÁS NATURAL
# Ciclos sazonais de injeção/retirada com HISTERESE de permeabilidade
# relativa (Killough) — atividades 5 e 6 do cronograma.
# Poço único operando como injetor (6 meses) e produtor (6 meses).
# ═══════════════════════════════════════════════════════════════════

include(joinpath(@__DIR__, "campo_pilar_common.jl"))

import JutulDarcy: KilloughHysteresis, ReservoirRelativePermeabilities,
                   add_relperm_parameters!, brooks_corey_relperm

const CAMINHO_SAIDA_AGSS = joinpath(PASTA_SCRIPT, "agss_campo_pilar.csv")

function rodar_cenario_agss(;
        anos            = 10,
        fp_injecao      = 1.10,   # BHP injeção = fp × P_ini
        fp_retirada     = 0.90,   # BHP retirada = fp × P_ini
        i_falha         = 20,
        k_falha_mD      = 100.0,
        cenario_geomec  = "CENTRAL",
        prof_injecao    = PROF_INJECAO_OTIMA,
        caminho_csv     = CAMINHO_CSV_PADRAO,
        nx = 30, ny = 1, nz = 20,
        Lx = 3000.0, Lz = 200.0,
        k_selante_mD   = 0.0001,
        pc_selante_MPa = 2.0,
        k_FSASP_mD     = 100.0)

    Ly = ny == 1 ? 1.0 : 1000.0
    dz = Lz / nz
    n_res_grid = round(Int, ESPESSURA_RESERVATORIO_M / dz)
    i_injetor  = max(3, round(Int, nx/10))

    camadas = carregar_camadas_reservatorio(caminho_csv, prof_injecao, 12, 10.0, cenario_geomec)

    P_ini     = GRADIENTE_PRESSAO_PA_POR_M * prof_injecao
    P_fratura = 0.0181 * prof_injecao * 1e6
    coesao_Pa   = mean(camadas.coesao)
    atrito_grau = mean(camadas.atrito)
    P_reativ    = P_ini + coesao_Pa / tand(atrito_grau)

    swirr = mean(camadas.swirr)
    sgr   = mean(camadas.sgr)
    pc_entrada = mean(camadas.pc_entrada)

    # Gás natural (metano) nas condições do reservatório (~33 MPa, ~87 °C):
    rho_agua, rho_gas = 1020.0, 190.0
    mu_agua,  mu_gas  = 0.5e-3, 0.02e-3

    _, domain, regiao_pc = build_domain(nx, ny, nz, Lx, Ly, Lz, prof_injecao,
        camadas.kx, camadas.kz, camadas.phi,
        k_selante_mD, k_FSASP_mD, k_falha_mD,
        0.02, 0.20, 0.08, i_falha, 1)

    sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()),
        reference_densities=[rho_agua, rho_gas])

    Poco = setup_vertical_well(domain, i_injetor, 1,
        name=:Poco, heel=1, toe=n_res_grid)

    model = setup_reservoir_model(domain, sys, wells=[Poco])

    # Densidades: gás muito mais compressível (essencial para estocagem)
    c = [1e-6, 3e-3] / si_unit(:bar)
    replace_variables!(model,
        PhaseMassDensities = ConstantCompressibilityDensities(
            p_ref = P_ini, density_ref = [rho_agua, rho_gas],
            compressibility = c))

    # ── Permeabilidade relativa com HISTERESE (Killough) ──
    # Drenagem (injeção): gás residual baixo. Embebição (retirada): gás
    # aprisionado = SGR do CSV → controla o "cushion gas" irrecuperável.
    sw_tab = collect(range(0.0, 1.0, length = 50))
    sg_tab = collect(range(0.0, 1.0, length = 50))
    krw   = PhaseRelativePermeability(sw_tab,
              brooks_corey_relperm.(sw_tab, n = 2.5, residual = swirr), label = :w)
    krg_d = PhaseRelativePermeability(sg_tab,
              brooks_corey_relperm.(sg_tab, n = 2.0, residual = 0.05), label = :g)
    krg_i = PhaseRelativePermeability(sg_tab,
              brooks_corey_relperm.(sg_tab, n = 3.0, residual = sgr), label = :g)

    relperm = ReservoirRelativePermeabilities(
        w = krw, g = (krg_d, krg_i),
        hysteresis_g = KilloughHysteresis())
    replace_variables!(model, RelativePermeabilities = relperm)
    add_relperm_parameters!(model)

    configurar_pc!(model, pc_entrada, swirr, pc_selante_MPa * 1e6, regiao_pc)

    parameters = setup_parameters(model)
    ajustar_viscosidades!(parameters, mu_agua, mu_gas)

    state0 = setup_reservoir_state(model,
        Pressure = P_ini, Saturations = [1.0, 0.0])

    I_ctrl = InjectorControl(
        BottomHolePressureTarget(P_ini * fp_injecao),
        [0.0, 1.0], density = rho_gas)
    W_ctrl = ProducerControl(
        BottomHolePressureTarget(P_ini * fp_retirada))

    # Aquífero lateral: pressão constante na borda direita do reservatório.
    # Sem isso o domínio fechado equilibra com o BHP e a injeção cessa.
    zc_bc(k) = (k - 0.5) * dz
    ks_res_bc  = [k for k in 1:nz if zc_bc(k) <= ESPESSURA_RESERVATORIO_M]
    celulas_bc = [(k-1)*(nx*ny) + nx for k in ks_res_bc]
    bc = flow_boundary_condition(celulas_bc, domain, P_ini, fractional_flow = [1.0, 0.0])

    forces_inj = setup_reservoir_forces(model, control = Dict(:Poco => I_ctrl), bc = bc)
    forces_ret = setup_reservoir_forces(model, control = Dict(:Poco => W_ctrl), bc = bc)

    # Ciclo anual: 6 meses injeção + 6 meses retirada, passos mensais
    dia = si_unit(:day)
    dt_mes = 30.4375 * dia
    ciclo_forces = vcat(fill(forces_inj, 6), fill(forces_ret, 6))
    forces = repeat(ciclo_forces, anos)
    dt = fill(dt_mes, 12 * anos)

    wd, states, t = simulate_reservoir(state0, model, dt,
        parameters = parameters, forces = forces, info_level = -1)

    # ── Métricas por ciclo (SÓ GÁS: :mass_rate incluiria a água produzida) ──
    q_gas  = wd.wells[:Poco][:Vapor_mass_rate]    # kg/s de gás
    q_agua = wd.wells[:Poco][:Liquid_mass_rate]   # kg/s de água
    m_gas  = q_gas  .* dt
    m_agua = q_agua .* dt

    resultados_ciclos = DataFrame(ciclo = Int[], injetado_kt = Float64[],
        retirado_kt = Float64[], eficiencia = Float64[], agua_produzida_kt = Float64[])
    for cicl in 1:anos
        passos = ((cicl-1)*12 + 1):(cicl*12)
        inj = sum(max.(m_gas[passos], 0.0)) / 1e6
        ret = -sum(min.(m_gas[passos], 0.0)) / 1e6
        agua = -sum(min.(m_agua[passos], 0.0)) / 1e6
        push!(resultados_ciclos, (cicl, round(inj, digits=3), round(ret, digits=3),
              round(inj > 0 ? ret/inj : NaN, digits=3), round(agua, digits=3)))
    end

    # ── Contenção e geomecânica ──
    z_c(k) = (k - 0.5) * dz
    ks_res   = [k for k in 1:nz if z_c(k) <= ESPESSURA_RESERVATORIO_M]
    ks_sel   = [k for k in 1:nz if ESPESSURA_RESERVATORIO_M < z_c(k) <= ESPESSURA_RESERVATORIO_M + ESPESSURA_SELANTE_M]
    ks_fsasp = [k for k in 1:nz if z_c(k) > ESPESSURA_RESERVATORIO_M + ESPESSURA_SELANTE_M]
    eh_falha(i) = abs(i - i_falha) <= 1
    idx_falha_res = [(k-1)*(nx*ny) + i_falha for k in ks_res]
    idx_selante   = [(k-1)*(nx*ny) + i for k in ks_sel for i in 1:nx if !eh_falha(i)]
    idx_FSASP     = vcat([collect(((k-1)*nx*ny + 1):(k*nx*ny)) for k in ks_fsasp]...)

    P_falha = -Inf; S_selante = 0.0; S_FSASP = 0.0
    for s in states
        Sf = s[:Saturations][2, :]
        P_falha   = max(P_falha,   maximum(s[:Pressure][idx_falha_res]))
        S_selante = max(S_selante, maximum(Sf[idx_selante]))
        S_FSASP   = max(S_FSASP,   maximum(Sf[idx_FSASP]))
    end

    geomec = analise_slip_tendency(caminho_csv, prof_injecao;
        coesao_Pa = coesao_Pa, atrito_grau = atrito_grau)

    return (ciclos = resultados_ciclos,
            P_falha_MPa = round(P_falha/1e6, digits=3),
            P_reativ_MPa = round(P_reativ/1e6, digits=3),
            P_reativ_slip_MPa = geomec.P_reativ_slip_MPa,
            reativacao = P_falha >= min(P_reativ, geomec.P_reativ_slip_MPa*1e6),
            margem_fratura_MPa = round((P_fratura - P_falha)/1e6, digits=2),
            S_selante = S_selante, S_FSASP = S_FSASP,
            fp_injecao = fp_injecao, fp_retirada = fp_retirada)
end

# ═══════════════════════════════════════════════
# EXECUÇÃO: caso base + sensibilidade de pressões operacionais
# ═══════════════════════════════════════════════

println("\n" * "═"^60)
println(" AGSS — ciclos sazonais de gás natural com histerese Killough")
println("═"^60)

todas_linhas = DataFrame()
for (fpi, fpr) in [(1.10, 0.90), (1.05, 0.90), (1.10, 0.85), (1.05, 0.95)]
    println("\n─── fp_inj=$fpi | fp_ret=$fpr ───")
    r = rodar_cenario_agss(fp_injecao = fpi, fp_retirada = fpr, anos = 10)
    ef_medio = mean(skipmissing(r.ciclos.eficiencia[2:end]))  # 1º ciclo enche cushion gas
    println("  Working gas médio (ciclos 2+): injetado=$(round(mean(r.ciclos.injetado_kt[2:end]),digits=1))kt  " *
            "retirado=$(round(mean(r.ciclos.retirado_kt[2:end]),digits=1))kt  eficiência=$(round(ef_medio,digits=3))")
    println("  P_falha=$(r.P_falha_MPa)MPa  P_reativ=$(r.P_reativ_MPa)/$(round(r.P_reativ_slip_MPa,digits=1))MPa  " *
            "reativação=$(r.reativacao)  S_selante=$(round(r.S_selante,digits=5))  S_FSASP=$(round(r.S_FSASP,digits=5))")
    ciclos = copy(r.ciclos)
    ciclos[!, :fp_injecao]  .= fpi
    ciclos[!, :fp_retirada] .= fpr
    ciclos[!, :P_falha_MPa] .= r.P_falha_MPa
    ciclos[!, :reativacao]  .= r.reativacao
    global todas_linhas = vcat(todas_linhas, ciclos)
end

CSV.write(CAMINHO_SAIDA_AGSS, todas_linhas)
println("\nTabela por ciclo salva em: $CAMINHO_SAIDA_AGSS")
