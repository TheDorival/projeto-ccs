# ═══════════════════════════════════════════════════════════════════
# CAMPO PILAR — CCS COM PVT REAL (sistema :co2brine do JutulDarcy)
# Atividade 2 do cronograma: caracterização de fluidos com PVT real.
# Inclui: solubilidade mútua CO2-brine (trapeamento por dissolução),
# propriedades dependentes de P e T (temperatura real do poço 1PIR),
# histerese de kr do gás (trapeamento residual) e período pós-injeção.
# ═══════════════════════════════════════════════════════════════════

include(joinpath(@__DIR__, "campo_pilar_common.jl"))

import JutulDarcy: KilloughHysteresis, ReservoirRelativePermeabilities,
                   add_relperm_parameters!, brooks_corey_relperm

const CAMINHO_SAIDA_CO2B = joinpath(PASTA_SCRIPT, "co2brine_campo_pilar.csv")

function rodar_cenario_co2brine(;
        anos_injecao    = 20,
        anos_pos        = 30,     # migração/dissolução após parar de injetar
        fp_injecao      = 1.10,
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
    T_res_C = mean(camadas.temperatura_C)   # temperatura real (~87 °C)

    _, domain, regiao_pc = build_domain(nx, ny, nz, Lx, Ly, Lz, prof_injecao,
        camadas.kx, camadas.kz, camadas.phi,
        k_selante_mD, k_FSASP_mD, k_falha_mD,
        0.02, 0.20, 0.08, i_falha, 1;
        temperatura_C = T_res_C)

    Injetor  = setup_vertical_well(domain, i_injetor, 1,
        name=:Injetor, heel=1, toe=n_res_grid)
    Produtor = setup_vertical_well(domain, nx-2, 1,
        name=:Produtor, heel=1, toe=n_res_grid)

    # ── Sistema CO2-brine com PVT real (k-value: dissolução mútua) ──
    model = setup_reservoir_model(domain, :co2brine,
        wells = [Injetor, Produtor],
        co2_physics = :kvalue)

    # ── kr com histerese do gás (trapeamento residual, SGR do CSV) ──
    so_tab = collect(range(0.0, 1.0, length = 50))
    sg_tab = collect(range(0.0, 1.0, length = 50))
    krog  = PhaseRelativePermeability(so_tab,
              brooks_corey_relperm.(so_tab, n = 2.5, residual = swirr), label = :og)
    krg_d = PhaseRelativePermeability(sg_tab,
              brooks_corey_relperm.(sg_tab, n = 2.0, residual = 0.05), label = :g)
    krg_i = PhaseRelativePermeability(sg_tab,
              brooks_corey_relperm.(sg_tab, n = 3.0, residual = sgr), label = :g)
    relperm = ReservoirRelativePermeabilities(
        og = krog, g = (krg_d, krg_i),
        hysteresis_g = KilloughHysteresis())
    replace_variables!(model, RelativePermeabilities = relperm)
    add_relperm_parameters!(model)

    configurar_pc!(model, pc_entrada, swirr, pc_selante_MPa * 1e6, regiao_pc)

    parameters = setup_parameters(model)

    # Compositional: fração molar global no lugar de saturações
    state0 = setup_reservoir_state(model,
        Pressure = P_ini,
        OverallMoleFractions = [1.0, 0.0])

    rho_CO2_ref = 650.0
    I_ctrl = InjectorControl(
        BottomHolePressureTarget(P_ini * fp_injecao),
        [0.0, 1.0], density = rho_CO2_ref)
    P_ctrl = ProducerControl(BottomHolePressureTarget(P_ini * 0.95))

    forces_inj = setup_reservoir_forces(model,
        control = Dict(:Injetor => I_ctrl, :Produtor => P_ctrl))
    forces_pos = setup_reservoir_forces(model,
        control = Dict(:Injetor => DisabledControl(), :Produtor => DisabledControl()))

    dia = si_unit(:day)
    dt_inj = fill(365.0 * dia, anos_injecao)
    dt_pos = fill(5 * 365.0 * dia, ceil(Int, anos_pos / 5))
    dt = vcat(dt_inj, dt_pos)
    forces = vcat(fill(forces_inj, length(dt_inj)), fill(forces_pos, length(dt_pos)))

    wd, states, t = simulate_reservoir(state0, model, dt,
        parameters = parameters, forces = forces, info_level = -1)

    # ── Métricas: pressão na falha, contenção, partição livre/dissolvido ──
    pv = pore_volume(model, parameters)

    z_c(k) = (k - 0.5) * dz
    ks_res   = [k for k in 1:nz if z_c(k) <= ESPESSURA_RESERVATORIO_M]
    ks_sel   = [k for k in 1:nz if ESPESSURA_RESERVATORIO_M < z_c(k) <= ESPESSURA_RESERVATORIO_M + ESPESSURA_SELANTE_M]
    ks_fsasp = [k for k in 1:nz if z_c(k) > ESPESSURA_RESERVATORIO_M + ESPESSURA_SELANTE_M]
    eh_falha(i) = abs(i - i_falha) <= 1
    idx_falha_res  = [(k-1)*(nx*ny) + i_falha for k in ks_res]
    idx_selante    = [(k-1)*(nx*ny) + i for k in ks_sel for i in 1:nx if !eh_falha(i)]
    idx_falha_selo = [(k-1)*(nx*ny) + i for k in ks_sel for i in 1:nx if  eh_falha(i)]
    idx_FSASP      = vcat([collect(((k-1)*nx*ny + 1):(k*nx*ny)) for k in ks_fsasp]...)

    hist = DataFrame(tempo_anos = Float64[], P_falha_MPa = Float64[],
        S_selante = Float64[], S_falha_selo = Float64[], S_FSASP = Float64[],
        massa_CO2_livre_kt = Float64[], massa_CO2_dissolvida_kt = Float64[],
        fracao_dissolvida = Float64[])

    tempo_acum = 0.0
    for (istep, s) in enumerate(states)
        tempo_acum += dt[istep] / (365.0 * dia)
        Sg  = s[:Saturations][2, :]
        Sl  = s[:Saturations][1, :]
        P   = s[:Pressure]
        rhoL = s[:PhaseMassDensities][1, :]
        rhoV = s[:PhaseMassDensities][2, :]
        Xco2 = s[:LiquidMassFractions][2, :]   # fração mássica de CO2 na brine

        m_livre = sum(rhoV .* Sg .* pv) / 1e6          # kt
        m_diss  = sum(Xco2 .* rhoL .* Sl .* pv) / 1e6  # kt
        frac = (m_livre + m_diss) > 0 ? m_diss / (m_livre + m_diss) : 0.0

        push!(hist, (round(tempo_acum, digits=1),
            round(maximum(P[idx_falha_res])/1e6, digits=3),
            round(maximum(Sg[idx_selante]), digits=6),
            round(maximum(Sg[idx_falha_selo]), digits=4),
            round(maximum(Sg[idx_FSASP]), digits=6),
            round(m_livre, digits=2), round(m_diss, digits=2),
            round(frac, digits=4)))
    end

    geomec = analise_slip_tendency(caminho_csv, prof_injecao;
        coesao_Pa = coesao_Pa, atrito_grau = atrito_grau)

    P_falha_max = maximum(hist.P_falha_MPa) * 1e6
    return (historico = hist,
            T_res_C = T_res_C,
            P_ini_MPa = round(P_ini/1e6, digits=2),
            P_reativ_MPa = round(P_reativ/1e6, digits=2),
            P_reativ_slip_MPa = round(geomec.P_reativ_slip_MPa, digits=2),
            reativacao = P_falha_max >= min(P_reativ, geomec.P_reativ_slip_MPa*1e6),
            margem_fratura_MPa = round((P_fratura - P_falha_max)/1e6, digits=2))
end

# ═══════════════════════════════════════════════
# EXECUÇÃO
# ═══════════════════════════════════════════════

println("\n" * "═"^60)
println(" CCS com PVT real (:co2brine, k-value) — 20 anos injeção + 30 pós")
println("═"^60)

r = rodar_cenario_co2brine()

println("\nT reservatório: $(round(r.T_res_C, digits=1)) °C | P_ini: $(r.P_ini_MPa) MPa")
println("Reativação: $(r.reativacao) | P_reativ=$(r.P_reativ_MPa)/$(r.P_reativ_slip_MPa) MPa | " *
        "margem fratura=$(r.margem_fratura_MPa) MPa")
println("\nEvolução (últimas 5 linhas):")
show(last(r.historico, 5), allrows=true, allcols=true)
println("\n\nFração dissolvida final: $(last(r.historico.fracao_dissolvida)*100) %")

CSV.write(CAMINHO_SAIDA_CO2B, r.historico)
println("Histórico salvo em: $CAMINHO_SAIDA_CO2B")
