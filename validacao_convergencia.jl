# ═══════════════════════════════════════════════════════════════════
# CAMPO PILAR — VALIDAÇÃO POR CONVERGÊNCIA DE MALHA (atividade 9)
# Refina a malha mantendo a geometria física (reservatório 120 m,
# selante 30 m, posições físicas de poço e falha) e verifica a
# estabilidade de P_falha, massa injetada e saturações.
# ═══════════════════════════════════════════════════════════════════

include(joinpath(@__DIR__, "campo_pilar_common.jl"))

const CAMINHO_SAIDA_VAL = joinpath(PASTA_SCRIPT, "validacao_convergencia.csv")

# posição física alvo: injetor a ~10% de Lx, falha a ~2000 m (célula 20/30)
malhas = [(30, 20), (60, 40), (90, 60)]

println("\n" * "═"^60)
println(" CONVERGÊNCIA DE MALHA — caso central (fp=1.10, k_falha=100 mD)")
println(" AVISO: malhas finas podem demorar bastante.")
println("═"^60)

resultados = []
for (nx, nz) in malhas
    i_falha   = round(Int, 20/30 * nx)
    println("\n─── Malha $(nx)×$(nz) (i_falha=$i_falha) ───")
    tempo = @elapsed r = rodar_cenario(
        nx = nx, nz = nz,
        i_falha = i_falha,
        k_falha_mD = 100.0, P_injecao_fator = 1.10,
        anos = 20, cenario_geomec = "CENTRAL")
    push!(resultados, (nx = nx, nz = nz, dx_m = 3000.0/nx, dz_m = 200.0/nz,
        P_falha_MPa = r.P_falha_MPa,
        massa_CO2_kt = round(r.massa_CO2_kg/1e6, digits=2),
        S_falha = r.S_falha, S_selante = r.S_selante, S_FSASP = r.S_FSASP,
        injetividade_kg_dia_MPa = round(r.injetividade_kg_dia_MPa, digits=1),
        risco = r.risco, tempo_s = round(tempo, digits=1)))
    println("  P_falha=$(r.P_falha_MPa)MPa  massa=$(round(r.massa_CO2_kg/1e6,digits=1))kt  " *
            "S_falha=$(r.S_falha)  [$(round(tempo,digits=0))s]")
end

df = DataFrame(resultados)
println("\n", df)

# variação relativa entre malha mais grossa e mais fina
vp = abs(df.P_falha_MPa[end] - df.P_falha_MPa[1]) / df.P_falha_MPa[1] * 100
vm = abs(df.massa_CO2_kt[end] - df.massa_CO2_kt[1]) / df.massa_CO2_kt[1] * 100
println("\nVariação grossa→fina: P_falha=$(round(vp,digits=1))%  massa=$(round(vm,digits=1))%")
println(vp < 5 && vm < 10 ?
    "✓ Resultados aproximadamente convergidos (P<5%, massa<10%)." :
    "⚠️ Variação ainda significativa — considere refinar mais ou reduzir dt.")

CSV.write(CAMINHO_SAIDA_VAL, df)
println("\nTabela salva em: $CAMINHO_SAIDA_VAL")
