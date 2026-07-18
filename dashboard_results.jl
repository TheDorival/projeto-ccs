using CSV, DataFrames, Statistics

df = CSV.read("sensibilidade_campo_pilar.csv", DataFrame)

P_fratura = 0.0181 * 2207.7
P_reativacao = 21.62 + 5.0 / tand(30.0)

println("═"^80)
println("  DASHBOARD — CAMPO PILAR CCS — BACIA SEAL")
println("═"^80)
println("  Dados carregados: $(nrow(df)) cenários\n")

n_alto = count(==("ALTO"), df.risco)
n_medio = count(==("MEDIO"), df.risco)
n_baixo = count(==("BAIXO"), df.risco)
n_reat = count(df.reativacao)

println("┌─────────────────────────────────────────────────────────────────────────────┐")
println("│                           📊 RESUMO GERAL                                   │")
println("├─────────────────────────────────────────────────────────────────────────────┤")
println("│  Cenários: $(nrow(df))  │  🔴 ALTO: $n_alto  │  🟡 MÉDIO: $n_medio  │  🟢 BAIXO: $n_baixo  │")
println("│  Reativações: $n_reat  │  P_fratura: $(round(P_fratura, digits=1)) MPa  │  P_reativação: $(round(P_reativacao, digits=1)) MPa  │")
println("└─────────────────────────────────────────────────────────────────────────────┘")

println("\n┌─────────────────────────────────────────────────────────────────────────────┐")
println("│                    📏 ANÁLISE POR DISTÂNCIA INJETOR-FALHA                    │")
println("├──────────┬──────────┬──────────┬──────────┬──────────────────────────────────┤")
println("│ Distância│   ALTO   │  MÉDIO   │  BAIXO   │  Conclusão                       │")
println("├──────────┼──────────┼──────────┼──────────┼──────────────────────────────────┤")

for dist in sort(unique(df.dist_falha_m))
    subset = filter(row -> row.dist_falha_m == dist, df)
    a = count(==("ALTO"), subset.risco)
    m = count(==("MEDIO"), subset.risco)
    b = count(==("BAIXO"), subset.risco)
    
    if b == nrow(subset)
        conclusao = "✅ 100% SEGURO"
    elseif a > 0
        conclusao = "🚨 RISCO ALTO"
    else
        conclusao = "⚠️  RISCO MODERADO"
    end
    
    println("│  $(rpad(string(Int(dist))*"m", 8))│    $(rpad(string(a), 6))│    $(rpad(string(m), 6))│    $(rpad(string(b), 6))│  $(rpad(conclusao, 32))│")
end
println("└──────────┴──────────┴──────────┴──────────┴──────────────────────────────────┘")

println("\n┌─────────────────────────────────────────────────────────────────────────────┐")
println("│                         🔥 MAPA DE CALOR DO RISCO                            │")
println("├─────────────────────────────────────────────────────────────────────────────┤")
println("│                     Permeabilidade da Falha (mD)                             │")
println("│  Distância    10        100       1000      10000                           │")
println("│  ────────────────────────────────────────────────                            │")

for dist in sort(unique(df.dist_falha_m))
    print("│  $(rpad(string(Int(dist))*"m", 10))")
    for k_f in sort(unique(df.k_falha_mD))
        subset = filter(row -> row.dist_falha_m == dist && row.k_falha_mD == k_f, df)
        if nrow(subset) > 0
            if any(==("ALTO"), subset.risco)
                print("🔴 ")
            elseif any(==("MEDIO"), subset.risco)
                print("🟡 ")
            else
                print("🟢 ")
            end
        end
    end
    println("  │")
end
println("└─────────────────────────────────────────────────────────────────────────────┘")

println("\n📊 DISTRIBUIÇÃO DE RISCO POR DISTÂNCIA\n")

for dist in sort(unique(df.dist_falha_m))
    subset = filter(row -> row.dist_falha_m == dist, df)
    n_b = count(==("BAIXO"), subset.risco)
    n_m = count(==("MEDIO"), subset.risco)
    total = nrow(subset)
    
    pct_b = floor(Int, n_b/total*40)
    pct_m = floor(Int, n_m/total*40)
    
    println("  Distância: $(Int(dist))m")
    println("  🟢 BAIXO: " * "█"^max(0, pct_b) * " $n_b/$total")
    println("  🟡 MÉDIO: " * "█"^max(0, pct_m) * " $n_m/$total\n")
end

println("📈 EFEITO DO FATOR DE PRESSÃO\n")
for fp in sort(unique(df.fator_pressao))
    subset = filter(row -> row.fator_pressao == fp, df)
    n_b = count(==("BAIXO"), subset.risco)
    n_m = count(==("MEDIO"), subset.risco)
    P_media = mean(subset.P_falha_MPa)
    # Correção: usar round(Int, fp*100) para evitar InexactError
    println("  fp = $fp ($(round(Int, fp*100))% Pi):  P_média = $(round(P_media, digits=2)) MPa  |  BAIXO: $n_b  MÉDIO: $n_m")
end

println("\n📋 TABELA DE DECISÃO PARA OPERAÇÃO CCS\n")
println("┌─────────────────────────────────────────────────────────────────────┐")
println("│  CONDIÇÃO                    │  AÇÃO                                 │")
println("├─────────────────────────────────────────────────────────────────────┤")
println("│  Distância > 1700m           │  ✅ OPERAR (risco baixo)              │")
println("│  Distância 1100-1700m        │  ⚠️  OPERAR COM MONITORAMENTO         │")
println("│  Distância < 1100m           │  🔍 AVALIAR (risco médio)             │")
println("│  P_injeção > 1.20×Pi         │  ⚠️  REDUZIR PRESSÃO                   │")
println("│  S_falha > 0.5               │  🚨 INTERROMPER INJEÇÃO               │")
println("│  S_FSASP > 0                 │  🚨 VAZAMENTO - EVACUAR              │")
println("└─────────────────────────────────────────────────────────────────────┘")

println("\n🗺️  MODELO CONCEITUAL — CORTE VERTICAL\n")
nx, nz = 30, 20
println("┌" * "─"^nx * "┐")
for k in nz:-1:1
    print("│")
    for i in 1:nx
        if k <= 12
            if i == 15; print("║")
            elseif i == 3; print("↓")
            elseif i == 28; print("↑")
            else; print("·"); end
        elseif k <= 15
            if i == 15; print("║")
            else; print("▓"); end
        else
            print("░")
        end
    end
    println("│")
end
println("└" * "─"^nx * "┘")
println(" · Arenito  ▓ Selante  ░ FSASP  ║ Falha  ↓ Injetor  ↑ Produtor")

println("\n" * "═"^80)
println("  💡 RECOMENDAÇÕES PARA O CAMPO PILAR")
println("═"^80)
println("""
  1. Distância mínima segura: 1700m da falha
  2. Pressão de injeção: até 1.20×Pi é segura
  3. Selante efetivo em TODOS os cenários
  4. NENHUM vazamento para o FSASP
  5. Sem reativação (Mohr-Coulomb) com ΔP < 8.66 MPa
""")
println("═"^80)
println("  ✅ Dashboard concluído!")
println("═"^80)