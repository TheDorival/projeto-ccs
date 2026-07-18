# ═══════════════════════════════════════════════════════════════════
# GRÁFICOS — CAMPO PILAR CCS (GLMakie)
# Compatível com o CSV novo (varreduras de falha/pressão e selo/FSASP)
# ═══════════════════════════════════════════════════════════════════

using CSV, DataFrames, Statistics, GLMakie

df = CSV.read(joinpath(@__DIR__, "sensibilidade_campo_pilar.csv"), DataFrame)

eh_default_selo(r) = r.k_selante_mD == 0.0001 && r.pc_selante_MPa == 2.0 && r.k_FSASP_mD == 100.0
df_falha = filter(eh_default_selo, df)
df_selo  = filter(r -> !eh_default_selo(r), df)

dists = sort(unique(df_falha.dist_falha_m))
ks    = sort(unique(df_falha.k_falha_mD))

# ═══════════════════════════════════════════════
# 1. MAPA DE CALOR — saturação máxima na falha (pior caso da célula)
# ═══════════════════════════════════════════════
fig1 = Figure(size=(800,600))
ax1 = Axis(fig1[1,1],
    xlabel="Distância do injetor à falha (m)",
    ylabel="Permeabilidade da falha (mD)",
    yticks=(1:length(ks), string.(Int.(ks))),
    xticks=(1:length(dists), string.(Int.(dists))),
    title="Saturação máxima de CO₂ na falha (pior caso)")

M = zeros(length(dists), length(ks))
for (j,k) in enumerate(ks), (i,d) in enumerate(dists)
    sub = df_falha[(df_falha.k_falha_mD .== k) .& (df_falha.dist_falha_m .== d), :]
    M[i,j] = nrow(sub) > 0 ? maximum(sub.S_falha) : NaN
end
hm = heatmap!(ax1, 1:length(dists), 1:length(ks), M, colormap=:plasma)
Colorbar(fig1[1,2], hm, label="Saturação")
save(joinpath(@__DIR__, "heatmap_sfalha.png"), fig1)
display(fig1)

# ═══════════════════════════════════════════════
# 2. BARRAS AGRUPADAS — risco por distância (dodge por índice)
# ═══════════════════════════════════════════════
fig2 = Figure(size=(800,500))
ax2 = Axis(fig2[1,1],
    xlabel="Distância injetor-falha (m)", ylabel="Número de cenários",
    xticks=(1:length(dists), string.(Int.(dists))),
    title="Classificação de risco por distância (varredura falha/pressão)")

baixo = [count(==("BAIXO"), df_falha[df_falha.dist_falha_m .== d, :risco]) for d in dists]
medio = [count(==("MEDIO"), df_falha[df_falha.dist_falha_m .== d, :risco]) for d in dists]
alto  = [count(==("ALTO"),  df_falha[df_falha.dist_falha_m .== d, :risco]) for d in dists]

x = 1:length(dists)
barplot!(ax2, x .- 0.25, Float64.(baixo), width=0.22, color=:green,  label="BAIXO")
barplot!(ax2, x,         Float64.(medio), width=0.22, color=:orange, label="MÉDIO")
barplot!(ax2, x .+ 0.25, Float64.(alto),  width=0.22, color=:red,    label="ALTO")
axislegend(ax2, position=:rt)
save(joinpath(@__DIR__, "barras_risco.png"), fig2)
display(fig2)

# ═══════════════════════════════════════════════
# 3. DISPERSÃO — pressão × saturação na falha, com limiares reais
# ═══════════════════════════════════════════════
fig3 = Figure(size=(800,600))
ax3 = Axis(fig3[1,1],
    xlabel="Pressão na falha (MPa)",
    ylabel="Saturação de CO₂ na falha",
    title="Pressão × saturação — todos os cenários da varredura 1")

cores = [:red, :orange, :green]
for (i,d) in enumerate(dists)
    sub = df_falha[df_falha.dist_falha_m .== d, :]
    scatter!(ax3, sub.P_falha_MPa, sub.S_falha,
        color=cores[i], label="$(Int(d)) m", markersize=10)
end
P_reativ_pess = minimum(df_falha.P_reativ_MPa)
vlines!(ax3, [P_reativ_pess], color=:black, linestyle=:dash,
    label="P_reativ pessimista ($(round(P_reativ_pess, digits=1)) MPa)")
axislegend(ax3, position=:lt)
save(joinpath(@__DIR__, "dispersao_pressao_sat.png"), fig3)
display(fig3)

# ═══════════════════════════════════════════════
# 4. CONTENÇÃO DO CAPEADOR — S_selante × (k_selante, Pc_selo)
# ═══════════════════════════════════════════════
if nrow(df_selo) > 0
    ksel = sort(unique(df_selo.k_selante_mD))
    pcs  = sort(unique(df_selo.pc_selante_MPa))
    fig4 = Figure(size=(800,500))
    ax4 = Axis(fig4[1,1],
        xlabel="Pc de entrada do capeador (MPa)",
        ylabel="S máx. no capeador intacto (log)",
        yscale=log10,
        title="Contenção do capead