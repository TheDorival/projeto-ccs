using CSV, DataFrames, GLMakie

# ═══════════════════════════════════════════════
# CARREGAR DADOS
# ═══════════════════════════════════════════════
df = CSV.read("sensibilidade_campo_pilar.csv", DataFrame)
dists = sort(unique(df.dist_falha_m))
ks    = sort(unique(df.k_falha_mD))

# ═══════════════════════════════════════════════
# 1. MAPA DE CALOR (Saturação na Falha)
# ═══════════════════════════════════════════════
fig1 = Figure(size=(800,600))
ax1 = Axis(fig1[1,1],
    xlabel="Distância do injetor à falha (m)",
    ylabel="Permeabilidade da falha (mD)",
    title="Saturação de CO₂ na Falha")

M = zeros(length(ks), length(dists))
for (i,k) in enumerate(ks), (j,d) in enumerate(dists)
    sub = df[(df.k_falha_mD .== k) .& (df.dist_falha_m .== d), :]
    if nrow(sub) > 0
        M[i,j] = sub.S_falha[1]
    end
end

hm = heatmap!(ax1, dists, ks, M, colormap=:plasma)
Colorbar(fig1[1,2], hm, label="Saturação")
save("heatmap_sfalha.png", fig1)
display(fig1)

# ═══════════════════════════════════════════════
# 2. GRÁFICO DE BARRAS AGRUPADAS (substitui empilhado)
# ═══════════════════════════════════════════════
fig2 = Figure(size=(800,500))
ax2 = Axis(fig2[1,1],
    xlabel="Distância (m)", ylabel="Número de cenários",
    title="Classificação de Risco por Distância")

baixo = [count(==("BAIXO"), df[df.dist_falha_m .== d, :risco]) for d in dists]
medio = [count(==("MEDIO"), df[df.dist_falha_m .== d, :risco]) for d in dists]
alto  = [count(==("ALTO"),  df[df.dist_falha_m .== d, :risco]) for d in dists]

# Barras agrupadas (dodge)
barplot!(ax2, Float32.(dists) .- 0.15, Float32.(baixo), width=0.12, color=:green, label="BAIXO")
barplot!(ax2, Float32.(dists),         Float32.(medio), width=0.12, color=:orange, label="MÉDIO")
barplot!(ax2, Float32.(dists) .+ 0.15, Float32.(alto),  width=0.12, color=:red, label="ALTO")
axislegend(ax2, position=:rt)
save("barras_risco.png", fig2)
display(fig2)

# ═══════════════════════════════════════════════
# 3. DISPERSÃO: Pressão × Saturação
# ═══════════════════════════════════════════════
fig3 = Figure(size=(800,600))
ax3 = Axis(fig3[1,1],
    xlabel="Pressão na falha (MPa)",
    ylabel="Saturação de CO₂ na falha",
    title="Pressão vs Saturação - Todos os Cenários")

cores = [:red, :orange, :green]
for (i,d) in enumerate(dists)
    sub = df[df.dist_falha_m .== d, :]
    scatter!(ax3, sub.P_falha_MPa, sub.S_falha,
        color=cores[i], label="$(Int(d)) m", markersize=12)
end
vlines!(ax3, [30.28], color=:black, linestyle=:dash, label="P_reat")
axislegend(ax3)
save("dispersao_pressao_sat.png", fig3)
display(fig3)

# ═══════════════════════════════════════════════
# 4. EXPLORADOR INTERATIVO COM SLIDERS
# ═══════════════════════════════════════════════
fig4 = Figure(size=(1000,600))

sg = SliderGrid(fig4[1:2,1],
    (label="Distância (m)", range=500:100:1700, startvalue=1100),
    (label="k_falha (mD)", range=[10,100,1000,10000], startvalue=100),
    (label="Fator de Pressão", range=1.05:0.01:1.20, startvalue=1.10),
    tellwidth=false)

ax_press = Axis(fig4[1,2], xlabel="Cenário", ylabel="Pressão na falha (MPa)")
ax_sat   = Axis(fig4[2,2], xlabel="Cenário", ylabel="Saturação na falha")

on(sg.sliders[1].value) do d
    on(sg.sliders[2].value) do k
        on(sg.sliders[3].value) do fp
            sub = df[(df.dist_falha_m .== d) .& (df.k_falha_mD .== k) .& (df.fator_pressao .== fp), :]
            if nrow(sub) > 0
                empty!(ax_press)
                empty!(ax_sat)
                scatter!(ax_press, [1], [sub.P_falha_MPa[1]], markersize=20, color=:blue)
                scatter!(ax_sat,   [1], [sub.S_falha[1]],   markersize=20, color=:red)
                ylims!(ax_press, 21, 35)
                ylims!(ax_sat,   0,  1)
                xlims!(ax_press, 0, 2)
                xlims!(ax_sat,   0, 2)
            end
        end
    end
end

save("explorador_cenarios.png", fig4)
display(fig4)

println("✅ Dashboard interativo gerado com sucesso!")