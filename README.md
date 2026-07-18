# Projeto CCS — Campo Pilar (Sub-bacia de Alagoas)

Modelagem numérica de armazenamento geológico de CO₂ (CCS) e gás natural (AGSS) no Campo Pilar, Sub-bacia de Alagoas (Bacia SEAL), usando dados reais do poço **1-PIR-0001-AL** e simulação de fluxo com **Julia + [Jutul/JutulDarcy](https://github.com/sintefmath/JutulDarcy.jl)**.

**Pergunta central:** a injeção de CO₂ no intervalo selecionado eleva a pressão na falha próxima a ponto de reativá-la ou comprometer a contenção pelo capeador?

**Resposta (resumo):** viável com pressão de injeção ≤ 1,10 × P_ini. Com fator 1,20 a falha reativa no cenário geomecânico pessimista. A contenção do capeador intacto depende criticamente da pressão de entrada capilar do selo (≥ 5 MPa → invasão nula). Detalhes em [`RELATORIO_INTERVALO_INJECAO.md`](RELATORIO_INTERVALO_INJECAO.md) e no relatório completo [`Relatorio_Campo_Pilar_CCS.docx`](Relatorio_Campo_Pilar_CCS.docx).

## Intervalo de injeção

A busca automática (janelas de 120 m, faixa 800–4000 m) seleciona **3330–3440 m na Fm. Penedo** (injeção a 3390 m): KX médio 79,7 mD, PHIE 16,8 %, 11/12 arenito, 87 °C, P_ini hidrostática 33,2 MPa — CO₂ supercrítico com densidade ~650 kg/m³. O selo regional é a Fm. Coqueiro Seco (315–2240 m, pelítica).

## Estrutura do código

| Arquivo | Função |
|---|---|
| `campo_pilar_common.jl` | Módulo comum (incluído pelos demais): seleção da janela ótima, carga dos dados do poço (KX, **KZ**, temperatura, geomecânica por cenário), domínio com tensor de permeabilidade `[kx, kx, kz]`, pressão capilar de Brooks-Corey **por região** (reservatório × capeador), viscosidades, análise de *slip tendency* e o cenário CCS base (`rodar_cenario`) |
| `simulacao_teste1.jl` | **Sensibilidade CCS**: 135 cenários — 108 de falha/pressão/geomecânica (3 distâncias × 4 k_falha × 3 fatores de pressão × 3 cenários geomec.) + 27 de selante/FSASP (k_selante × Pc_selo × k_FSASP). Gera `sensibilidade_campo_pilar.csv` |
| `simulacao_agss.jl` | **AGSS**: ciclos sazonais de gás natural (6 meses injeção / 6 retirada, 10 anos, passos mensais) com **histerese de Killough** e aquífero lateral de pressão constante. Gera `agss_campo_pilar.csv` |
| `simulacao_co2brine.jl` | **CCS com PVT real**: sistema `:co2brine` (compositional k-value, correlações Salo 2024), dissolução mútua CO₂-salmoura, temperatura real do poço, 20 anos de injeção + 30 de pós-injeção. Gera `co2brine_campo_pilar.csv` |
| `validacao_convergencia.jl` | **Validação**: convergência de malha (30×20 → 60×40 → 90×60) com geometria física preservada. Gera `validacao_convergencia.csv` |
| `dashboard_results.jl` | Dashboard textual do CSV de sensibilidade, com recomendações derivadas dos dados |
| `graficos_makie.jl` | Gráficos (GLMakie): heatmap de saturação na falha, risco por distância, pressão × saturação, contenção do capeador e explorador interativo |
| `malha_1PIR_enriquecida.csv` | Dados do poço 1PIR: 253 blocos de 10 m (315–3980 m) com VSH, PHIE, KX, KZ, fácies, geomecânica (pessimista/central/otimista), Pc, Swirr, Sgr, temperatura, formação |
| `CODIGO JUTULDARCY V1.txt` | Versão histórica original do script (referência) |

## Modelo

Seção 2D vertical 3000 × 200 m (padrão 30 × 1 × 20 células): reservatório de 120 m com propriedades reais bloco a bloco, selante sintético de 30 m (o selo local não é observável — lacuna de dados 2480–3240 m no poço), FSASP acima, e falha vertical conduíte parametrizada. As unidades são definidas por **profundidade física**, permitindo refinar a malha sem alterar a geometria.

Geomecânica com dois critérios em paralelo: Mohr-Coulomb simplificado (P_reativ = P_ini + c/tanφ; pessimista 38,6 / central 51,2 MPa) e *slip tendency* no plano da falha (σv integrado da densidade bulk do poço, σh = K0·σv′, mergulho 60° → ST = 0,17; P_crit = 80 MPa). Risco: **ALTO** (CO₂ acima do selo ou capeador intacto invadido), **MÉDIO** (reativação, CO₂ na falha ou margem de fratura < 5 MPa), **BAIXO**.

## Como rodar

```bash
# 1ª vez: instalar todas as dependências (Project.toml/Manifest.toml)
julia --project=. -e "using Pkg; Pkg.instantiate()"

# GLMakie (fora do Project.toml; necessário apenas para graficos_makie.jl)
# No PowerShell use: julia --project=. -e 'using Pkg; Pkg.add(\"GLMakie\")'
# ou no REPL: julia --project=.  →  tecla ]  →  add GLMakie
julia --project=. -e 'using Pkg; Pkg.add("GLMakie")'

# Sensibilidade CCS (valida a base; ~15 min na 1ª execução por precompilação)
julia --project=. simulacao_teste1.jl

# Demais cenários (independentes entre si)
julia --project=. simulacao_agss.jl
julia --project=. simulacao_co2brine.jl
julia --project=. validacao_convergencia.jl

# Pós-processamento
julia --project=. dashboard_results.jl
julia --project=. graficos_makie.jl   # requer GLMakie
```

Requisitos: Julia ≥ 1.10 (testado em 1.12), pacotes do `Project.toml`/`Manifest.toml` + CSV e DataFrames. No PowerShell, prefira o REPL (`julia --project=.` e tecla `]`) para comandos Pkg com strings.

## Resultados principais (rodada 18/07/2026)

- **Reativação**: 4/135 cenários, todos PESSIMISTA com fator de pressão 1,20 (anos 3–8). Nenhuma com fator ≤ 1,10; nenhuma por slip tendency.
- **Contenção**: Pc de entrada do capeador é o parâmetro decisivo (5 MPa → invasão zero; 1–2 MPa → S até 0,015). k do FSASP irrelevante. Prioridade experimental: ensaio MICP no selo.
- **Capacidade/injetividade** (2D, por metro de largura): 6–95 kt em 20 anos conforme o fator; II ≈ 390–1950 kg/dia/MPa.
- **AGSS**: eficiência de working gas 44–77 % (máx. com BHP 1,10/0,85 × P_ini); sem reativação (P_falha ≤ 35,3 MPa).
- **PVT real**: capeador intacto contém 100 % do CO₂; migração só pela falha-conduíte; 9,8 % dissolvido aos 50 anos (crescente).
- **Convergência**: P_falha converge (1,1 %); malha grossa superestima a invasão do capeador em ~2 ordens de grandeza — usar 60×40 para resultados finais.

## Limitações

Modelo 2D (massas por metro de largura); selo local não observado (lacuna 2480–3240 m); geomecânica desacoplada (acoplamento previsto no cronograma mar–jun/2026); gradiente de fratura e K0 assumidos; sem trapeamento mineral.
