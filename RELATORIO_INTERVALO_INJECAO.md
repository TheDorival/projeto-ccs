# Viabilidade de estocagem de CO₂ — Campo Pilar (Sub-bacia de Alagoas)
## Análise do melhor intervalo de injeção — poço 1PIR-0001-AL

Data: 2026-07-18

## 1. Melhor intervalo de profundidade: 3330–3440 m (Fm. Penedo)

A varredura de todas as janelas de 12 blocos × 10 m sobre a malha enriquecida do 1PIR (315–3980 m) identifica a zona 3200–3450 m da Formação Penedo como a de melhor qualidade de reservatório. A janela ótima é **3330,0–3440,1 m**, com profundidade de injeção de **3390,05 m**:

| Critério | Janela ótima (3330–3440 m) | Melhor janela na faixa antiga 800–2500 m |
|---|---|---|
| KX médio | **79,7 mD** | 33,4 mD |
| PHIE médio | **16,8 %** | 9,9 % |
| Arenito (FACIES=1) | 11/12 | 10/12 |
| VSH médio | 0,17 | 0,27 |
| Formação | Penedo | Penedo (topo) |
| Temperatura | ~87 °C | ~64 °C |

Justificativas adicionais:

- **CO₂ supercrítico garantido**: a 3390 m, P ≈ 33,2 MPa e T ≈ 87 °C, muito acima do ponto crítico (7,38 MPa / 31 °C). Densidade real do CO₂ nessas condições ≈ 650–660 kg/m³, consistente com o valor usado no código (650).
- **Selo regional**: a Fm. Coqueiro Seco (315–2240 m) é predominantemente pelítica (VSH médio ~0,7, KX ~10⁻⁴ mD nos folhelhos) e atua como selo regional acima da Penedo.
- **Margem à fratura ampla**: P_fratura (0,0181 MPa/m) ≈ 61,4 MPa a 3390 m, ~28 MPa acima da pressão máxima de injeção simulada.
- A restrição anterior do código (800–2500 m) selecionaria uma janela com menos da metade da permeabilidade e da porosidade — o CSV de sensibilidade existente (prof. 3390,05 m) já havia sido gerado com o intervalo correto.

**Mudança aplicada ao código** (`simulacao_teste1.jl`): `profundidade_max` da busca ampliada de 2500 → 4000 m. Validado por réplica da lógica em Python: o seletor passa a escolher exatamente 3330–3440,1 m → prof. 3390,05 m, idêntico ao usado em `sensibilidade_campo_pilar.csv`.

## 2. Viabilidade a 3390 m — resultados da simulação corrigida (108 cenários)

Rodada de 18/07/2026 com o código corrigido (P_ini = 33,20 MPa hidrostático, Pc ativa, curvas Brooks-Corey com Swirr/Sgr efetivamente aplicadas, API JutulDarcy atual). 20 anos, 3 distâncias de falha × 4 permeabilidades de falha × 3 fatores de pressão × 3 cenários geomecânicos:

- **Reativação da falha: 7 de 108 cenários — todos PESSIMISTA com fator 1,20** (pico no ano 3–5), confirmando a previsão analítica (margem ≈ c/tanφ − (fp−1)·P_ini ≈ −2 MPa nesse caso). Nenhuma reativação com fator ≤ 1,10, em nenhum cenário geomecânico.
- **Risco**: 0 BAIXO / 86 MÉDIO / 22 ALTO. O CO₂ alcança a falha em todos os cenários (S_falha > 0,01 → MÉDIO). Os casos ALTO (migração pela falha até acima do selo, S_FSASP > 10⁻⁴) concentram-se em falha pouco permeável (10 mD) com fator ≥ 1,10 e em fator 1,20 com k_falha ≤ 100 mD.
- **Massa injetada em 20 anos**: ~7–32 kt (fator 1,05), ~15–66 kt (1,10), ~40–150 kt (1,20). Margem à fratura sempre > 20 MPa (P_fratura ≈ 61,4 MPa).

**Conclusão**: a estocagem a 3390 m é viável com **fator de pressão ≤ 1,10** (nenhuma reativação); a 1,05 o risco de migração pela falha também é menor. Fator 1,20 é inaceitável: reativa a falha no cenário geomecânico pessimista e promove migração de CO₂ acima do selo. A distância injetor–falha e a permeabilidade da falha modulam a massa armazenável, mas o fator de pressão é o controle dominante do risco.

## 4. Modelo completo (rodada de 18/07/2026, versão com KZ, Pc por região, AGSS, PVT real e validação)

O modelo foi estendido para cobrir todas as atividades do cronograma: tensor de permeabilidade com KZ real do poço, pressão capilar por região (capeador com entrada própria), injetividade, slip tendency, ciclos AGSS com histerese Killough, sistema CO₂-brine com PVT real e validação por convergência de malha. Principais resultados:

**CCS (135 cenários):** o controle dominante do risco segue sendo o fator de pressão; reativação só no PESSIMISTA a 1,20. A varredura do capeador mostrou que a **pressão de entrada capilar do selo é o parâmetro decisivo de contenção**: com Pc_selo = 5 MPa não há invasão do capeador intacto em nenhum caso; com 1–2 MPa há invasão incipiente (S até 0,015) que cresce com k_selo. O k do FSASP é irrelevante para a contenção (S_FSASP = 0 em toda a varredura do selo). Recomenda-se medir a Pc de entrada do capeador (ensaio MICP) — é o dado de laboratório de maior valor para reduzir a incerteza.

**AGSS (gás natural, ciclos 6+6 meses, 10 anos, histerese Killough):** operação viável sem reativação (P_falha ≤ 35,3 MPa ≪ P_reativ 51 MPa). Eficiência de working gas de 44–77% por ciclo, maximizada com maior amplitude de pressão (1,10/0,85 → 77%); o primeiro ciclo forma o cushion gas (gás aprisionado ≈ SGR). Valores de massa são por metro de largura do modelo 2D.

**CCS com PVT real (:co2brine, 20 anos + 30 pós-injeção):** capeador intacto contém totalmente o CO₂ (S_selante = 0); a migração acima do selo ocorre exclusivamente pela falha-conduíte (S_FSASP ≈ 0,18, decrescente no pós-injeção). Trapeamento por dissolução: 9,8% da massa aos 50 anos, crescente — mecanismo de segurança de longo prazo que o modelo imiscível não capturava.

**Validação (malhas 30×20 → 60×40 → 90×60):** P_falha converge (variação 1,1%); a massa decresce monotonicamente com refinamento (27,5 → 24,7 → 23,6 kt, diferenças decrescentes) e a invasão do capeador é superestimada na malha grossa (S_selante 10⁻³ → 10⁻⁵). Para resultados finais de publicação, usar a malha 60×40 (custo baixo após compilação).

## 5. Limitações e próximos passos

1. **Lacunas de dados no poço**: 970–1380 m e 2480–3240 m sem registros. A segunda lacuna fica logo acima da janela de injeção; o selo intra-Penedo local não é observável nos dados — o modelo usa um selo sintético de 30 m (10⁻⁴ mD). Recomenda-se caracterizar o intervalo 2480–3240 m (sísmica/poços vizinhos) para confirmar o selo local.
2. O critério de reativação é simplificado (não resolve tensor de tensões nem orientação da falha); a análise geomecânica acoplada está prevista no cronograma (mar–jun/2026).
3. Modelo 2D (ny=1), imiscível, sem dissolução/trapeamento mineral — conservador quanto a capacidade, otimista quanto a migração lateral.
