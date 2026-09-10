# poisselect: Poisson-Selektionsmodell für Zähldaten

Das R-Paket `poisselect` schätzt ein Heckman-artiges Selektionsmodell für Zähldaten mittels Maximum Likelihood. Es ist besonders geeignet, wenn die Wahrscheinlichkeit, Daten zu beobachten (Selektion), mit den tatsächlichen Zählwerten korreliert.

## Inhalt der Abgabe

| Datei                                       | Inhalt                                                                                                                        |
| ------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| `poisselect_0.1.0.tar.gz`                   | das gebaute Paket, so installieren wie unten gezeigt                                                                          |
| `poisselect/`                               | der Quellbaum des Pakets (identisch mit dem Inhalt des `tar.gz`, plus `data-raw/`)                                            |
| `README.md`                                 | diese Datei                                                                                                                   |
| `Pseudocode_poisselect.pdf` (Quelle `.qmd`) | Programmentwurf nach dem Top-Down-Workflow der Vorlesung                                                                      |
| `KI-Erklaerung.pdf` (Quelle `.md`)          | Kennzeichnung der KI-Nutzung gemäß Abschnitt 4 der KI-Richtlinie                                                              |
| `beispiel_end_to_end.R`                     | Skript, das Schätzung, `print`, `summary`, beide Plots, alle `predict`-Typen sowie `coef`, `vcov`, `logLik`, `AIC` durchläuft |

## Installation

```r
remotes::install_local("poisselect_0.1.0.tar.gz",
                       dependencies = TRUE,
                       type = "source")
```

_Abhängigkeiten: `checkmate`, Basispakete: `stats, graphics, utils`, `testthat (nur für Tests notwendig)`._

## Schnellstart

```r
library(poisselect)

# Beispiel: Schätzung auf simulierten Daten
fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = sim_selection)

summary(fit)      # Koeffizienten, SE, rho, AIC
plot(fit)         # Diagnoseplots
predict(fit)      # Erwartungswert E[Y|x]
predict(fit, type="pselect") # Selektionswahrscheinlichkeit
```

## Modellbeschreibung

Das Modell besteht aus einer Outcome-Gleichung für die Zählvariable $y$ und einer Selektionsgleichung für den binären Indikator $s$:

- **Outcome:** $\ln \mu_i = x_i' \beta + \epsilon_i, \quad y_i \sim \text{Poisson}(\mu_i)$
- **Selektion:** $s_i^* = z_i' \gamma + u_i, \quad s_i = \mathbf{1}\{s_i^* > 0\}$
- **Fehlerterme:** $(\epsilon_i, u_i) \sim N_2(0, \begin{pmatrix} \sigma^2 & \rho \sigma \\ \rho \sigma & 1 \end{pmatrix})$

Der Parameter $\rho$ misst die Korrelation zwischen den Fehlertermen und damit die Stärke der Selektionsverzerrung.

## Wichtige Hinweise

- **Schnittstelle:** Die Funktion `poisselect()` folgt dem Standard von `lm()` oder `glm()`.
- **NA-Werte:** `NA`-Werte in `y` sind zulässig, sofern `s = 0`.
- **Exclusion Restriction:** Es sollte mindestens eine Variable in der Selektionsgleichung geben, die nicht in der Outcome-Gleichung vorkommt.
- **Genauigkeit:** Standardmäßig werden $K=20$ Gauß-Hermite-Knoten verwendet. Bei sehr hohen Zählwerten kann die Anzahl $K$ erhöht werden.
- **$\sigma$ und $\rho$:** Die Schätzung erfolgt über $\log(\sigma)$ und $\text{atanh}(\rho)$, um die Parameterbereiche korrekt zu berücksichtigen.

---

_Hinweis zur KI-Nutzung:_ Details zur Verwendung von KI-Werkzeugen finden Sie in der Datei `KI-Erklaerung.md`.

_Hinweise zur Dokumentation:_ Die exportierten Funktionen sind über die Hilfeseiten erreichbar (`?poisselect`, `?predict.poisselect` usw.). Die internen Hilfsfunktionen sind ohne eigene Hilfeseite direkt in den R-Dateien mit roxygen-Kommentaren dokumentiert; der Programmentwurf und die Begründung der numerischen Verfahren stehen in `Pseudocode_poisselect.pdf`.
