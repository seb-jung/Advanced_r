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

## Paketdateien und Funktionen

Eine Datei je Teilaufgabe des Top-Down-Entwurfs. **Fett** = exportiert (mit Hilfeseite), alle anderen Funktionen sind intern.

| Datei (`R/`)         | Funktionen                                                                                                                                                       | Zweck                                                                                     |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| `poisselect.R`       | **`poisselect()`**                                                                                                                                               | Hauptfunktion, ruft die Teilschritte der Reihe nach auf                                   |
| `check_input.R`      | `check_arguments()`, `check_formula()`, `check_variables_available()`, `check_start_values()`, `check_model_data()` mit `check_selection_indicator()`, `check_outcome_counts()`, `report_ignored_outcomes()`, `check_design_matrix()`, `check_sample_size()`, `check_exclusion_restriction()` | Input-Checks: erst die Argumente, dann die abgeleiteten Daten                              |
| `model_data.R`       | `build_model_data()`, `build_equation_data()`, `as_selection_indicator()`                                                                                        | Formeln zu `y`, `s`, Designmatrizen, `terms`, Faktorstufen                                |
| `gauss_hermite.R`    | `build_gauss_hermite()`                                                                                                                                          | Knoten und Gewichte (Golub-Welsch)                                                        |
| `loglik.R`           | `compute_loglik()`, `compute_node_pieces()`, `compute_log_terms()`, `add_log_weight()`, `log_sum_exp_rows()`, `compute_loglik_gradient()`, `compute_negative_loglik()`, `compute_negative_gradient()` | Log-Likelihood auf der Log-Skala, analytischer Gradient, Zielfunktion für `optim()` |
| `parameters.R`       | `split_parameters()`, `pack_parameters()`, `build_parameter_names()`                                                                                             | Reparametrisierung `log(sigma)`, `atanh(rho)`                                             |
| `start_values.R`     | `compute_start_values()`, `fit_probit_start()`, `fit_glm_start()`, `estimate_start_sigma()`                                                                      | Startwerte aus Probit- und Poisson-GLM, Momentenschätzer für `sigma`                      |
| `optimise.R`         | `maximise_loglik()`, `warn_about_fit()`                                                                                                                          | `optim()` (BFGS), Warnungen bei Nichtkonvergenz und Randlösungen                          |
| `standard_errors.R`  | `compute_standard_errors()`, `invert_information()`, `extract_standard_errors()`                                                                                 | Hesse-Matrix, Inverse, Delta-Methode                                                      |
| `new_poisselect.R`   | `new_poisselect()`                                                                                                                                               | Konstruktor des S3-Objekts                                                                |
| `methods_print.R`    | **`print.poisselect()`**, `print_coefficient_block()`, `describe_convergence()`                                                                                  | Kurzausgabe                                                                               |
| `methods_summary.R`  | **`summary.poisselect()`**, **`print.summary.poisselect()`**, `build_coefficient_table()`, `print_coefficient_table()`, `describe_selection_bias()`               | Koeffiziententabellen mit z- und p-Werten, `rho`, Log-Likelihood, AIC                     |
| `methods_plot.R`     | **`plot.poisselect()`**, `plot_count_distribution()`, `compute_count_distribution()`, `plot_rho_profile()`, `compute_rho_profile()`                              | Plot 1: beobachtete vs. modell-implizierte Counts; Plot 2: Log-Likelihood entlang `rho`   |
| `methods_predict.R`  | **`predict.poisselect()`**, `build_prediction_matrix()`                                                                                                          | `link`, `response`, `pselect`, auch für `newdata`                                         |
| `methods_extract.R`  | **`coef.poisselect()`**, **`vcov.poisselect()`**, **`logLik.poisselect()`**                                                                                      | Extraktoren; `AIC()`/`BIC()` laufen über `logLik()`                                       |
| `simulate.R`         | **`simulate_poisselect()`**, `check_simulate_arguments()`                                                                                                        | Daten aus dem Modell simulieren                                                           |
| `data.R`             | Datensatz **`sim_selection`**                                                                                                                                    | Beispieldaten (n = 800) mit bekannten Parametern                                          |
| `poisselect-package.R` | keine                                                                                                                                                          | Paketdokumentation, `@importFrom`                                                         |

Tests liegen in `tests/testthat/` (6 Dateien, 340 Assertions).

## Wichtige Hinweise

- **Schnittstelle:** Die Funktion `poisselect()` folgt dem Standard von `lm()` oder `glm()`.
- **NA-Werte:** `NA`-Werte in `y` sind zulässig, sofern `s = 0`.
- **Exclusion Restriction:** Es sollte mindestens eine Variable in der Selektionsgleichung geben, die nicht in der Outcome-Gleichung vorkommt.
- **Genauigkeit:** Standardmäßig werden $K=20$ Gauß-Hermite-Knoten verwendet. Bei sehr hohen Zählwerten kann die Anzahl $K$ erhöht werden.
- **$\sigma$ und $\rho$:** Die Schätzung erfolgt über $\log(\sigma)$ und $\text{atanh}(\rho)$, um die Parameterbereiche korrekt zu berücksichtigen.

---

_Hinweis zur KI-Nutzung:_ Details zur Verwendung von KI-Werkzeugen finden Sie in der Datei `KI-Erklaerung.md`.

_Hinweise zur Dokumentation:_ Die exportierten Funktionen sind über die Hilfeseiten erreichbar (`?poisselect`, `?predict.poisselect` usw.). Die internen Hilfsfunktionen sind ohne eigene Hilfeseite direkt in den R-Dateien mit roxygen-Kommentaren dokumentiert; der Programmentwurf und die Begründung der numerischen Verfahren stehen in `Pseudocode_poisselect.pdf`.
