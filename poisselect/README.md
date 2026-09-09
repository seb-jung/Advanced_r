# poisselect: Poisson-Selektionsmodell für Zähldaten

Hausarbeit im Kurs _Advanced Data Analysis with R_.
Verfasser: Sebastian Jung, Finn Metzler.

Das Paket schätzt ein Heckman-artiges Selektionsmodell für Zähldaten per
Maximum Likelihood. Eine Poisson-Outcome-Gleichung mit lognormaler
unbeobachteter Heterogenität wird gemeinsam mit einer
Probit-Selektionsgleichung geschätzt; die Korrelation `rho` der beiden
Fehlerterme korrigiert die Verzerrung, der ein einfaches Poisson-GLM auf der
selbstselektierten Teilstichprobe unterliegt.

Die Schätzung ist vollständig selbst implementiert: Gauß-Hermite-Knoten
(Golub-Welsch), Log-Likelihood auf der Log-Skala, analytischer Gradient,
Maximierung mit `optim()`, Standardfehler aus der Hesse-Matrix. Es wird kein
fertiges Selektionsmodell-Paket verwendet.

## Inhalt der Abgabe

| Datei                                       | Inhalt                                                                                                                        |
| ------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| `poisselect_0.1.0.tar.gz`                   | das gebaute Paket, so installieren wie unten gezeigt                                                                          |
| `poisselect/`                               | der Quellbaum des Pakets (identisch mit dem Inhalt des `tar.gz`, plus `LICENSE.md`, `README.md`, `data-raw/`)                 |
| `README.md`                                 | diese Datei                                                                                                                   |
| `Pseudocode_poisselect.pdf` (Quelle `.qmd`) | Programmentwurf nach dem Top-Down-Workflow der Vorlesung                                                                      |
| `KI-Erklaerung.pdf` (Quelle `.md`)          | Kennzeichnung der KI-Nutzung gemäß Abschnitt 4 der KI-Richtlinie                                                              |
| `beispiel_end_to_end.R`                     | Skript, das Schätzung, `print`, `summary`, beide Plots, alle `predict`-Typen sowie `coef`, `vcov`, `logLik`, `AIC` durchläuft |
| `BEFUND.md`, `MERGE_ENTSCHEIDUNGEN.md`      | Prüfprotokoll und Begründung der Zusammenführung zweier Vorfassungen                                                          |

## Installation

```r
remotes::install_local("poisselect_0.1.0.tar.gz",
                       dependencies = TRUE,
                       type = "source")
```

Getestet mit R 4.5.1 (macOS, arm64). Laufzeit-Abhängigkeiten: `checkmate`
(>= 2.1.0) sowie die Basispakete `stats`, `graphics` und `utils`. Für die
Tests zusätzlich `testthat` (>= 3.0.0). `dependencies = TRUE` installiert
`checkmate` bei Bedarf von CRAN nach; dafür ist eine Internetverbindung nötig.

## Schnellstart

```r
library(poisselect)

# Mitgelieferter simulierter Datensatz mit bekannten Parametern
# (beta = c(0.5, 0.8, -0.4), gamma = c(0.3, 0.5, 0.7), sigma = 0.6, rho = 0.5)
fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = sim_selection)

fit                                   # Kurzübersicht
summary(fit)                          # Koeffiziententabellen, rho, LL, AIC
plot(fit)                             # beide Diagnoseplots

head(predict(fit))                    # E[Y | x] = exp(x'beta + sigma^2 / 2)
head(predict(fit, type = "link"))     # x'beta
head(predict(fit, type = "pselect"))  # Phi(z'gamma)

coef(fit); vcov(fit); logLik(fit); AIC(fit)

# Eigene Daten aus dem Modell ziehen
set.seed(42)
daten <- simulate_poisselect(n = 2000, rho = 0.6)
poisselect(y ~ x1 + x2, s ~ x1 + z1, data = daten, K = 40)
```

## Modell

```
Outcome:    ln(mu_i) = x_i'beta + eps_i,     y_i ~ Poisson(mu_i)
Selektion:  s_i*     = z_i'gamma + u_i,      s_i = 1{s_i* > 0}
(eps_i, u_i) ~ N2(0, [[sigma^2, rho sigma], [rho sigma, 1]])
```

`y_i` wird nur für `s_i = 1` beobachtet, `s_i` ist immer beobachtet. Das
Integral der Log-Likelihood wird per Gauß-Hermite-Quadratur mit `K` Knoten
approximiert (Default `K = 20`) und vollständig auf der Log-Skala mit dem
Log-Sum-Exp-Trick ausgewertet.

## Paketdateien

Der Code folgt dem Programmentwurf in `Pseudocode_poisselect.pdf`: eine
Hauptfunktion, die nur benannte Teilaufgaben aufruft, und je eine Datei pro
Teilaufgabe.

| Datei                                      | Inhalt                                                                                                  |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------------- |
| `R/poisselect.R`                           | Hauptfunktion `poisselect()`, ruft die Teilschritte in Reihenfolge auf                                  |
| `R/check_input.R`                          | Alle Input-Checks, zweistufig: erst die Argumente, dann die abgeleiteten Daten                          |
| `R/model_data.R`                           | Formeln zu `y`, `s`, Designmatrizen, `terms`, Faktorstufen; Aufteilung in selektiert / nicht selektiert |
| `R/gauss_hermite.R`                        | Knoten und Gewichte per Golub-Welsch-Eigenwertmethode                                                   |
| `R/loglik.R`                               | Log-Likelihood auf der Log-Skala, ihr analytischer Gradient, Zielfunktion für `optim()`                 |
| `R/parameters.R`                           | Reparametrisierung `log(sigma)`, `atanh(rho)`; Parameternamen                                           |
| `R/start_values.R`                         | Startwerte aus Probit- und Poisson-GLM, Momentenschätzer für `sigma`, Prüfung auf perfekte Trennung     |
| `R/optimise.R`                             | `optim()`-Aufruf (BFGS) und Warnungen bei Nichtkonvergenz oder Randlösungen                             |
| `R/standard_errors.R`                      | Hesse-Matrix per `optimHess()`, Inverse, Delta-Methode                                                  |
| `R/new_poisselect.R`                       | Konstruktor des S3-Objekts                                                                              |
| `R/methods_print.R`, `R/methods_summary.R` | `print()`, `summary()` und `print.summary()`                                                            |
| `R/methods_plot.R`                         | beide Plots und ihre Berechnungen                                                                       |
| `R/methods_predict.R`                      | `predict()` mit `link`, `response`, `pselect`                                                           |
| `R/methods_extract.R`                      | `coef()`, `vcov()`, `logLik()`                                                                          |
| `R/simulate.R`                             | `simulate_poisselect()`                                                                                 |
| `R/data.R`                                 | Dokumentation des Datensatzes `sim_selection`                                                           |
| `R/poisselect-package.R`                   | Paketdokumentation, `@importFrom`                                                                       |
| `data/sim_selection.rda`                   | Beispieldatensatz (n = 800), erzeugt von `data-raw/sim_selection.R`                                     |
| `tests/testthat/`                          | 6 Testdateien, 340 Assertions                                                                           |

Exportiert sind nur `poisselect()`, `simulate_poisselect()` und die acht
S3-Methoden; alle Hilfsfunktionen sind intern (`@noRd`).

## Was zu beachten ist

**Schnittstelle.** `poisselect()` nimmt zwei zweiseitige Formeln und einen
`data.frame`, wie `lm()`/`glm()`. Der Selektionsindikator ist die linke Seite
von `selection` (0/1 oder logisch), die Zählvariable die linke Seite von
`outcome`. Faktoren, Interaktionen und `predict(newdata = ...)` funktionieren
über `model.matrix()`. Das Argument für die Knotenzahl heißt `K` wie in der
Aufgabenstellung; das ist die einzige Abweichung vom `snake_case`-Styleguide
und an der betreffenden Zeile als Ausnahme markiert.

**NA-Muster.** `y` darf und soll `NA` sein, wo `s = 0`. Ist `y` dort trotzdem
beobachtet, wird es ignoriert und eine `message()` weist darauf hin. Die
Kovariablen der Selektionsgleichung müssen für alle Einheiten vollständig
sein, die der Outcome-Gleichung nur für die selektierten Einheiten.

**Exclusion Restriction.** Enthält die Selektionsgleichung keine Variable,
die nicht auch in der Outcome-Gleichung steht (der Intercept zählt nicht),
gibt es eine Warnung: das Modell ist dann nur über die Funktionalform
identifiziert. Plot 2 macht das sichtbar, eine flache Kurve bedeutet ein
schwach identifiziertes `rho`.

**Genauigkeit der Quadratur.** `K = 20` ist für kleine und mittlere
Zählwerte genau (Vergleich mit `integrate()`: Fehler unter 1e-4 pro Einheit).
Bei Zählwerten in den Hunderten zusammen mit großem `sigma` wird die feste
Knotenzahl grob; die approximierte Likelihood kann dann mehrere lokale Maxima
haben. In solchen Fällen mit `K = 40` oder `K = 80` nachrechnen und
Log-Likelihood und Schätzwerte vergleichen (siehe `?poisselect`).

**`sigma` und `rho`.** `optim()` arbeitet unrestringiert auf `log(sigma)`
und `atanh(rho)`, sodass `sigma > 0` und `|rho| < 1` konstruktionsbedingt
gelten. Kovarianzmatrix und Standardfehler werden per Delta-Methode auf
`sigma` und `rho` zurücktransformiert; `vcov(fit)` bezieht sich auf die
Originalparameter. Die Zeile für `rho` in `summary()` ist der direkte Test auf
Selektionsverzerrung.

**`predict(type = "response")`** liefert die unbedingte Populationserwartung
`exp(x'beta + sigma^2 / 2)`, nicht das selektierte Mittel; der Faktor
`exp(sigma^2 / 2)` ist der Erwartungswert des lognormalen Fehlers.

**Warnungen.** Nichtkonvergenz von `optim()`, `|rho| > 0.99`, `sigma`
praktisch 0, fehlende Exclusion Restriction, perfekte Trennung in der
Selektionsgleichung sowie eine singuläre oder indefinite Hesse-Matrix werden
als Warnung gemeldet; das Objekt wird trotzdem zurückgegeben
(`fit$converged`, `fit$convergence`), betroffene Standardfehler sind `NA`.

**Fehler.** Falsche Formeln, fehlende Variablen, ungültiges `K`, `start` oder
`control`, ein Indikator außerhalb von {0, 1}, `NA` im Indikator, fehlende
oder nicht-ganzzahlige Zählwerte bei `s = 1`, ausschließlich Nullen,
kollineare Designmatrizen, zu kleine Stichproben, `NA` in Kovariablen und
unpassendes `newdata` werden mit einer Meldung abgefangen, die das Argument,
das Problem und wo möglich die betroffene Variable nennt.

## Tests und Prüfstatus

`devtools::test()` bzw. `R CMD check` führt 340 Assertions in 63
`test_that()`-Blöcken aus (Laufzeit unter 30 s):

| Datei                  | Prüft                                                                                                                                                                                                                                                                     |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `test-check_input.R`   | jede Fehlermeldung einzeln: Formeln, `data`, `K`, `control`, `start`, Indikator, Zählwerte, Kovariablen, Kollinearität, Fallzahl, Exclusion Restriction, perfekte Trennung, Faktoren und Interaktionen, `predict`- und `plot`-Argumente                                   |
| `test-gauss_hermite.R` | Knoten und Gewichte gegen publizierte Werte (K = 2, 3, 5), Summe der Gewichte, exakte Integration von Polynomen bis Grad 2K - 1, Symmetrie, Normalerwartung                                                                                                               |
| `test-loglik.R`        | Log-Likelihood gegen eine naive Referenzimplementierung der Formel, Zerlegung bei `rho = 0`, Grenzfall `sigma -> 0`, Endlichkeit bei `y = 5000`, analytischer Gradient gegen zentrale Differenzen                                                                         |
| `test-estimation.R`    | Parameterrückgewinnung für `rho > 0`, `= 0`, `< 0`, Überlegenheit gegenüber dem naiven GLM, Überdeckung der Wald-Intervalle, Invarianz gegenüber `K`, Startwerte, `control`, Gradient am Optimum, Rand- und Konvergenzwarnungen, Reproduzierbarkeit und Zeilenreihenfolge |
| `test-methods.R`       | Struktur und Werte aller S3-Methoden, Delta-Methode, `AIC`/`BIC`, `predict` mit `newdata` und Faktorstufen, Plots, modell-implizierte Verteilung gegen naive Auswertung der Formel                                                                                        |
| `test-simulate.R`      | Eigenschaften des Generators und des Datensatzes                                                                                                                                                                                                                          |

| Prüfung                                             | Ergebnis                                         |
| --------------------------------------------------- | ------------------------------------------------ |
| `R CMD check --as-cran`                             | 0 Errors, 0 Warnings, 3 umgebungsbedingte Notes¹ |
| `testthat`                                          | 340/340, keine Skips                             |
| `lintr::lint_package()` mit Default-Regeln          | 0 Befunde                                        |
| `cyclocomp` (Maximum über alle Funktionen)          | 6                                                |
| `remotes::install_local(...)` in frische Bibliothek | erfolgreich                                      |

¹ _New submission_ (jedes noch nicht auf CRAN veröffentlichte Paket), _unable
to verify current time_ (Prüfrechner ohne Zeitserver-Zugriff) und _Skipping
checking HTML validation_ (kein aktuelles HTML Tidy installiert). Keine der
drei betrifft den Paketinhalt.

## Hinweis zur KI-Unterstützung

Die Entwicklung erfolgte mit erheblicher Unterstützung von KI-Werkzeugen. Die
Kennzeichnung gemäß Abschnitt 4 der KI-Richtlinie liegt der Abgabe als eigenes
Dokument bei: `KI-Erklaerung.pdf` (Quelle `KI-Erklaerung.md`).
