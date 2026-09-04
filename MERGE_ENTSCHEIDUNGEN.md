# Merge-Entscheidungen für `final_version/poisselect`

Grundlage: die Aufgabenstellung (`Hausarbeit_poisselect.pdf`) und die
ausgeführten Prüfungen aus `BEFUND.md`. „Eigene" = `eigene_version/poisselect`
(Sebastian), „Kollege" = `Kollege_version/Advanced_r/poisselect` (Finn).
„Neu" = im Zuge der Zusammenführung von Fable 5.1 geschrieben, weil keine der
beiden Fassungen das Problem gelöst hatte.

## Struktur und Schnittstelle

| Punkt | Eigene | Kollege | Übernommen | Warum |
|---|---|---|---|---|
| Aufruf | `poisselect(outcome, selection, data, K, start, control)` | gleich, aber `k` | Eigene (`K`) | Die Aufgabenstellung nennt das Argument `K`; eine Korrektorin, die `K = 10` tippt, bekommt beim Kollegen „unused argument". Die einzige Stilabweichung ist per `# nolint` an genau einer Zeile dokumentiert, keine `.lintr` nötig. |
| Zulässiges `K` | 1 bis 200 | 2 bis 500 | 2 bis 200 | Mit `K = 1` sitzt der einzige Knoten bei 0, `sigma` fällt aus der Likelihood heraus (Befund: singuläre Hesse-Matrix). Ein klarer Fehler ist besser als ein Fit mit `NA`-Standardfehlern. Oberhalb von ~50 Knoten bringt mehr `K` nichts mehr (integrate()-Vergleich). |
| Dateiaufteilung | 20 Dateien, 69 Funktionen | 11 Dateien, 32 Funktionen | 17 Dateien, 57 Funktionen | Eine Datei je Teilaufgabe des Top-Down-Entwurfs (Eigene), aber ohne die Kleinsthelfer der eigenen Fassung (`add_log_weights`, `describe_convergence` bleibt, `replace_nonfinite`, `split_standard_errors`, `check_control_names`, `compute_implied_probability`, `choose_max_count` wurden eingeschmolzen). |
| Ergebnisobjekt | `coefficients$outcome/selection`, `sigma`, `rho` oben, SE auch für `sigma`/`rho`, volle `vcov` | `coefficients` mit vier Einträgen, `vcov` nur für `(beta, gamma)` | Eigene | Die Delta-Methode liefert SE für `sigma` und `rho` praktisch gratis; die Zeile für `rho` in `summary()` ist der direkte Test auf Selektionsverzerrung. Überdeckung der `rho`-Intervalle im Monte-Carlo: 0.89 bis 0.93. |
| Namen in `coef()`/`vcov()` | `outcome_x1` | `outcome:x1` | Eigene | `:` liest sich in R wie ein Interaktionsterm. |
| `coef(which = )` | ja | nein | Eigene | Gleichungsweise Koeffizienten unter ihren Variablennamen sind in Tests und Beispielen bequemer. |
| Beispieldatensatz `sim_selection` | ja | nein | Eigene | Beispiele laufen ohne `set.seed`-Vorspann und sind exakt reproduzierbar. |

## Numerik

| Punkt | Eigene | Kollege | Übernommen | Warum |
|---|---|---|---|---|
| Gauß-Hermite (Golub-Welsch) | Index-Matrizen für beide Nebendiagonalen | gleich | Neu (vereinfacht) | Obere Nebendiagonale füllen und `jacobi + t(jacobi)`: drei Zeilen, ohne Indexakrobatik. Ergebnis identisch (1e-15 zu publizierten Werten und zu `statmod`). |
| Log-Likelihood | `compute_node_pieces()` + `row_log_sum_exp()` mit `max.col`-Trick und `-Inf`-Sonderbehandlung | `sweep()` + `apply(a, 1, max)` | Eigene Struktur, Kollegen-Vereinfachung | `apply(a, 1, max)` ist das, was die Vorlesung nahelegt; der `max.col`-Trick war nur clever. Die `-Inf`-Sonderbehandlung ist überflüssig: ein `NaN` läuft in die 1e10-Strafe der Zielfunktion. Beide Fassungen stimmten mit der naiven Referenz auf 1e-13 überein. |
| Gewichte addieren | `rep(log_weight, each = n)` (spaltenweise Speicherung ausgenutzt) | `sweep()` | Neu: `matrix(..., byrow = TRUE)` | Am leichtesten zu erklären. |
| Gradient | numerisch (optim-intern) | numerisch | **Neu: analytisch** | Befund 1.10: mit numerischem Gradienten hielten beide Optimierer bei großen Zählwerten vorzeitig an (Lücke zum Optimum 33 bzw. 24 Log-Likelihood-Einheiten). Mit analytischem Gradienten: Lücke 0, 3- bis 4-mal schneller. Gegen zentrale Differenzen geprüft (Abweichung 3e-7). |
| Optimierer | BFGS, `reltol = 1e-10`, `maxit = 1000` | BFGS, Default-`reltol`, `maxit = 500` | Eigene Toleranzen | Die Hesse-Matrix wird am zurückgegebenen Punkt genommen; eine straffe Toleranz ist bei vektorisierter Likelihood billig. |
| Startwerte | Probit + Poisson-GLM, Momentenschätzer für `sigma`, Intercept-Korrektur | Probit + Poisson-GLM mit `tryCatch`-Fallback, `sigma = 1` | Eigene Werte + Kollegen-Absicherung | Der Momentenschätzer startet nahe am wahren `sigma`. `suppressWarnings(glm.fit())` und Ersetzen nicht-endlicher Koeffizienten durch 0 stammen von der Fallback-Idee des Kollegen; die durchsickernden `glm.fit`-Warnungen (Befund 1.8, Fall 14) sind damit weg. |
| Hesse-Matrix | `numDeriv::hessian()` (Richardson) | `optim(hessian = TRUE)` | **Kollege, ohne numDeriv** | Befund 1.6: `optimHess` weicht relativ < 1e-4 von der Richardson-Referenz ab, mit analytischem Gradienten sogar < 1e-6. Eine Abhängigkeit weniger auf dem frischen R der Korrektur. |
| Rücktransformation der SE | Delta-Methode, diagonale Jacobi | nur der invariante `(beta, gamma)`-Block | Eigene | Liefert SE für `sigma` und `rho` und die volle `vcov`. Test: `vcov` = skalierte Inverse der Hesse-Matrix. |
| Singuläre / indefinite Hesse-Matrix | getrennte Warnungen, `NA` nur für betroffene Varianzen | eine Warnung, alles `NA` | Eigene | Feiner: bei `K`-Problemen oder Randlösungen bleiben die verwertbaren SE erhalten. |
| Randlösungen | keine Warnung | Warnung bei `|rho| > 0.99`, `sigma < 1e-4` | Kollege | Befund 1.10, Szenario G: `rho` läuft bei kleinem `n` gegen 1. Der Nutzer muss das sehen. |

## Input-Checks

| Punkt | Eigene | Kollege | Übernommen | Warum |
|---|---|---|---|---|
| Zweistufig (Argumente, dann Daten) | ja | ja | beide | Vorlesung „fail fast, fail appropriately". |
| `data` fehlt | Basis-R-Meldung | Basis-R-Meldung | Neu: eigene Meldung | Kriterium „informative Meldung statt kryptischer Basis-R-Meldung". |
| `.` in der Formel | eigene Meldung | „variable `.` missing" | Eigene | Kollegen-Meldung war kryptisch (Befund 1.8, Fall 35). |
| `control` | benannte Liste | Namen ⊂ optim-Controls | Kollege | Tippfehler wie `maxiter` werden abgefangen; bei der eigenen Fassung nur eine optim-Warnung (Fall 53). |
| Selektionsindikator | logisch erlaubt, Faktor mit Hinweis abgelehnt | logisch erlaubt, Faktor generisch abgelehnt | Eigene | Die Meldung nennt den Grund (Level-Reihenfolge) und die Lösung. |
| `y` beobachtet bei `s = 0` | `message()` | still ignoriert | Eigene | Simulierte Daten tragen `y` oft vollständig; stilles Ignorieren verschleiert Vertauschungen von `s`. |
| `NA` in `x` bei nicht-selektierten Einheiten | Fehler | erlaubt | **Kollege** | Die Outcome-Kovariablen nicht-selektierter Einheiten gehen nirgends in die Likelihood ein; reale Daten (Lohngleichung bei Nicht-Erwerbstätigen) haben genau dieses Muster. Der Rangcheck läuft dann auf `x_selected`. |
| Alle `y = 0` | still (Intercept -27, SE 64600) | Hesse-Warnung | Neu: Fehler | Der ML-Schätzer existiert nicht; ein Fit ist Unsinn (Fall 62). |
| Perfekte Trennung | `glm.fit`-Warnungen sickern durch | still (`gamma = 408`, SE 9990) | Neu: Warnung | Die Probit-Startschätzung hat dann Devianz numerisch 0; das wird geprüft und benannt (Fall 14). |
| Exclusion Restriction | Warnung; Intercept zählt als exklusiv | Warnung; Intercept zählt nicht | Kollege | `y ~ 0 + x1`, `s ~ x1` hat keine echte Exclusion Restriction. |
| Fallzahl | `n > npar`, `n_sel > n_beta` | `n_sel >= p + 2`, `n >= p+q+2` | Mischung | `n > npar` und `n_sel > n_beta + 1` (beta und sigma leben auf den selektierten Einheiten); vor dem Rangcheck, damit bei 2 selektierten Einheiten nicht „rank deficient" gemeldet wird (Fall 43 beim Kollegen). |
| `newdata` mit `NA` | Meldung nennt `x1TRUE` | ebenso | Neu | Prüfung auf dem Modellframe, bevor `model.matrix` aus `NA` eine Dummy-Spalte macht (Fall 66). |

## Methoden

| Punkt | Eigene | Kollege | Übernommen | Warum |
|---|---|---|---|---|
| `print` | Fallzahlen, LL, AIC, K, Konvergenz, drei Blöcke | ähnlich, ohne AIC | Eigene | vollständiger. |
| `summary` | zwei Tabellen + Tabelle `sigma`/`rho`, Richtung der Verzerrung im Klartext | zwei Tabellen, `sigma`/`rho` ohne SE | Eigene | Die Aufgabenstellung verlangt `rho`, LL und AIC; die SE-Zeile für `rho` ist der Selektionstest. |
| Plot 1 | gruppiertes Balkendiagramm, Log-Skala, Kappung beim 99 %-Quantil | Balken + Punkte, natürliche Skala | Eigene | Log-Skala wie in der Likelihood (DRY über `compute_node_pieces`), Kappung verhindert, dass ein Ausreißer alles zusammenquetscht. Mit naiver Auswertung der Formel der Aufgabenstellung getestet. |
| Plot 2 | Gitter ±0.99 inkl. Schätzwert, Markierung von Schätzwert und `rho = 0` | `rho_grid`-Argument, Punkt am Gittermaximum | Eigene | Der Schätzwert liegt garantiert auf dem Gitter, das Maximum der Kurve ist damit der Fit selbst (getestet). |
| `predict` | drei Typen, `newdata` mit `terms`/`xlevels` | gleich | beide | identische Logik. |
| `simulate_poisselect` | festes Design `x1, x2, z1`, `y_complete` | flexibles `p`, `q` mit `w`-Variablen | Eigene | Einfacher zu erklären; `y_complete` macht die Verzerrung sichtbar. |

## Paketmetadaten

| Punkt | Entscheidung | Warum |
|---|---|---|
| Autoren | Sebastian Jung (aut, cre), Finn Metzler (aut) | Beide Beteiligten; der Platzhalter „Vorname Nachname" ist ersetzt. **Bitte prüfen**, ob die Schreibweisen und E-Mail-Adressen stimmen und wer Maintainer sein soll. |
| Lizenz | GPL-3, `LICENSE.md` mit dem GPL-3-Text (im `.Rbuildignore`, wie `usethis::use_gpl3_license()` es anlegt) | Das Kurs-Beispiel `Beispiel_Paket.R` verwendet `use_gpl3_license()`; GPL-3 braucht in `DESCRIPTION` keine Zusatzdatei, also keine Check-NOTE. |
| `Roxygen: list(markdown = TRUE)` | gesetzt | Ohne diesen Eintrag standen in der gerenderten Hilfe der eigenen Fassung 103 Zeilen mit wörtlichen Backticks (Befund 1.5). |
| `README.md` im `.Rbuildignore` | ja | Ohne pandoc erzeugt `R CMD check` sonst eine NOTE (Befund beim Kollegen). Die README liegt auf oberster Ebene der Abgabe. |
| Abhängigkeiten | `checkmate (>= 2.1.0)`, `stats`, `graphics`, `utils` | `numDeriv` entfällt (siehe Hesse-Matrix); `assert_formula()` gibt es ab checkmate 2.1.0. |

## Tests

Vereinigung beider Suiten plus die Lücken aus Befund 1.8:

- aus der eigenen Fassung: naive Referenz der Log-Likelihood, `rho = 0`-Zerlegung, `sigma -> 0`-Grenzfall, Delta-Methode, Profil-Peak, Count-Verteilung, Faktorstufen in `predict`, alle Argumentchecks, Wald-Überdeckung;
- aus der Kollegen-Fassung: Fit mit Zählwerten > 150 ohne Warnung, Reproduzierbarkeit (gleiche Daten zweimal), Randlösungs-Warnungen, `control`-Namen, `NA` in `x` bei nicht-selektierten Zeilen erlaubt;
- neu: analytischer Gradient gegen zentrale Differenzen, Gradient ≈ 0 am Optimum, Zeilenpermutation, perfekte Trennung, `y` alle 0, `data` fehlt, `K = 1`, Intercept zählt nicht als Exclusion Restriction, modell-implizierte Verteilung gegen naive Auswertung der Formel, Kappung des Supports.

Ergebnis: 6 Dateien, 63 Blöcke, 340 Assertions, 0 Fehler, 0 Skips, 0 Warnungen.
