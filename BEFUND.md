# Befund: Prüfung, Vergleich und Zusammenführung der beiden `poisselect`-Fassungen

Stand: 3. September 2026. Prüfumgebung: R 4.5.1 (aarch64-apple-darwin20),
checkmate 2.3.4, testthat 3.3.2, roxygen2 8.1.0, lintr 3.4.0, kein pandoc,
kein aktuelles HTML Tidy. Alle Zahlen in diesem Dokument stammen aus
tatsächlich ausgeführten Läufen; nichts wurde aus dem Code abgeleitet.

Bezeichnungen: **Eigene** = `eigene_version/poisselect` (byte-identisch mit
`eigene_version/abgabe/poisselect`, per `diff -r` geprüft), **Kollege** =
`Kollege_version/Advanced_r/poisselect`, **Final** = `final_version/poisselect`.
Die Nummern P1 bis P22 verweisen auf das Pflichtenheft in Abschnitt 0.

---

## 0. Pflichtenheft aus der Aufgabenstellung

| Nr | Anforderung | Prüfmittel |
|---|---|---|
| P1 | Approximierte Log-Likelihood exakt in der vorgegebenen Form (`1/sqrt(pi)`, `mu_ik = exp(x'beta + sqrt(2) sigma t_k)`, `eta_ik = (z'gamma + sqrt(2) rho t_k)/sqrt(1 - rho^2)`) | naive Referenz, `integrate()` |
| P2 | Auswertung vollständig auf der Log-Skala, Log-Sum-Exp über die Knoten | `y` bis 5412, extreme Selektionsindizes |
| P3 | Gauß-Hermite-Knoten und -Gewichte zu `exp(-t^2)`, selbst berechnet (Golub-Welsch) | publizierte Werte, `statmod`, Polynome bis Grad `2K-1` |
| P4 | `K` als Argument der Schätzfunktion, Default 20 | Signatur |
| P5 | Maximierung per `optim` über alle Parameter mit Beachtung von `sigma > 0`, `-1 < rho < 1`, geeignete Startwerte | Code, Rückgewinnung, Optimierer-Diagnose |
| P6 | Standardfehler von `beta`, `gamma` aus der Hesse-Matrix am Optimum | Hesse-Vergleich, Überdeckung |
| P7 | S3-Objekt `poisselect` mit `beta, gamma, sigma, rho`, SE, Kovarianzmatrix, Log-Likelihood, Konvergenzstatus, `n`, `n_selected`, `call` | `str(fit)` |
| P8 | `print`: Fallzahlen, Log-Likelihood, Konvergenz, Koeffizienten | Output |
| P9 | `summary`: Tabellen beider Gleichungen mit Estimate, SE, z, zweiseitigem p; `rho`, Log-Likelihood, AIC | Output |
| P10 | Plot 1 nach der vorgegebenen Formel, gemittelt über selektierte Einheiten | Summe der Verteilung, naive Referenz |
| P11 | Plot 2: Log-Likelihood entlang `rho` | Output |
| P12 | `predict` mit `link`, `response` (Default, `exp(x'beta + sigma^2/2)`), `pselect` | Handrechnung |
| P13 | `coef`, `vcov`, `logLik` | Output, `AIC()` |
| P14 | Eigenimplementierung, keine Selektionsmodell-Pakete | `Imports` |
| P15 | Ausführliche Input-Checks mit informativen Fehlermeldungen | Batterie mit 73 Fällen |
| P16 | Dokumentation mit roxygen2 | `roxygenise()` ohne Diff, gerenderte Hilfe |
| P17 | Tests mit testthat inkl. Fehlerbehandlung und Rückgewinnung | Testlauf |
| P18 | CRAN-Checks | `R CMD check --as-cran` |
| P19 | Installierbar per `remotes::install_local(..., dependencies = TRUE, type = "source")` | frische Bibliothek |
| P20 | Zip mit Pseudocode-PDF und README | Dateien |
| P21 | Kurskonventionen: Kommentare, Top-Down, Style Guide, DRY, `lintr` | `lintr`, `cyclocomp`, Lesen |
| P22 | Reproduzierbarkeit | gleicher Seed, Zeilenreihenfolge |

Hinweis zur Prüfumgebung: die Systembibliothek dieses Rechners ist ein
einziges Verzeichnis (Basis- und Zusatzpakete), deshalb konnte die
Abhängigkeitsauflösung von `dependencies = TRUE` nicht gegen eine komplett
leere Bibliothek geprüft werden; installiert wurde jeweils in eine frische,
vorangestellte Bibliothek. Außerdem liegt in der Systembibliothek eine alte
`poisselect`-Installation (Autor „Vorname Nachname"); sie wurde nicht
angefasst, verfälscht aber `lintr`s `object_usage_linter`, wenn das zu
prüfende Paket nicht vorher geladen ist (siehe 1.4).

---

## 1. Funktionsprüfung

### 1.1 Build und Installation (P19)

| | Eigene | Kollege | Final |
|---|---|---|---|
| `R CMD build` | OK | OK | OK |
| `remotes::install_local(tar.gz, dependencies = TRUE, type = "source")` in frische Bibliothek | OK, Exporte `poisselect`, `sim_selection`, `simulate_poisselect` | OK, Exporte `poisselect`, `simulate_poisselect` | OK |

### 1.2 `R CMD check --as-cran` (P18)

| | Eigene | Kollege | Final |
|---|---|---|---|
| ERROR / WARNING | 0 / 0 | 0 / 0 | 0 / 0 |
| NOTE „New submission" | ja (umgebungsbedingt) | ja | ja |
| NOTE „unable to verify current time" | ja (kein Zeitserver-Zugriff) | ja | ja |
| NOTE „Skipping checking HTML validation" (HTML Tidy zu alt) | ja | ja | ja |
| NOTE „README.md cannot be checked without pandoc" | nein (README im `.Rbuildignore`) | **ja** | nein |
| Tests im Check | 15 s OK | 16 s OK | OK |

Alle verbleibenden NOTEs entstehen durch die Prüfumgebung, nicht durch den
Paketinhalt. Die README-NOTE des Kollegen verschwindet nur auf Rechnern mit
pandoc.

### 1.3 Tests (P17)

| | Eigene | Kollege | Final |
|---|---|---|---|
| Dateien / Blöcke / Assertions | 6 / 78 / 296 | 7 / 36 / 131 | 6 / 63 / 340 |
| Fehlschläge / Skips / Warnungen | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| Laufzeit | 14 s | 12 s | unter 30 s |
| naive Referenz der Log-Likelihood | ja | nein (nur Grenzfälle) | ja |
| Delta-Methode geprüft | ja | entfällt | ja |
| analytischer Gradient geprüft | entfällt | entfällt | ja |
| Reproduzierbarkeit | nein | ja | ja (plus Zeilenpermutation) |

### 1.4 Stil (P21)

| | Eigene | Kollege | Final |
|---|---|---|---|
| `lint_package()` mit Paket-Konfiguration | 0 (`.lintr` erlaubt CamelCase) | 0 | 0 |
| `lint_package()` mit Default-Regeln | **5** (`object_name_linter`, alle wegen `K` in 4 R-Dateien und einem Test) | 0 | 0 (eine Zeile mit `# nolint: object_name_linter.` für `K`) |
| `cyclocomp` Maximum | 7 | 8 | 6 |
| Funktionen / R-Dateien | 69 / 20 | 32 / 11 | 57 / 17 |

Achtung bei eigener Nachprüfung: `lintr::lint_package()` in einer frischen
Sitzung meldet auf diesem Rechner 9 Scheinbefunde („no visible global function
definition"), weil `object_usage_linter` die alte installierte Version als
Namespace nimmt. Mit `pkgload::load_all()` vor dem Lint sind es 0.

### 1.5 Dokumentation (P16)

| | Eigene | Kollege | Final |
|---|---|---|---|
| `roxygenise()` erzeugt Diff in `man/` oder `NAMESPACE` | nein | nein | nein |
| alle Exporte dokumentiert, Beispiele laufen | ja | ja | ja |
| gerenderte Hilfe (`Rd2txt`) | **103 Zeilen mit wörtlichen Backticks**, `[stats::optim()]`-Links als Klartext, weil `Roxygen: list(markdown = TRUE)` in `DESCRIPTION` fehlt | sauber | sauber |

### 1.6 Numerische Korrektheit (P1, P2, P3)

Alle drei Fassungen liefern hier identische Zahlen, weil die Likelihood in
allen dreien dieselbe Formel vektorisiert.

| Prüfung | Ergebnis (Eigene = Kollege = Final) |
|---|---|
| Knoten und Gewichte gegen publizierte Werte, `K = 2, 3, 5` | max. Abweichung 4e-15 (Knoten), 2e-15 (Gewichte) |
| gegen `statmod::gauss.quad(K, "hermite")`, `K = 10, 20, 50` | max. 9e-15 |
| `sum(w) = sqrt(pi)` | Abweichung 2e-15 |
| Polynome bis Grad `2K - 1` exakt (relativ zur Summandengröße), `K = 2, 3, 5, 20` | max. 2e-14; Grad `2K` ist wie erwartet nicht exakt |
| Log-Likelihood gegen naive Schleife über Einheiten und Knoten (natürliche Skala, 80 Einheiten, `K = 5` und `20`) | Differenz 6e-14 bzw. 9e-14 |
| Mutationstest: Referenz mit `sigma / sqrt(2)` | Differenz 4.26, der Vergleich hat also Zähne |
| `rho = 0`: Zerlegung in Probit + Poisson-Mischung; `sigma = 1e-8`: Poisson-GLM | Differenz 9e-14 |
| Zählwerte 137 bis 1644 | endlich, Differenz zur Log-Skala-Referenz 9e-13 |
| Selektionsindex mit `gamma = (0, 12, 12)` | endlich, Differenz 2e-13 |
| kompletter Fit mit `y` bis 5412 | läuft durch, `converged = TRUE` |

**Grenzen der Quadratur (Vergleich mit `integrate()`, eine Einheit, `zg = 0.4`):**

| `y` | `sigma` | `rho` | Fehler `K = 5` | Fehler `K = 20` | Fehler `K = 40` |
|---|---|---|---|---|---|
| 0 | 0.3 | 0.5 | 3e-5 | 2e-15 | 3e-15 |
| 3 | 0.3 | 0.5 | 3e-4 | 3e-14 | 3e-15 |
| 25 | 0.3 | 0.5 | 6e-2 | 2e-5 | 4e-10 |
| 3 | 0.8 | 0.5 | 4e-2 | 4e-5 | 3e-8 |
| 25 | 0.8 | 0.5 | 1.5e-1 | **1.7e-1** | 5e-3 |
| 25 | 0.8 | 0.9 | 1.7e-1 | **1.8e-1** | 2e-3 |

Bei `y = 25` und `sigma = 0.8` ist der Fehler der Log-Likelihood-Beiträge mit
`K = 20` schon 0.17 pro Einheit. Das ist keine Eigenschaft einer der
Implementierungen, sondern der vorgeschriebenen festen Quadratur: der
Integrand wird mit wachsendem `y * sigma` immer schmaler. Folge (siehe 1.10):
mehrere lokale Maxima der approximierten Likelihood bei großen Zählwerten
oder großem `sigma`.

### 1.7 Parameterrückgewinnung (P5, P6), 100 Replikationen, `n = 2000`, `K = 20`

Beide Pakete auf identischen Datensätzen (eigener Generator); Bias, RMSE,
mittlerer SE und empirische Überdeckung der 95 %-Wald-Intervalle. Eigene und
Kollege stimmen bis auf die dritte Nachkommastelle überein, Final ist mit der
eigenen Fassung identisch (gleiche Optima, nur schneller). Auszug:

| wahres `rho` | Parameter | Bias | RMSE | mittl. SE | Überdeckung |
|---|---|---|---|---|---|
| 0 | `beta_x1` | 0.002 | 0.030 | 0.032 | 0.95 |
| 0 | `sigma` | -0.002 | 0.027 | 0.026 | 0.94 |
| 0 | `rho` | -0.022 | 0.153 | 0.153 | 0.93 |
| 0.5 | `beta_0` | 0.004 | 0.060 | 0.056 | 0.97 |
| 0.5 | `beta_x1` | -0.007 | 0.034 | 0.029 | 0.87 |
| 0.5 | `rho` | 0.002 | 0.140 | 0.135 | 0.93 |
| 0.8 | `sigma` | -0.007 | 0.033 | 0.028 | 0.90 |
| 0.8 | `rho` | 0.013 | 0.095 | 0.092 | 0.89 |

Konvergenz 300/300 in allen Fassungen, keine Warnungen. Überdeckung über alle
Parameter zwischen 0.87 und 0.98; die leichte Unterdeckung bei `beta_x1`
(`rho = 0.5`) und `rho` (`rho = 0.8`) passt zur Quadraturungenauigkeit bei
`K = 20`. Der Kollege liefert keine SE für `sigma` und `rho` (nicht
gefordert). Laufzeit pro Fit: Eigene 1.0 bis 1.4 s, Kollege 1.7 bis 1.9 s,
Final 0.25 bis 0.29 s.

### 1.8 Fehlerbehandlung (P15), 73 Fälle, identische Batterie

Kriterium: informative Meldung statt Absturz, stillem Unsinn oder kryptischer
Basis-R-Meldung. Nur die Fälle mit Unterschieden oder Schwächen:

| Fall | Eigene | Kollege | Final |
|---|---|---|---|
| 11 `y` beobachtet bei `s = 0` | `message`, Werte ignoriert | still ignoriert | `message` |
| 14 perfekte Trennung (`s = 1{z1 > 0}`) | `glm.fit`-Warnungen sickern durch, `gamma_z1 = 408` mit SE 9990, Hesse-Warnung | **still**, `gamma_z1 = 408`, SE 5090 | Warnung „separate ... perfectly" plus Hesse-Warnung |
| 18 `K = 1` | läuft, Warnung singuläre Hesse, SE(`sigma`) = NA | Fehler „>= 2" | Fehler „>= 2" |
| 26 `K = 1000` | Fehler (max 200) | Fehler (max 500) | Fehler (max 200) |
| 27 `data` fehlt | Basis-R „argument "data" is missing" | gleich | eigene Meldung |
| 35 `.` in der Formel | eigene Meldung | „variables ... missing from data: ." | eigene Meldung |
| 42 `n = 12` | Hesse-Warnung, 2 SE = NA | **still**, `sigma = 0.009`, `rho = 0.87` | Hesse-Warnung |
| 43 nur 2 selektierte Einheiten | „observed for only 2 unit(s)" | „collinear columns" (irreführend) | „observed for only 2 unit(s)" |
| 47 `NA` in `x` nur bei `s = 0` | **Fehler**, Meldung behauptet fälschlich, nicht-selektierte Einheiten seien Teil beider Gleichungen | läuft | läuft |
| 53 `control = list(bogus = 1)` | nur `optim`-Warnung „unknown names" | Fehler mit erlaubter Namensliste | Fehler mit erlaubter Namensliste |
| 62 alle selektierten `y = 0` | **still**, `beta_0 = -27.3`, SE 64600 | Hesse-Warnung, alle SE NA | Fehler „outcome is 0 for every selected unit" |
| 66 `newdata` mit `NA` | Meldung nennt Spalte `x1TRUE` | gleich | Meldung nennt Variable `x1` |
| 69 neues Faktorlevel in `newdata` | Basis-R „factor grp has new level" | gleich | gleich (akzeptabel) |

Alle übrigen 60 Fälle: beide Ursprungsfassungen und Final geben eine
informative Meldung oder laufen wie vorgesehen durch (Faktoren,
Interaktionen, `y ~ 1`, logischer Indikator, tibble, `K = 150`,
`sigma`-Randfall mit `sigma = 0.075` und SE 0.14).

### 1.9 Reproduzierbarkeit (P22)

| | Eigene | Kollege | Final |
|---|---|---|---|
| `set.seed` im Paketcode | nur in Beispielen | nur in Beispielen | nur in Beispielen |
| gleiche Daten zweimal | identisch | identisch | identisch |
| Zeilenpermutation, max. Koeffizientenabweichung | 5e-11 | 2e-10 | 1e-9 |
| Fit verändert `.Random.seed` | nein | nein | nein |

### 1.10 Kreuzvergleich (P5) auf identischen Datensätzen

Log-Likelihood am jeweils berichteten Optimum, ausgewertet mit derselben
(verifizierten) Likelihood-Funktion; „Lücke" = Abstand zum besten Optimum, das
mit BFGS-Neustarts, Nelder-Mead und analytischem Gradienten gefunden wurde.

| Szenario | Eigene | Kollege | Lücke Eigene | Lücke Kollege | Final |
|---|---|---|---|---|---|
| A `rho = 0.5`, `n = 1000` | -1743.7749 | -1743.7749 | 0 | 0 | -1743.7749 |
| B `rho = 0`, `n = 1000` | -1696.2432 | -1696.2432 | 0 | 0 | -1696.2432 |
| C `rho = 0.8`, `n = 2000` | -3606.676 | -3607.482 | 0 | 0.81 | -3606.676 |
| D `rho = -0.6`, `sigma = 0.9`, `n = 800` | -1332.649 | -1332.172 | 0.48 | 0 | -1332.649 |
| E `beta_0 = 6` (Zählwerte bis ~1600) | -1456.576 | -1480.534 | 0.0001 | **24.0** | -1456.576 |
| F `beta_0 = 7.5` (Zählwerte bis 5412) | -1142.942 | -1096.135 | **46.8** | 0 | -1110.179 |
| G `n = 150` | -253.571 | -253.426 | 0.22 | 0.07 | -253.827 (Rand `rho = 1`, Warnung) |

Zwei getrennte Ursachen, beide mit ausgeführten Läufen belegt:

1. **Vorzeitiger Stopp durch den numerischen Gradienten** (E, F): mit
   analytischem Gradienten ist die Lücke in E und F exakt 0, sowohl von der
   eigenen als auch von der Kollegen-Startlösung aus. Der numerische Gradient
   von `optim` ist bei Intercepts um 7 und Zählwerten in den Tausenden zu
   ungenau.
2. **Mehrere lokale Maxima der groben Quadratur** (C, D, G): von der eigenen
   Startlösung, der Kollegen-Startlösung (`sigma = 1`, `rho = 0`) und vom
   wahren Parameter aus landet BFGS mit `K = 20` an verschiedenen Punkten.
   Mit `K = 80` konvergieren alle drei Startpunkte in C, D und G auf dieselbe
   Lösung. Das ist die Grenze der vorgeschriebenen Methode (siehe 1.6), keine
   Implementierungsschwäche einer Seite; welche Seite „recht" hat, entscheidet
   hier die Startlösung, nicht der Code.

Für F liefert Final mit `K = 80` bzw. `150` `sigma = 0.25`, `rho = 0.27` bzw.
`0.20` (wahr 0.3, 0.4) und Log-Likelihood -1039 statt -1110: die in `?poisselect`
dokumentierte Empfehlung, bei großen Zählwerten `K` zu erhöhen, ist damit belegt.

**Standardfehler (P6), Hesse-Verfahren im Vergleich** (Referenz: `numDeriv`
Richardson mit `r = 8`; relative Abweichung der SE, drei Datensätze):

| Verfahren | max. rel. Abweichung | Laufzeit |
|---|---|---|
| `numDeriv::hessian()` Default (Eigene) | 6e-7 bis 6e-5 | 0.10 bis 0.65 s |
| `optimHess()`, `ndeps = 1e-3` (Kollege, über `optim(hessian = TRUE)`) | 5e-6 bis 9e-5 | 0.08 bis 0.58 s |
| `optimHess()` mit analytischem Gradienten (Final) | Delta-Methode-Test auf 1e-8 genau gegen dieselbe Formel | vergleichbar |

Die Genauigkeit von `optimHess` ist für Standardfehler mehr als ausreichend,
`numDeriv` als Abhängigkeit deshalb verzichtbar.

### 1.11 Befundtabellen je Version

**Eigene Version**

| Punkt | Status | Beleg | Schwere |
|---|---|---|---|
| P19 Installation | OK | 1.1 | |
| P18 Check | OK (3 Umgebungs-NOTEs) | 1.2 | |
| P17 Tests | OK, 296 Assertions | 1.3 | |
| P21 lintr Default | Stil: 5 Befunde `K` | 1.4 | Minor |
| P16 gerenderte Hilfe | Bug: Markdown nicht aktiviert, 103 Zeilen mit Backticks | 1.5 | Minor (Doku), fällt beim Lesen der Hilfe sofort auf |
| P1, P2, P3 Numerik | OK | 1.6 | |
| P5 Optimierer | Bug: vorzeitiger Stopp bei großen Zählwerten (Lücke 47) | 1.10 | Major, aber nur außerhalb des `K = 20`-Gültigkeitsbereichs |
| P6, P7 Rückgewinnung, Überdeckung | OK | 1.7 | |
| P15 Fall 47 | Spec-Abweichung: `NA` in `x` bei `s = 0` abgelehnt, Meldung inhaltlich falsch | 1.8 | Major (falsche Datensätze werden abgelehnt) |
| P15 Fall 62 | Bug: alle `y = 0` läuft still mit `beta_0 = -27` | 1.8 | Major |
| P15 Fall 14 | Stil: `glm.fit`-Warnungen sickern durch | 1.8 | Minor |
| P15 Fall 53 | fehlende Prüfung der `control`-Namen | 1.8 | Minor |
| DESCRIPTION | Bug: Autor „Vorname Nachname" | Datei | Major (Paket-Metadaten) |
| P22 | OK | 1.9 | |

**Version Kollege**

| Punkt | Status | Beleg | Schwere |
|---|---|---|---|
| P19 Installation | OK | 1.1 | |
| P18 Check | OK, aber README-NOTE ohne pandoc | 1.2 | Minor |
| P17 Tests | OK, 131 Assertions; keine naive Referenz der Likelihood, kein Test der Faktorstufen in `predict` | 1.3 | Major (Testabdeckung, „ausführlich (!)" gefordert) |
| P21 lintr Default | OK | 1.4 | |
| P16 Doku | OK | 1.5 | |
| P4 Argumentname | Spec-Abweichung: `k` statt `K` | Signatur | Minor, aber ein Tippfehler-Risiko für die Korrektur |
| P1, P2, P3 Numerik | OK | 1.6 | |
| P5 Optimierer | Bug: vorzeitiger Stopp bei großen Zählwerten (Lücke 24) | 1.10 | Major (wie oben) |
| P6, P7 | OK; keine SE für `sigma`, `rho` (nicht gefordert) | 1.7 | Kosmetik |
| P15 Fall 14 | Bug: perfekte Trennung still, `gamma = 408` | 1.8 | Major |
| P15 Fall 42 | fehlende Warnung bei degeneriertem Fit | 1.8 | Minor |
| P15 Fall 35 | kryptische Meldung bei `.` | 1.8 | Minor |
| P15 Fall 11 | `y` bei `s = 0` still ignoriert | 1.8 | Minor |
| P15 Fall 43 | irreführende Meldung „collinear" bei zu wenigen Einheiten | 1.8 | Minor |
| README | behauptet 0 NOTEs; stimmt nur mit pandoc | Datei | Kosmetik |
| Pseudocode | Reihenfolge „bottom-up, testgetrieben" statt Top-Down der Vorlesung | Datei | Minor |
| P22 | OK | 1.9 | |

---

## 2. Vergleich und Fehleranalyse

### 2.1 Feature-Matrix

| Anforderung | Eigene | Kollege | Final |
|---|---|---|---|
| P1 Log-Likelihood | `compute_loglik()` + `compute_node_pieces()` | `compute_poisselect_loglik()` + `compute_node_matrices()` | wie Eigene, vereinfachtes Log-Sum-Exp |
| P2 Log-Skala | `row_log_sum_exp()` mit `max.col` | `apply(a, 1, max)` | `apply(a, 1, max)` |
| P3 Golub-Welsch | Index-Matrizen | Index-Matrizen | obere Nebendiagonale + Transponierte |
| P4 `K` | `K`, 1 bis 200, `.lintr` | `k`, 2 bis 500 | `K`, 2 bis 200, `# nolint` |
| P5 Optimierung | BFGS, numerischer Gradient, `reltol 1e-10` | BFGS, numerischer Gradient, Default-`reltol` | BFGS, **analytischer Gradient**, `reltol 1e-10` |
| P5 Restriktionen | `log(sigma)`, `atanh(rho)` | gleich | gleich |
| P5 Startwerte | Probit, Poisson-GLM, Momentenschätzer, Intercept-Korrektur | Probit, Poisson-GLM, `sigma = 1`, `tryCatch` | Eigene + `suppressWarnings`, Trennungs-Check |
| P6 Hesse | `numDeriv` | `optim(hessian = TRUE)` | `optimHess()` des analytischen Gradienten |
| P6 SE-Rücktransformation | Delta-Methode, volle `vcov` | nur `(beta, gamma)`-Block | Delta-Methode, volle `vcov` |
| P7 Objekt | 17 Elemente, SE für `sigma`, `rho` | 12 Elemente | 18 Elemente |
| P8 `print` | vollständig | ohne AIC | vollständig |
| P9 `summary` | drei Tabellen, Richtung der Verzerrung | zwei Tabellen | drei Tabellen |
| P10 Plot 1 | gruppierte Balken, Log-Skala, 99 %-Kappung | Balken + Punkte, natürliche Skala | wie Eigene |
| P11 Plot 2 | Gitter ±0.99 mit Schätzwert | `rho_grid`-Argument | wie Eigene |
| P12 `predict` | drei Typen, `newdata` | gleich | gleich, bessere `NA`-Meldung |
| P13 Extraktoren | `coef(which)`, `vcov`, `logLik` mit `df`/`nobs` | `coef`, `vcov`, `logLik` | wie Eigene |
| P14 Abhängigkeiten | checkmate, numDeriv | checkmate | checkmate |
| P15 Checks | 16 Checkfunktionen | 6 Checkfunktionen | 12 Checkfunktionen |
| P16 Doku | roxygen, Markdown nicht aktiviert | roxygen, Markdown | roxygen, Markdown |
| P17 Tests | 296 | 131 | 340 |
| P20 Abgabe | Zip mit README, Pseudocode-PDF (Python), KI-Erklärung | Zip mit README, Pseudocode-PDF (Quarto) | README, Pseudocode (Quarto), KI-Erklärung, Befund, Merge-Entscheidungen |
| Datensatz | `sim_selection` | keiner | `sim_selection` |
| Simulator | festes Design | flexibles Design | festes Design |

### 2.2 Lesbarkeit, Top-Down, DRY, Code Smells

- **Eigene:** sauberer Top-Down-Aufbau (Hauptfunktion ruft neun benannte
  Schritte), aber Überfragmentierung: 69 Funktionen, darunter Einzeiler wie
  `add_log_weights()`, `replace_nonfinite()`, `split_standard_errors()`,
  `check_control_names()`. Zwei „clevere" Stellen ohne Nutzen
  (`max.col`-Trick, `rep(..., each = )` für die Gewichte). Roxygen-Blöcke mit
  redundanten `@return invisible(TRUE)`-Sätzen an jedem internen Helfer.
- **Kollege:** kompakter, gute Kommentare mit Begründung, aber
  `check_inputs.R` mischt Prüfung und Datenaufbau (`build_model_data()` ruft
  die Checks selbst), `poisselect.R` enthält vier Helfer, die in eigene
  Dateien gehören (Zielfunktion, Kovarianz, Konstruktor, Warnungen); `::`
  überall statt `@importFrom` (zulässig, aber lauter). `cyclocomp` bis 8.
- **Final:** Top-Down-Skelett der eigenen Fassung, Helferzahl auf 57 reduziert,
  Kommentare auf „warum" beschränkt, `cyclocomp` maximal 6, keine Variable
  überschattet einen Basisnamen (die Parameterliste heißt `parameters`, nicht
  `par`).

### 2.3 Fehlerliste mit Schwere

| # | Fassung | Befund | Schwere | Behoben in Final |
|---|---|---|---|---|
| 1 | beide | vorzeitiger Stopp des Optimierers bei großen Zählwerten (numerischer Gradient) | Blocker für solche Daten, sonst unsichtbar | ja, analytischer Gradient |
| 2 | Eigene | Autor-Platzhalter in `DESCRIPTION` | Major | ja |
| 3 | Eigene | `NA` in Outcome-Kovariablen nicht-selektierter Einheiten wird abgelehnt, Meldung falsch | Major | ja |
| 4 | Eigene | alle `y = 0` läuft still durch | Major | ja, Fehler |
| 5 | Kollege | perfekte Trennung still | Major | ja, Warnung |
| 6 | Kollege | Testabdeckung: keine Referenzimplementierung, keine Delta-/Faktortests | Major | ja, Vereinigung |
| 7 | Eigene | gerenderte Hilfe zeigt Markdown-Syntax | Minor | ja |
| 8 | Eigene | `lintr` Default: 5 Treffer | Minor | ja |
| 9 | Kollege | README-NOTE ohne pandoc | Minor | ja |
| 10 | Kollege | `k` statt `K` | Minor | ja |
| 11 | Kollege | kryptische `.`-Meldung, stille Fälle 11 und 42, irreführende Meldung 43 | Minor | ja |
| 12 | Eigene | `glm.fit`-Warnungen sickern durch; `control`-Namen ungeprüft | Minor | ja |
| 13 | beide | `newdata`-`NA`-Meldung nennt `x1TRUE` | Kosmetik | ja |
| 14 | beide | Lizenz uneinheitlich (GPL-3 vs. MIT) | Kosmetik | GPL-3 |
| 15 | beide | `K = 20` bei großen Zählwerten zu grob, lokale Maxima | inhärent, dokumentiert | Doku-Abschnitt, Empfehlung `K = 40/80` |

Die Zeilen „Meine Version / Version Kollege / besser ist / weil" stehen je
Teilproblem in `MERGE_ENTSCHEIDUNGEN.md`.

### 2.4 Die vier vorgegebenen Punkte

- **Autor-Platzhalter:** ersetzt durch Sebastian Jung (aut, cre) und Finn
  Metzler (aut). Schreibweise und Maintainer-Rolle bitte prüfen.
- **Lizenz:** GPL-3, weil das Kurs-Beispiel `Beispiel_Paket.R` mit
  `use_gpl3_license()` arbeitet; `LICENSE.md` liegt bei und ist im
  `.Rbuildignore` (Standardvorgehen von `usethis`, keine Check-NOTE).
- **numDeriv:** entfällt. Beleg: `optimHess` weicht in den SE relativ höchstens
  9e-5 von der Richardson-Referenz ab, mit analytischem Gradienten noch
  weniger; die Monte-Carlo-Überdeckung ist mit beiden Verfahren identisch.
- **Autorenschaft:** beide Personen in `DESCRIPTION`, README, Pseudocode und
  KI-Erklärung; der Anteil der Werkzeuge ist in der KI-Erklärung ausgewiesen.

---

## 3. Abnahme der zusammengeführten Fassung

| Kriterium | Nachweis |
|---|---|
| `R CMD check --as-cran` | 0 ERROR, 0 WARNING, 3 NOTEs (New submission; unable to verify current time; HTML Tidy), alle umgebungsbedingt und in 1.2 begründet |
| Tests | 6 Dateien, 63 Blöcke, 340 Assertions, 0 Fehler, 0 Skips, 0 Warnungen |
| `lintr` Default-Regeln | 0 Befunde (mit geladenem Paket, siehe 1.4); keine `.lintr` |
| `roxygenise()` | kein Diff |
| Installation aus dem `tar.gz` in frische Bibliothek | erfolgreich, `library()` und Exporte geprüft |
| Ende-zu-Ende-Skript `beispiel_end_to_end.R` aus der frischen Bibliothek | läuft komplett durch: Fit, `print`, `summary`, beide Plots (PDF), drei `predict`-Typen mit und ohne `newdata`, `coef`, `vcov`, `logLik`, `AIC`, `BIC`, Vergleich mit naivem GLM, `K`-Sensitivität |
| Fehlerbatterie (73 Fälle) | alle Fälle informativ, siehe 1.8 |
| Monte-Carlo | wie 1.7, 300/300 konvergiert, Laufzeit 0.25 s pro Fit |
| Numerik | wie 1.6, identisch zu den Ursprungsfassungen |
