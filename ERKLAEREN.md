# Das musst du erklären können

Die Stellen im finalen Paket, die im Kolloquium am ehesten hinterfragt werden,
mit einer Antwort, die du so verwenden kannst. Datei und Funktion stehen
jeweils dabei.

1. **Warum `sqrt(2) * t_k` und `1/sqrt(pi)`?** (`loglik.R`, Kopfkommentar)
   Das Integral läuft gegen die Standardnormaldichte `phi(v)`, Gauß-Hermite
   ist aber für die Gewichtsfunktion `exp(-t^2)` definiert. Die Substitution
   `v = sqrt(2) t` macht aus `exp(-v^2/2)` genau `exp(-t^2)`; übrig bleibt
   der konstante Faktor `dv / sqrt(2 pi) = sqrt(2) dt / sqrt(2 pi) = dt /
   sqrt(pi)`. Deshalb steht in jeder Rate `sqrt(2) sigma t_k` und vor der
   Summe `1/sqrt(pi)`. Ein Test mit `sigma / sqrt(2)` statt `sigma` zeigt,
   dass die naive Referenz den Fehler sofort bemerken würde.

2. **Woher kommt `eta = (z'gamma + sqrt(2) rho t_k) / sqrt(1 - rho^2)`?**
   (`compute_node_pieces()`) Gegeben der standardisierte Outcome-Fehler `v`
   ist der Selektionsfehler `u` bedingt normalverteilt mit Erwartungswert
   `rho v` und Varianz `1 - rho^2`. Die Wahrscheinlichkeit `u > -z'gamma` ist
   deshalb `Phi((z'gamma + rho v) / sqrt(1 - rho^2))`; am Knoten wird `v`
   durch `sqrt(2) t_k` ersetzt.

3. **Was macht der Log-Sum-Exp-Trick und warum ist er nötig?**
   (`log_sum_exp_rows()`) Für `y = 500` ist `dpois(y, mu)` etwa 1e-300 und
   unterläuft beim Multiplizieren mit den Gewichten zu 0, `log(0) = -Inf`.
   Ich bilde deshalb jeden Summanden auf der Log-Skala, ziehe das
   Zeilenmaximum `A_i` ab, so dass der größte Summand `exp(0) = 1` ist,
   summiere und addiere `A_i` wieder. Mathematisch ist das identisch,
   numerisch bleibt alles im darstellbaren Bereich.

4. **Warum `pnorm(x, lower.tail = FALSE, log.p = TRUE)` statt `log(1 -
   pnorm(x))`?** (`compute_loglik()`) Für `z'gamma = 10` ist `pnorm(10)` in
   doppelter Genauigkeit exakt 1, also wäre `1 - pnorm` gleich 0 und der
   Logarithmus `-Inf`. `pnorm` kann den oberen Schwanz und dessen Logarithmus
   direkt berechnen, ohne die Auslöschung.

5. **Wie funktioniert Golub-Welsch?** (`build_gauss_hermite()`) Die monischen
   Hermite-Polynome erfüllen `pi_{k+1}(t) = t pi_k(t) - (k/2) pi_{k-1}(t)`.
   Die symmetrische Tridiagonalmatrix mit Nullen auf der Diagonale und
   `sqrt(k/2)` daneben hat als Eigenwerte die Nullstellen von `pi_K`, das sind
   die Knoten. Die Gewichte sind `mu_0` mal das Quadrat der ersten Komponente
   des normierten Eigenvektors, wobei `mu_0 = Integral exp(-t^2) dt =
   sqrt(pi)`. Geprüft gegen publizierte Werte für `K = 2, 3, 5` auf 1e-15 und
   über die exakte Integration aller Polynome bis Grad `2K - 1`.

6. **Warum `log(sigma)` und `atanh(rho)`?** (`parameters.R`) `optim` mit
   BFGS kennt keine Schranken. `exp` und `tanh` bilden die ganze reelle Achse
   auf `sigma > 0` und `-1 < rho < 1` ab, die Restriktionen gelten also
   konstruktionsbedingt. `atanh` ist Fishers z-Transformation, der
   Standardweg für Korrelationen.

7. **Wie kommen die Standardfehler von `sigma` und `rho` zustande?**
   (`compute_standard_errors()`) Die Inverse der Hesse-Matrix gibt die
   Kovarianz von `log(sigma)` und `atanh(rho)`. Die Delta-Methode
   transformiert mit der Jacobi-Matrix der Rücktransformation, die hier
   diagonal ist: `d sigma / d log(sigma) = sigma` und `d rho / d atanh(rho) =
   1 - rho^2`. `V = J V_u J'` reduziert sich damit auf ein elementweises
   Skalieren. Für `beta` und `gamma` ändert sich nichts, dort ist die
   Jacobi-Matrix die Einheitsmatrix. Test: `vcov(fit)` entspricht genau der
   so skalierten Inversen.

8. **Warum ein analytischer Gradient, wenn `optim` auch ohne auskommt?**
   (`compute_loglik_gradient()`) Bei der Prüfung der beiden Vorfassungen hielt
   BFGS mit dem numerischen Gradienten bei Zählwerten in den Tausenden
   vorzeitig an, bis zu 47 Log-Likelihood-Einheiten unter dem Optimum. Der
   Gradient der Log-Sum-Exp ist der mit `p_ik = exp(a_ik - l_i)` gewichtete
   Mittelwert der Knotengradienten; die Knotengradienten sind die einer
   Poisson-Log-Dichte, `(y - mu) x`, und von `log Phi(eta)`, nämlich `phi/Phi`
   mal die Ableitung von `eta`. Er ist gegen zentrale Differenzen getestet
   (Abweichung unter 1e-6) und macht den Fit drei- bis viermal schneller.

9. **Warum `optimHess` statt `numDeriv`?** (`compute_standard_errors()`)
   `optimHess` differenziert den analytischen Gradienten einmal numerisch,
   das ist genauer als die Log-Likelihood zweimal zu differenzieren.
   Gemessen: relative Abweichung der Standardfehler zur
   Richardson-Referenz unter 1e-4 schon ohne Gradient. Eine Abhängigkeit
   weniger für die Installation auf einem frischen R.

10. **Wie sind die Startwerte gewählt?** (`compute_start_values()`) Beide
    Gleichungen werden zunächst getrennt geschätzt, Probit für `gamma`,
    Poisson-GLM auf den selektierten Einheiten für `beta`, das ist das Modell
    unter `rho = 0`, also startet `rho` bei 0. `sigma` kommt aus der
    Überdispersion: unter dem Modell ist `Var(Y|x) = mu + mu^2 (exp(sigma^2)
    - 1)`. Das Poisson-GLM schätzt `log E[Y|x] = x'beta + sigma^2/2`, deshalb
    wird der Intercept um `sigma^2/2` nach unten korrigiert.

11. **Warum reicht `K = 20` nicht immer?** (`?poisselect`, Abschnitt
    „Accuracy of the quadrature") Der Integrand `p(y | mu(v)) Phi(.) phi(v)`
    wird mit wachsendem `y * sigma` immer schmaler; eine feste Knotenzahl
    kann eine schmale Spitze zwischen den Knoten verfehlen. Vergleich mit
    `integrate()`: bei `y = 25`, `sigma = 0.8` ist der Fehler mit `K = 20`
    0.17 pro Einheit, mit `K = 40` 0.005. Die approximierte Likelihood kann
    dann mehrere lokale Maxima haben; mit `K = 80` verschwinden sie in allen
    geprüften Fällen. Deshalb die Empfehlung, bei großen Zählwerten mit
    `K = 40` oder `80` nachzurechnen.

12. **Warum ist `predict(type = "response")` gleich `exp(x'beta +
    sigma^2/2)` und nicht `exp(x'beta)`?** (`predict.poisselect()`) Die
    Aufgabenstellung verlangt die unbedingte Populationserwartung `E[Y|x]`.
    Mit `mu = exp(x'beta + eps)` und `eps ~ N(0, sigma^2)` ist `E[exp(eps)]
    = exp(sigma^2/2)`, der Erwartungswert einer Lognormalverteilung. Das
    selektierte Mittel `E[Y|x, s = 1]` wäre wieder verzerrt.

13. **Wie entsteht Plot 1 und warum der Nenner?** (`compute_count_distribution()`)
    Für eine selektierte Einheit ist die bedingte Verteilung von `y` gegeben
    `s = 1` der Quotient aus der gemeinsamen Wahrscheinlichkeit von `y = m`
    und `s = 1` und der Selektionswahrscheinlichkeit; beides sind
    Quadratursummen mit denselben Knoten, im Nenner fehlt nur der
    Poisson-Faktor. Der Nenner hängt nicht von `m` ab und wird einmal
    berechnet; alles läuft auf der Log-Skala. Die Kurve ist das Mittel über
    die selektierten Einheiten und summiert sich bis `m = 40` auf über 0.98.

14. **Warum darf `y` bei `s = 0` fehlen, `x` nur bei `s = 0` und `z` nie?**
    (`check_model_data()`) Nicht-selektierte Einheiten gehen nur über
    `ln(1 - Phi(z'gamma))` in die Likelihood ein; ihre Outcome-Werte und
    Outcome-Kovariablen kommen dort nicht vor. Die Selektionskovariablen
    braucht jede Einheit. Deshalb wird `x` nur auf den selektierten Zeilen
    auf Vollständigkeit und Rang geprüft.

15. **Was ist die Exclusion Restriction und warum nur eine Warnung?**
    (`check_exclusion_restriction()`) Ohne eine Variable, die nur die
    Selektion beeinflusst, ist das Modell allein über die Nichtlinearität von
    `Phi` identifiziert; das ist formal gültig, praktisch instabil (Plot 2
    wird flach, `rho` bekommt riesige Standardfehler). Der Fit ist dann nicht
    falsch, nur schwach, deshalb Warnung statt Fehler. Der Intercept zählt
    nicht, er verschiebt die Selektionswahrscheinlichkeit für alle gleich.

16. **Wie wird perfekte Trennung erkannt?** (`fit_probit_start()`) Wenn
    eine Kombination der Selektionskovariablen die beiden Gruppen perfekt
    trennt, hat das Probit keinen endlichen ML-Schätzer; die IRLS-Iteration
    treibt die Koeffizienten ins Unendliche und die Devianz gegen 0. Die
    Startwert-Probit-Schätzung ist also ein billiger Detektor: Devianz
    numerisch 0 heißt Trennung, und dann warnt das Paket, statt `gamma = 408`
    mit Standardfehler 10000 auszugeben.

17. **Warum heißt das Argument `K` und nicht `k`, wo doch snake_case gilt?**
    (`poisselect.R`) Die Aufgabenstellung nennt die Knotenzahl `K`; wer die
    Aufgabe liest und `K = 10` tippt, soll nicht „unused argument" bekommen.
    Es ist die einzige Abweichung, sie steht an genau einer Zeile mit
    `# nolint: object_name_linter.` und einem Kommentar; intern heißt die
    Größe `n_nodes`.
