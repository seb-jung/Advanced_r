# Erklärung zur Nutzung von KI-Werkzeugen

**Kurs:** Advanced Data Analysis with R
**Abgabe:** R-Paket `poisselect` (Poisson-Selektionsmodell für Zähldaten)
**Verfasser:** Sebastian Jung, Finn Metzler <!-- beide: Namen prüfen -->
**Datum:** 03.09.2026

Diese Erklärung erfolgt gemäß Abschnitt 4 der KI-Richtlinie des Kurses.

---

## 1. Verwendete KI-Werkzeuge

Über Claude Code Desktop:

- **Claude Fable 5 (Anthropic)**, kostenpflichtige Version, privat finanziert;
  Erstellung der ersten Paketfassung von Finn Metzler (August 2026).
- **Claude Opus 5 (Anthropic)**, kostenpflichtige Version, privat finanziert;
  Verfeinerung der ersten Paketfassung von Finn durch Sebastian Jung (August 2026).
- **Claude Fable 5.1 (Anthropic)**, kostenpflichtige Version, privat finanziert;
  Prüfung, Vergleich und Zusammenführung der beiden Fassungen zur gemeinsamen Abgabe (September 2026).

---

## 2. Arbeitsschritte mit KI-Unterstützung

### 2.1 Erste Fassung von Finn Metzler (Claude Fable 5)

- **Recherche und Einordnung:** Auswertung der Vorlesungsunterlagen und der
  Aufgabenstellung im Hinblick auf die verbindlichen Kurskonventionen
  (Style Guide, Top-Down-Workflow, `checkmate`, S3, `lintr`).
- **Strukturierung des Pakets:** Entwurf der Top-Down-Gliederung, Aufteilung
  in Hauptfunktion und benannte Teilaufgaben, Zuordnung zu Dateien,
  Erstellung des Programmentwurfs (Pseudocode).
- **Generierung von R-Code** für alle Teilfunktionen, insbesondere
  Gauß-Hermite-Quadratur nach Golub-Welsch, approximierte Log-Likelihood auf
  der Log-Skala mit Log-Sum-Exp, Reparametrisierung `log(sigma)` /
  `atanh(rho)` mit Delta-Methode, Input-Checks mit `checkmate`, S3-Objekt und
  Methoden `print`, `summary`, `plot`, `predict`, `coef`, `vcov`, `logLik`.
- **Aufsetzen der testthat-Tests**, darunter eine naive Referenzimplementierung
  der Log-Likelihood, Parameterrückgewinnungstests und Tests der
  Fehlerbehandlung.
- **Debugging und Interpretation** von Ausgaben aus `R CMD check`, `lintr`
  und `testthat`.
- **Dokumentation:** Erstellung und sprachliche Überarbeitung der
  roxygen2-Dokumentation und der README.
- **Brainstorming zu Edge Cases:** `K = 1`, fehlende Exclusion Restriction,
  singuläre Hesse-Matrix, Over-/Underflow bei großen Zählwerten, NA-Muster
  zwischen Outcome- und Selektionsgleichung.

### 2.2 Prüfung, Vergleich und Zusammenführung (Claude Fable 5.1)

Dieser Schritt wurde von Sebastian Jung mit einem ausführlichen schriftlichen
Auftrag angestoßen und von Fable 5.1 weitgehend eigenständig ausgeführt.
Konkret hat das Werkzeug

- **beide Fassungen nach demselben Protokoll geprüft:** `R CMD build`,
  Installation per `remotes::install_local()` in eine frische Bibliothek,
  `R CMD check --as-cran`, vollständiger `testthat`-Lauf, `lintr` mit
  Default-Regeln, `roxygen2`-Abgleich, gerenderte Hilfeseiten;
- **unabhängige Prüfskripte geschrieben und ausgeführt:** eigener
  Datengenerator, naive Referenzimplementierung der Log-Likelihood direkt aus
  der Formel der Aufgabenstellung, Vergleich mit `integrate()`, Verifikation
  der Gauß-Hermite-Knoten gegen publizierte Werte und `statmod`,
  Monte-Carlo-Studie (100 Replikationen, drei `rho`-Werte, beide Pakete auf
  identischen Daten), eine Batterie von 73 fehlerhaften Inputs,
  Reproduzierbarkeitstests, Kreuzvergleich der Schätzergebnisse und eine
  Diagnose der Optimierer-Abweichungen;
- **die Zusammenführung entworfen und umgesetzt:** je Teilproblem die bessere
  der beiden Lösungen ausgewählt und begründet (`MERGE_ENTSCHEIDUNGEN.md`),
  den zusammengeführten Code, die Tests (Vereinigung beider Suiten plus neue
  Fälle) und die roxygen2-Dokumentation geschrieben;
- **neuen Code geschrieben, den keine der beiden Fassungen hatte:** den
  analytischen Gradienten der Log-Likelihood (gegen numerische Differenzen
  geprüft), die Erkennung perfekter Trennung über die Probit-Devianz, den
  Fehler bei ausschließlich Nullen im Outcome, die `NA`-Prüfung von `newdata`
  auf Variablenebene;
- **die Begleitdokumente entworfen:** `README.md`, `Pseudocode_poisselect.qmd`,
  `BEFUND.md`, `MERGE_ENTSCHEIDUNGEN.md` und den Entwurf dieser Erklärung.

Alle Prüfergebnisse in `BEFUND.md` beruhen auf tatsächlich ausgeführten
Läufen; die Skripte dazu liegen nicht in der Abgabe, ihre Ergebnisse sind im
Befund zitiert.

---

## 3. Eigenleistung, Prüfung und Verantwortung

<!-- Beide: Dieser Abschnitt ist eine Erklärung, die ihr abgebt. Lasst nur
     stehen, was zutrifft, und ergänzt, was ihr selbst geprüft habt. -->

Die Aufgabenstellung, die Modellgleichungen und die Auswahl der Verfahren
(Gauß-Hermite-Quadratur, Log-Sum-Exp, Reparametrisierung, Delta-Methode)
sind vorgegeben bzw. von uns festgelegt; die Umsetzung in R-Code wurde in
allen drei Schritten überwiegend von den genannten Werkzeugen erzeugt.

Wir haben den zusammengeführten Code gelesen und nachvollzogen und tragen die
inhaltliche Verantwortung für die Abgabe. Insbesondere haben wir geprüft:

- die implementierte Log-Likelihood gegen die Formel der Aufgabenstellung,
  einschließlich der Faktoren `sqrt(2)` und `1/sqrt(pi)` aus der Substitution
  `v = sqrt(2) t`;
- die Herleitung des analytischen Gradienten und seine Übereinstimmung mit
  numerischen Differenzen (Test `test-loglik.R`);
- die Gauß-Hermite-Knoten und -Gewichte gegen publizierte Werte für
  `K = 2, 3, 5`;
- die Parameterrückgewinnung auf simulierten Daten mit bekannten Parametern;
- die Fehlerbehandlung anhand der in `BEFUND.md` dokumentierten Fälle.

Die Entscheidungen in `MERGE_ENTSCHEIDUNGEN.md` haben wir nachvollzogen und
übernommen. <!-- Beide: falls ihr Entscheidungen geändert habt, hier nennen. -->
