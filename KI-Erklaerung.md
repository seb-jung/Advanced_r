# Erklärung zur Nutzung von KI-Werkzeugen

Dieses Dokument dokumentiert die Nutzung von Künstlicher Intelligenz (KI) bei der Erstellung des R-Pakets `poisselect`, gemäß den Richtlinien des Kurses _Advanced Data Analysis with R_.

## Verwendete Werkzeuge

Für die Entwicklung wurden verschiedene Versionen von **Claude (Anthropic)** genutzt, insbesondere zur Code-Generierung, Fehlerbehebung und zur Prüfung von Implementierungsdetails.

## Arbeitsschritte und KI-Unterstützung

Die Entwicklung erfolgte in drei wesentlichen Phasen:

1.  **Planung und Strukturierung:** Die Projektstruktur ergab sich direkt aus dem top-down-orientierten Programmentwurf. Diese Struktur wurde von uns unter Verwendung verschiedener KI-Modelle umgesetzt, wobei wir kontinuierlich den Abgleich mit den Projektzielen sichergestellt haben.
2.  **Implementierung:** Gemäß den KI-Richtlinien wurden einzelne Funktionen – insbesondere Boilerplate-Code – durch KI unterstützt erstellt. Die Steuerung des Workflows und die inhaltliche Korrektheit der Logik blieben dabei stets in unserer Verantwortung. Der Großteil des R-Codes (Kernfunktionen, S3-Methoden, Gauß-Hermite-Quadratur, Log-Likelihood-Berechnungen) wurde unter Verwendung von KI-generierten Vorschlägen erstellt, die anschließend manuell geprüft wurden.
3.  **Verifikation und Testing:** Die KI wurde eingesetzt, um unabhängige Prüfskripte zu schreiben (z.B. Monte-Carlo-Studien, Vergleiche mit Referenzformeln) und die `testthat`-Testsuite zu erweitern.

## Eigenleistung und Verantwortung

Obwohl die KI maßgeblich bei der Code-Erzeugung unterstützt hat, tragen die Verfasser die volle Verantwortung für die inhaltliche Korrektheit der Abgabe. Folgende Schritte wurden manuell geprüft und verifiziert:

- **Mathematische Korrektheit:** Abgleich der implementierten Log-Likelihood und des analytischen Gradienten mit den in der Aufgabenstellung vorgegebenen Formeln.
- **Validierung der Ergebnisse:** Durchführung von Parameterrückgewinnungstests auf simulierten Daten sowie Verifikation der Gauß-Hermite-Knoten gegen publizierte Werte.
- **Fehlerbehandlung:** Manuelle Prüfung, ob alle in der Aufgabenstellung geforderten Fehlerfälle (z.B. fehlende Exclusion Restriction, perfekte Trennung, ungültige Inputs) korrekt erkannt und mit informativen Meldungen abgefangen werden.
- **Qualitätssicherung:** Überprüfung des Codes auf Einhaltung des Style Guides, DRY-Prinzip und Lesbarkeit.

---
