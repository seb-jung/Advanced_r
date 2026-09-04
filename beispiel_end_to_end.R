# Ende-zu-Ende-Beispiel fuer das Paket poisselect.
#
# Durchlaeuft Schaetzung, print, summary, beide Plots, alle drei predict-Typen
# sowie coef, vcov, logLik, AIC und BIC. Ausfuehren mit
#   Rscript beispiel_end_to_end.R
# nach der Installation des Pakets. Die Plots landen in
# beispiel_end_to_end_plots.pdf im Arbeitsverzeichnis.

library(poisselect)

# 1. Schaetzung auf dem mitgelieferten Datensatz (wahre Werte: siehe
#    ?sim_selection). K = 20 ist der Default.
fit <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = sim_selection)

# 2. print und summary
print(fit)
summary(fit)

# 3. Beide Plots
pdf("beispiel_end_to_end_plots.pdf", width = 10, height = 5)
plot(fit)
plot(fit, which = 1, max_count = 10)
plot(fit, which = 2, n_grid = 81)
dev.off()

# 4. Alle drei predict-Typen, auf den Fit-Daten und auf neuen Daten
head(predict(fit))                        # E[Y | x] = exp(x'beta + sigma^2/2)
head(predict(fit, type = "link"))         # x'beta
head(predict(fit, type = "pselect"))      # Phi(z'gamma)
neu <- data.frame(x1 = c(-1, 0, 1), x2 = 0, z1 = c(-1, 0, 1))
cbind(neu,
      response = predict(fit, newdata = neu),
      link = predict(fit, newdata = neu, type = "link"),
      pselect = predict(fit, newdata = neu, type = "pselect"))

# 5. Extraktoren
coef(fit)
coef(fit, which = "outcome")
round(vcov(fit), 4)
logLik(fit)
AIC(fit)
BIC(fit)

# 6. Vergleich mit dem naiven Poisson-GLM auf den selektierten Einheiten:
#    dessen Intercept ist nach oben verzerrt (rho > 0).
naiv <- glm(y ~ x1 + x2, family = poisson, data = subset(sim_selection, s == 1))
rbind(wahr = c(0.5, 0.8, -0.4),
      poisselect = coef(fit, which = "outcome"),
      naives_glm = coef(naiv))

# 7. Eigene Daten simulieren und die Sensitivitaet gegenueber K pruefen
set.seed(1)
daten <- simulate_poisselect(n = 1500, rho = 0.6)
fit20 <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = daten)
fit40 <- poisselect(y ~ x1 + x2, s ~ x1 + z1, data = daten, K = 40)
rbind(K20 = coef(fit20), K40 = coef(fit40))
c(loglik_K20 = logLik(fit20), loglik_K40 = logLik(fit40))
