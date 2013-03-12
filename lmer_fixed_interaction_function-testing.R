# Testing ###########
library(lme4)
data(cake)
cake1 <- lmer(
  angle ~ recipe * temperature + (1|recipe:replicate),
  data=cake
)

GetVarInfo(cake1, "recipeB:temperature.L")
recB_tempL_FI <- FixedInteraction(cake1, "recipeB:temperature.L")


# Try with custom-specified contrasts to see whether
# having unnamed contrasts breaks things
cake$recipe2 = cake$recipe
contrasts(cake$recipe2) = cbind(c(1, 0, 0), c(0, 0, 1))

cake3 = lmer(
  angle ~ recipe2 * temperature + (1|recipe:replicate),
  data=cake
)

GetVarInfo(cake3, "recipe1:temperature.L")

# Try with temperature as numeric instead
cake2 <- lmer(
  angle ~ recipe * temp + (1|recipe:replicate),
  data=cake)
recB_temp_FI <- FixedInteraction(cake2, "recipeB:temp")