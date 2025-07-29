# LSQR

R_lsqr = []       # Espace pour les résidus aux moindres-carrés
X_lsqr = []       # Espace pour les valeurs de x à chaque itération
AR_lsqr = []      # Espace pour les résidus d'optimalité
Niter_lsqr = []   # Espace pour le nombre d'itération
S_lsqr = []       # Espace pour le statut du problème

for p in problems 

  X = []  

  callback = workspace -> begin
    push!(X, copy(workspace.x))
    return false
  end

  # Extraction du problème sans enregistrement permanent
  A, b = get_mm(p)
  (x, stats) = lsqr(A, b, callback = callback ; history=true)

  # Ajout des données du problèmes en cours aux vecteurs alloués
  push!(R_lsqr, stats.residuals)
  push!(AR_lsqr, stats.Aresiduals)
  push!(X_lsqr, X)
  push!(Niter_lsqr, stats.niter)
  push!(S_lsqr, stats.solved)
end

# LSMR

R_lsmr = []       # Espace pour les résidus aux moindres-carrés
X_lsmr = []       # Espace pour les valeurs de x à chaque itération
AR_lsmr = []      # Espace pour les résidus d'optimalité
Niter_lsmr = []   # Espace pour le nombre d'itération
S_lsmr = []       # Espace pour le statut du problème

for p in problems 

  X = []  

  callback = workspace -> begin
    push!(X, copy(workspace.x))
    return false
  end

  # Extraction du problème sans enregistrement permanent
  A, b = get_mm(p)
  (x, stats) = lsmr(A, b, callback = callback ; history=true)

  # Ajout des données du problèmes en cours aux vecteurs alloués
  push!(R_lsmr, stats.residuals)
  push!(AR_lsmr, stats.Aresiduals)
  push!(X_lsmr, X)
  push!(Niter_lsmr, stats.niter)
  push!(S_lsmr, stats.solved)
end

# CRÉATION DU TABLEAU

# Récupère les données demandées (première et dernière itération)
r0_lsqr = first.(R_lsqr)
r0_lsmr = first.(R_lsmr)
rk_lsqr = last.(R_lsqr)
rk_lsmr = last.(R_lsmr)
a0_lsqr = first.(AR_lsqr)
a0_lsmr = first.(AR_lsmr)
ak_lsqr = last.(AR_lsqr)
ak_lsmr = last.(AR_lsmr)

# Création du dataframe colonne par colonne
df = DataFrame()

df[:,:Meth] = ["LSQR", "LSQR", "LSQR", "LSQR", "LSMR", "LSMR", "LSMR", "LSMR"]
df[:, :Prob] = vcat(problems, problems)
df[:, :Nit] = vcat(Niter_lsqr, Niter_lsmr)
df[:, :"||r0||"] = vcat(round.(r0_lsqr,sigdigits=4), round.(r0_lsmr,sigdigits=4))
df[:, :"||rk||"] = vcat(round.(rk_lsqr,sigdigits=3), round.(rk_lsmr,sigdigits=3))
df[:, :"||A*r0||"] = vcat(round.(a0_lsqr,sigdigits=5), round.(a0_lsmr,sigdigits=5))
df[:, :"||A*rk||"] = vcat(round.(ak_lsqr,sigdigits=2), round.(ak_lsmr,sigdigits=2))
df[:, :Status] = vcat(S_lsqr, S_lsmr)

# Affichage de la table
pretty_table(df, linebreaks=true)