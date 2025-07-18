# CRAIG

R_craig = []      # Espace pour les résidus aux moindres-carrés
X_craig = []      # Espace pour les valeurs de x à chaque itération
Niter_craig = []  # Espace pour le nombre d'itération
S_craig = []      # Espace pour le statut du problème

for p in problems 

  X = []  

  callback = workspace -> begin
    push!(X, copy(workspace.x))
    return false
  end

  # Extraction du problème sans enregistrement permanent
  A = get_mm(p)[1]

  A = A'
  u = rand(size(A,2))
  b = A * u           # b est dans l'image de A

  (x, y, stats) = craig(A, b, callback = callback ; history=true)

  # Ajout des données du problèmes en cours aux vecteurs alloués
  push!(R_craig, stats.residuals)
  push!(X_craig, X)
  push!(Niter_craig, stats.niter)
  push!(S_craig, stats.solved)
end

# LSQR

R_lsqr = []      # Espace pour les résidus aux moindres-carrés
X_lsqr = []      # Espace pour les valeurs de x à chaque itération
Niter_lsqr = []  # Espace pour le nombre d'itération
S_lsqr = []      # Espace pour le statut du problème

for p in problems 

  X = []  

  callback = workspace -> begin
    push!(X, copy(workspace.x))
    return false
  end

  # Extraction du problème sans enregistrement permanent
  A = get_mm(p)[1]

  A = A'
  u = rand(size(A,2))
  b = A * u           # b est dans l'image de A

  (x, stats) = lsqr(A, b, callback = callback ; history=true)

  # Ajout des données du problèmes en cours aux vecteurs alloués
  push!(R_lsqr, stats.residuals)
  push!(X_lsqr, X)
  push!(Niter_lsqr, stats.niter)
  push!(S_lsqr, stats.solved)
end

# LSMR

R_lsmr = []      # Espace pour les résidus aux moindres-carrés
X_lsmr = []      # Espace pour les valeurs de x à chaque itération
Niter_lsmr = []  # Espace pour le nombre d'itération
S_lsmr = []      # Espace pour le statut du problème

for p in problems 

  X = []  

  callback = workspace -> begin
    push!(X, copy(workspace.x))
    return false
  end

  # Extraction du problème sans enregistrement permanent
  A = get_mm(p)[1]

  A = A'
  u = rand(size(A,2))
  b = A * u           # b est dans l'image de A

  (x, stats) = lsmr(A, b, callback = callback ; history=true)

  # Ajout des données du problèmes en cours aux vecteurs alloués
  push!(R_lsmr, stats.residuals)
  push!(X_lsmr, X)
  push!(Niter_lsmr, stats.niter)
  push!(S_lsmr, stats.solved)
end

# CRÉATION DU TABLEAU

r0_lsqr = first.(R_lsqr)
r0_lsmr = first.(R_lsmr)
r0_craig = first.(R_craig)
rk_lsqr = last.(R_lsqr)
rk_lsmr = last.(R_lsmr)
rk_craig = last.(R_lsqr)

x0_lsqr = norm.(first.(X_lsqr))
x0_lsmr = norm.(first.(X_lsmr))
x0_craig = norm.(first.(X_craig))
xk_lsqr = norm.(last.(X_lsqr))
xk_lsmr = norm.(last.(X_lsmr))
xk_craig = norm.(last.(X_lsqr))

# Création du dataframe colonne par colonne
df = DataFrame()

df[:,:Meth] = ["LSQR", "LSQR", "LSQR", "LSQR", "LSMR", "LSMR", "LSMR", "LSMR", "CRAIG", "CRAIG", "CRAIG", "CRAIG"]
df[:, :Prob] = vcat(problems, problems, problems)
df[:, :Nit] = vcat(Niter_lsqr, Niter_lsmr, Niter_craig)
df[:, :"||x0||"] = vcat(round.(x0_lsqr, sigdigits=3), round.(x0_lsmr, sigdigits=3), round.(x0_craig, sigdigits=3))
df[:, :"||xk||"] = vcat(round.(xk_lsqr, sigdigits=3), round.(xk_lsmr, sigdigits=3), round.(xk_craig, sigdigits=3))
df[:, :"||r0||"] = vcat(round.(r0_lsqr, sigdigits=3), round.(r0_lsmr, sigdigits=3), round.(r0_craig, sigdigits=3))
df[:, :"||rk||"] = vcat(round.(rk_lsqr, sigdigits=3), round.(rk_lsmr, sigdigits=3), round.(rk_craig, sigdigits=3))
df[:, :S] = vcat(S_lsqr, S_lsmr, S_craig)

# Affichage de la table
pretty_table(df, linebreaks=true)
