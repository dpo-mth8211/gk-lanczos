# MINRES

R_minr = []       # Espace pour les résidus aux moindres-carrés
X_minr = []       # Espace pour les valeurs de x à chaque itération
Niter_minr = []   # Espace pour le nombre d'itération
S_minr = []       # Espace pour le statut du problème

for p in problems 

  X = []  

  callback = workspace -> begin
    push!(X, copy(workspace.x))
    return false
  end

  # Extraction du problème sans enregistrement permanent
  A, b = get_mm(p)

  # Création du problème de point de selle avec perturbation
  n, m = size(A)
  δ = 1.0e-8
  K = [sparse(1.0I,n,n) A; A' -δ*I]  
  rhs = [b; zeros(m)];

  (x, stats) = minres(K, rhs, callback = callback ; history=true)

  # Ajout des données du problèmes en cours aux vecteurs alloués
  push!(R_minr, stats.residuals)
  push!(X_minr, X)
  push!(Niter_minr, stats.niter)
  push!(S_minr, stats.solved)
end

# MINRES-QLP

R_minrq = []       # Espace pour les résidus aux moindres-carrés
X_minrq = []       # Espace pour les valeurs de x à chaque itération
Niter_minrq = []   # Espace pour le nombre d'itération
S_minrq = []       # Espace pour le statut du problème

for p in problems 

  X = []  

  callback = workspace -> begin
    push!(X, copy(workspace.x))
    return false
  end

  # Extraction du problème sans enregistrement permanent
  A, b = get_mm(p)

  # Création du problème de point de selle avec perturbation
  n, m = size(A)
  δ = 1.0e-8
  K = [sparse(1.0I,n,n) A; A' -δ*I]  
  rhs = [b; zeros(m)];

  (x, stats) = minres_qlp(K, rhs, callback = callback ; history=true)

  # Ajout des données du problèmes en cours aux vecteurs alloués
  push!(R_minrq, stats.residuals)
  push!(X_minrq, X)
  push!(Niter_minrq, stats.niter)
  push!(S_minrq, stats.solved)
end

# MINARES

R_mina = []       # Espace pour les résidus aux moindres-carrés
X_mina = []       # Espace pour les valeurs de x à chaque itération
Niter_mina = []   # Espace pour le nombre d'itération
S_mina = []       # Espace pour le statut du problème

for p in problems 

  X = []  

  callback = workspace -> begin
    push!(X, copy(workspace.x))
    return false
  end

  # Extraction du problème sans enregistrement permanent
  A, b = get_mm(p)

  # Création du problème de point de selle avec perturbation
  n, m = size(A)
  δ = 1.0e-8
  K = [sparse(1.0I,n,n) A; A' -δ*I]  
  rhs = [b; zeros(m)];

  (x, stats) = minares(K, rhs, callback = callback ; history=true)

  # Ajout des données du problèmes en cours aux vecteurs alloués
  push!(R_mina, stats.residuals)
  push!(X_mina, X)
  push!(Niter_mina, stats.niter)
  push!(S_mina, stats.solved)
end

# SYMMLQ

R_sym = []       # Espace pour les résidus aux moindres-carrés
X_sym = []       # Espace pour les valeurs de x à chaque itération
Niter_sym = []   # Espace pour le nombre d'itération
S_sym = []       # Espace pour le statut du problème

for p in problems 

  X = []  

  callback = workspace -> begin
    push!(X, copy(workspace.x))
    return false
  end

  # Extraction du problème sans enregistrement permanent
  A, b = get_mm(p)

  # Création du problème de point de selle avec perturbation
  n, m = size(A)
  δ = 1.0e-8
  K = [sparse(1.0I,n,n) A; A' -δ*I]  
  rhs = [b; zeros(m)];

  (x, stats) = symmlq(K, rhs, callback = callback ; history=true)

  # Ajout des données du problèmes en cours aux vecteurs alloués
  push!(R_sym, stats.residuals)
  push!(X_sym, X)
  push!(Niter_sym, stats.niter)
  push!(S_sym, stats.solved)
end

# CRÉATION DU TABLEAU

# Récupère les m premières composantes de X, correspondant à r du problèmes de moindre norme
function pr(X)
  Pr = []
  push!(Pr, [x[1:1033] for x in X[1]])
  push!(Pr, [x[1:1850] for x in X[2]])
  push!(Pr, [x[1:1033] for x in X[3]])
  push!(Pr, [x[1:1850] for x in X[4]])
  return Pr
end

Pr_minr = pr(X_minr)
Pr_minrq = pr(X_minrq)
Pr_mina = pr(X_mina)
Pr_sym = pr(X_sym)

r0_minr = first.(R_minr)
r0_minrq = first.(R_minrq)
r0_mina = first.(R_mina)
r0_sym = first.(R_sym)
rk_minr = last.(R_minr)
rk_minrq = last.(R_minrq)
rk_mina = last.(R_mina)
rk_sym = last.(R_sym)

pr0_minr = norm.(first.(Pr_minr))
pr0_minrq = norm.(first.(Pr_minrq))
pr0_mina = norm.(first.(Pr_mina))
pr0_sym = norm.(first.(Pr_sym))
prk_minr = norm.(last.(Pr_minr))
prk_minrq = norm.(last.(Pr_minrq))
prk_mina = norm.(last.(Pr_mina))
prk_sym = norm.(last.(Pr_sym))

# Création du dataframe colonne par colonne
df = DataFrame()

df[:,:Meth] = ["MINRES", "MINRES", "MINRES", "MINRES", "MIN-QLP", "MIN-QLP", "MIN-QLP", "MIN-QLP", "MINARES", "MINARES", "MINARES", "MINARES", "SYMMLQ", "SYMMLQ", "SYMMLQ", "SYMMLQ"]
df[:, :Prob] = vcat(problems, problems, problems, problems)
df[:, :Nit] = vcat(Niter_minr, Niter_minrq, Niter_mina, Niter_sym)
df[:, :"||r0-s||"] = vcat(round.(r0_minr, sigdigits=4), round.(r0_minrq, sigdigits=4), round.(r0_mina, sigdigits=4), round.(r0_sym, sigdigits=4))
df[:, :"||rk-s||"] = vcat(round.(rk_minr,sigdigits=3), round.(rk_minrq,sigdigits=3), round.(rk_mina,sigdigits=3), round.(rk_sym,sigdigits=3))
df[:, :"||r0||"] = vcat(round.(pr0_minr, sigdigits=4), round.(pr0_minrq, sigdigits=4), round.(pr0_mina, sigdigits=4), round.(pr0_sym, sigdigits=4))
df[:, :"||rk||"] = vcat(round.(prk_minr,sigdigits=3), round.(prk_minrq,sigdigits=3), round.(prk_mina,sigdigits=3), round.(prk_sym,sigdigits=3))
df[:, :S] = vcat(S_minr, S_minrq, S_mina, S_sym)

# Affichage de la table
pretty_table(df, linebreaks=true)