# Résolution avec QR pour obtenir une approximation de x* pour calculer l'erreur
XQR = []
solve_QR(A, b) = qr(A) \ b

for p in problems
  A, b = get_mm(p)
  x_qr = solve_QR(A, b);
  push!(XQR, x_qr)
end

# Calcul de l'erreur à chaque itération
E_lsqr = []
for i in 1:length(problems)
  e = [norm(X_lsqr[i][j] .- XQR[i]) for j in 1:length(X_lsqr[i])]
  push!(E_lsqr, e)
end

E_lsmr = []
for i in 1:length(problems)
  e = [norm(X_lsmr[i][j] .- XQR[i]) for j in 1:length(X_lsmr[i])]
  push!(E_lsmr, e)
end

# Création des graphiques de chaque problème
p1 = plot([R_lsqr[1], R_lsmr[1], AR_lsqr[1], AR_lsmr[1], E_lsqr[1], E_lsmr[1]], xlims=(0,600), ylims=(0,8E3), label=["||rₖ|| - Lsqr" "||rₖ|| - Lsmr" "||A*rₖ|| - Lsqr" "||A*rₖ|| - Lsmr" "Eₖ - Lsqr" "Eₖ - Lsmr"], xlabel="Nb itérations",  title="illc1033");
p2 = plot([R_lsqr[2], R_lsmr[2], AR_lsqr[2], AR_lsmr[2], E_lsqr[2], E_lsmr[2]], xlims=(0,1100), ylims=(0,5E3), label=["||rₖ|| - Lsqr" "||rₖ|| - Lsmr" "||A*rₖ|| - Lsqr" "||A*rₖ|| - Lsmr" "Eₖ - Lsqr" "Eₖ - Lsmr"], xlabel="Nb itérations", title="illc1850");
p3 = plot([R_lsqr[3], R_lsmr[3], AR_lsqr[3], AR_lsmr[3], E_lsqr[3], E_lsmr[3]], label=["||rₖ|| - Lsqr" "||rₖ|| - Lsmr" "||A*rₖ|| - Lsqr" "||A*rₖ|| - Lsmr" "Eₖ - Lsqr" "Eₖ - Lsmr"], xlabel="Nb itérations", title="well1033");
p4 = plot([R_lsqr[4], R_lsmr[4], AR_lsqr[4], AR_lsmr[4], E_lsqr[4], E_lsmr[4]], xlims=(0,250), ylims=(0,1E4), label=["||rₖ|| - Lsqr" "||rₖ|| - Lsmr" "||A*rₖ|| - Lsqr" "||A*rₖ|| - Lsmr" "Eₖ - Lsqr" "Eₖ - Lsmr"], xlabel="Nb itérations", title="well1850");

