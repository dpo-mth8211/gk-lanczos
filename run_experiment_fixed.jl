function run_experiment(matrix_name::String, experiment_type::String, algos::Vector{String})
    A, bmat = get_mm(matrix_name)
    b = bmat[:, 1]
    n, m = size(A)

    if experiment_type == "least_squares"
        println("Problème aux moindres carrés avec $algos pour $matrix_name")
        x_sol = A \ b
        err_map = Dict{String, Vector{Float64}}()
        stats_map = Dict{String, Any}()
        
        # on applique pour chaque methode 
        for algo in algos
            err_map[algo] = Float64[]
            solver = getproperty(Krylov, Symbol(algo))
            callback = ws -> (push!(err_map[algo], norm(x_sol - ws.x)); false)
            x, stats = solver(A, b; history=true, callback=callback)
            stats_map[algo] = stats
        end
        
        # on crée le tableau de résultats
        headers = ["", algos...]
        col_names = ["# itérations", "‖r‖ (initial)", "‖r‖ (final)", "‖Aᵗr‖ (initial)", "‖Aᵗr‖ (final)"]
        data = col_names
        
        # on extrait des statistiques pour chaque algorithme
        for algo in algos
            stats = stats_map[algo]
            stats_vec = [stats.niter, stats.residuals[1], stats.residuals[end],
                       stats.Aresiduals[1], stats.Aresiduals[end]]
            data = hcat(data, stats_vec)
        end

        # on affiche le tableau
        pretty_table(data; header=headers, tf=tf_markdown)

        # on affiche ici le statut
        for algo in algos
            println("Statut final pour $algo : ", string(stats_map[algo].status))
        end

        # on réalise le tracé des graphiques
        gr()
        p = plot(; yaxis=:log, xlabel="Itération", legend=:topright)
        for algo in algos
            stats = stats_map[algo]
            plot!(p, stats.residuals, label="‖r‖ $algo", lw=2)
            plot!(p, stats.Aresiduals, label="‖Aᵗr‖ $algo", lw=2, linestyle=:dash)
            plot!(p, err_map[algo], label="‖err‖ $algo", lw=2, linestyle=:dot)
        end
        savefig(p, "LS_$matrix_name.pdf")

    elseif experiment_type in ["saddle_primal", "saddle_adjoint"]
        is_primal = experiment_type == "saddle_primal"
        system_type = is_primal ? "primal" : "adjoint"
        println("Système de point de selle ($system_type) avec $algos pour $matrix_name")

        # On construit le système de point de selle
        δ = 0.0
        if is_primal
            K = [sparse(I, n, n) A; A' -δ * I]
            rhs = [b; zeros(m)]
        else
            K = [sparse(I, m, m) A'; A -δ * I]
            rhs = [b[1:m]; zeros(n)]
        end

        Aresiduals = Dict{String, Vector{Float64}}()
        stats_map = Dict{String, Any}()

        for algo in algos
            Aresiduals[algo] = Float64[]
            solver = getproperty(Krylov, Symbol(algo))
            callback = ws -> begin
                x_block = is_primal ? ws.x[1:n] : ws.x[1:m]
                res = is_primal ? A' * x_block : A * x_block
                push!(Aresiduals[algo], norm(res))
                return false
            end
            x, stats = solver(K, rhs; history=true, callback=callback)
            stats_map[algo] = stats
        end

        headers = ["", algos...]
        col_names = ["# itérations", "‖r‖ (initial)", "‖r‖ (final)", "‖Aᵗr‖ (initial)", "‖Aᵗr‖ (final)"]
        data = col_names

        for algo in algos
            stats = stats_map[algo]
            if algo == "symmlq"
                stats_vec = [stats.niter, stats.residuals[1], stats.residuals[end], "-", "-"]
            else
                stats_vec = [stats.niter, stats.residuals[1], stats.residuals[end],
                          stats.Aresiduals[1], stats.Aresiduals[end]]
            end
            data = hcat(data, stats_vec)
        end

        pretty_table(data; header=headers, tf=tf_markdown)
        for algo in algos
            println("Statut final pour $algo : ", string(stats_map[algo].status))
        end

        gr()
        p = plot(; yaxis=:log, xlabel="Itération", legend=:topright)
        for algo in algos
            stats = stats_map[algo]
            xaxis = 1:length(Aresiduals[algo])
            plot!(p, stats.residuals, label="‖r‖ $algo", lw=2)
            plot!(p, xaxis, Aresiduals[algo], label="‖Ar LS‖ $algo", lw=2, linestyle=:dash)
        end
        savefig(p, "AUG_$(system_type)_$matrix_name.pdf")

    elseif experiment_type == "minimum_norm"
        println("Problème de moindre norme avec $algos pour $matrix_name")

        # On prépare le problème de moindre norme
        A = A'
        b = b[1:size(A, 1)]
        x_sol = lq(Matrix(A)) \ b

        errors = Dict{String, Vector{Float64}}()
        norms = Dict{String, Vector{Float64}}()
        stats_map = Dict{String, Any}()

        # On applique les méthodes
        for algo in algos
            errors[algo] = Float64[]
            norms[algo] = Float64[]
            solver = getproperty(Krylov, Symbol(algo))
            callback = ws -> begin
                push!(errors[algo], norm(x_sol - ws.x))
                push!(norms[algo], norm(ws.x))
                return false
            end
            
            # Traitement spécial pour craig qui retourne un tuple (x, y, stats)
            if algo == "craig"
                x, y, stats = solver(A, b; history=true, callback=callback)
            else
                x, stats = solver(A, b; history=true, callback=callback)
            end
            stats_map[algo] = stats
        end

        headers = ["", algos...]
        col_names = ["# itérations", "‖err‖ (initial)", "‖err‖ (final)", "‖x‖ (initial)", "‖x‖ (final)"]
        data = col_names

        for algo in algos
            stats = stats_map[algo]
            stats_vec = [stats.niter,
                       errors[algo][1], errors[algo][end],
                       norms[algo][1], norms[algo][end]]
            data = hcat(data, stats_vec)
        end

        pretty_table(data; header=headers, tf=tf_markdown)
        for algo in algos
            println("Statut final pour $algo : ", string(stats_map[algo].status))
        end

        gr()
        p = plot(; yaxis=:log, xlabel="Itération", legend=:topright)
        for algo in algos
            plot!(p, errors[algo], label="‖err‖ $algo", lw=2)
            plot!(p, norms[algo], label="‖x‖ $algo", lw=2, linestyle=:dash)
        end
        
        savefig(p, "MoindreNorme_$matrix_name.pdf")
    else
        error("Type d'expérience traitée ici inconnu : $experiment_type")
    end
end
