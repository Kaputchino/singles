# This file contains methods to solve an instance (heuristically or with CPLEX)
using GLPK 
import MathOptInterface as MOI

include("generation.jl")

TOL = 0.00001

"""
Solve an instance with CPLEX
"""
function cplexSolve(grid::Matrix{Int})

    # Create the model
    model = Model(GLPK.Optimizer)
    n = size(grid, 1)


    # TODO
    println("In file resolution.jl, in method cplexSolve(), TODO: fix input and output, define the model")
    @variable(model, blacked[1:n, 1:n], Bin)   # 1 noir

        # 1) Unicité sur chaque ligne
    for i in 1:n
        for v in unique(grid[i, :])
            idx = findall(j -> grid[i,j] == v, 1:n)
            @constraint(model, sum(1 - blacked[i,j] for j in idx) ≤ 1)
        end
    end

    # 1) Unicité sur chaque colonne
    for j in 1:n
        for v in unique(grid[:, j])
            idx = findall(i -> grid[i,j] == v, 1:n)
            @constraint(model, sum(1 - blacked[i,j] for i in idx) ≤ 1)
        end
    end

    # 2) Aucune paire de noirs adjacents
    for i in 1:n, j in 1:n-1
        @constraint(model, blacked[i,j] + blacked[i,j+1] ≤ 1)
    end
    for i in 1:n-1, j in 1:n
        @constraint(model, blacked[i,j] + blacked[i+1,j] ≤ 1)
    end

    # 3) Toute case blanche a au moins un voisin blanc (contrainte tronqué de convexité)
    for i in 1:n, j in 1:n
        voisin = Tuple{Int,Int}[]
        i>1 && push!(voisin, (i-1,j))
        i<n && push!(voisin, (i+1,j))
        j>1 && push!(voisin, (i,j-1))
        j<n && push!(voisin, (i,j+1))
        @constraint(model,
            blacked[i,j] + sum(1 - blacked[p,q] for (p,q) in voisin) ≥ 1)
    end

    @objective(model, Min, 0)          
    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())

    if status == MOI.OPTIMAL || status == MOI.FEASIBLE_POINT
        B_bool = Matrix{Bool}([value(blacked[i,j]) > 0.5 for i in 1:n, j in 1:n])
        print_black_only(B_bool)
        return true, time() - start, B_bool
    else
        println("Pas de solution trouvée par CPLEX (statut = $status)")
        B_bool = fill(false, n, n) 
        return false, time() - start, B_bool
    end 
end

"""
Heuristically solve an instance
"""
function heuristicSolve()

    # TODO
    println("In file resolution.jl, in method heuristicSolve(), TODO: fix input and output, define the model")
    
end 

"""
Solve all the instances contained in "../data" through CPLEX and heuristics

The results are written in "../res/cplex" and "../res/heuristic"

Remark: If an instance has previously been solved (either by cplex or the heuristic) it will not be solved again

A changer usine à gaz
"""
function solveDataSet()

    dataFolder = "./data/"
    resFolder = "./res/"
    #TODO know true size


    
    
    # Array which contains the name of the resolution methods
    #resolutionMethod = ["cplex"]
    resolutionMethod = ["cplex", "heuristique"]

    # Array which contains the result folder of each resolution method
    resolutionFolder = resFolder .* resolutionMethod

    # Create each result folder if it does not exist
    for folder in resolutionFolder
        if !isdir(folder)
            mkdir(folder)
        end
    end
            
    global isOptimal = false
    global solveTime = -1

    # For each instance
    # (for each file in folder dataFolder which ends by ".txt")
    for file in filter(x->occursin(".txt", x), readdir(dataFolder))  
        
        println("-- Resolution of ", file)
        mat = readInputFile(dataFolder * file)

        # TODO
        println("In file resolution.jl, in method solveDataSet(), TODO: read value returned by readInputFile()")
        
        # For each resolution method
        for methodId in 1:size(resolutionMethod, 1)
            
            outputFile = resolutionFolder[methodId] * "/" * file

            # If the instance has not already been solved by this method
            if !isfile(outputFile)
                
                fout = open(outputFile, "w")  

                resolutionTime = -1
                isOptimal = false
                
                # If the method is cplex
                if resolutionMethod[methodId] == "cplex"
                    
                    # TODO 
                    println("In file resolution.jl, in method solveDataSet(), TODO: fix cplexSolve() arguments and returned values")
                    
                    # Solve it and get the results
                    isOptimal, resolutionTime,solved_blacked = cplexSolve(mat)
                    writeSolution(outputFile,mat,solved_blacked,isOptimal, resolutionTime)
                    # If a solution is found, write it
                    if isOptimal
                        # TODO
                        println("In file resolution.jl, in method solveDataSet(), TODO: write cplex solution in fout") 
                    end

                # If the method is one of the heuristics
                else
                    
                    isSolved = false

                    # Start a chronometer 
                    startingTime = time()
                    
                    # While the grid is not solved and less than 30 seconds are elapsed
                    while !isOptimal && resolutionTime < 30
                        
                        # TODO 
                        println("In file resolution.jl, in method solveDataSet(), TODO: fix heuristicSolve() arguments and returned values")
                        
                        # Solve it and get the results
                        solved_blacked, isOptimal = solveByGreedyScoredHeuristic(mat)

                        println("nothing crash that far")
                        resolutionTime = time() - startingTime

                        writeSolution(outputFile,mat,solved_blacked,isOptimal, resolutionTime)
                        println("write didn't crash the whole program YEAH")
                        # Stop the chronometer
                        
                    end

                    # Write the solution (if any)
                    if isOptimal

                        # TODO
                        println("In file resolution.jl, in method solveDataSet(), TODO: write the heuristic solution in fout")
                        
                    end 
                end

                #println(fout, "solveTime = ", resolutionTime) 
                #println(fout, "isOptimal = ", isOptimal)

                # TODO
                println("In file resolution.jl, in method solveDataSet(), TODO: write the solution in fout") 
                close(fout)
            end


            # Display the results obtained with the method on the current instance
            solveTime, isOptimal = parse_result_file(outputFile)
            grid = read_solution_grid(outputFile)
            println(resolutionMethod[methodId], " optimal: ", isOptimal)
            println(resolutionMethod[methodId], " time: " * string(round(solveTime, sigdigits=2)) * "s\n")
        end         
    end 
end

# Vérifie si une position est dans la grille
function inbounds(i, j, h, w)
    return 1 ≤ i ≤ h && 1 ≤ j ≤ w
end



function solveByHeuristic(grid::Matrix{Int64}; max_seconds::Float64 = 2.0)::Tuple{Matrix{Bool}, Bool}
    n, m = size(grid)
    is_blacked = fill(false, n, m)  # force Matrix{Bool}

    println("search for a solution...")

    start_time = time()
    max_seconds = max_seconds
    has_solution = backtrack!(grid, is_blacked, start_time, max_seconds)

    return is_blacked, has_solution
end


function backtrack!(
    grid::Matrix{Int64},
    is_blacked::Matrix{Bool},
    start_time::Float64,
    max_seconds::Float64
)::Bool
    if (time() - start_time) > max_seconds
        println("⏱Temps limite dépassé.")
        return false
    end
    println("backtrack call: time = $(round(time() - start_time, digits=2))s")

    n, m = size(grid)

    # Construire une liste des cases avec doublons
    score = Dict{Tuple{Int, Int}, Int}()
    for i in 1:n, j in 1:m
        if is_blacked[i, j]
            continue
        end

        val = grid[i, j]
        row_dupes = count(k -> k != j && grid[i, k] == val && !is_blacked[i, k], 1:m)
        col_dupes = count(k -> k != i && grid[k, j] == val && !is_blacked[k, j], 1:n)

        total_dupes = row_dupes + col_dupes
        if total_dupes > 0
            score[(i, j)] = total_dupes
        end
    end

    # Si plus aucun doublon, tester la validité finale
    if isempty(score)
        return is_valid_black(is_blacked) && is_white_connected(is_blacked)
    end

    # Trier les cases à essayer
    sorted_cells = sort(collect(keys(score)), by = pos -> -score[pos])

    for (i, j) in sorted_cells
        #  Essayer en noir
        is_blacked[i, j] = true
        if is_valid_black(is_blacked) && is_white_connected(is_blacked)
            if backtrack!(grid, is_blacked, start_time, max_seconds)
                return true
            end
        end
        is_blacked[i, j] = false

        #  Essayer en blanc
        if backtrack!(grid, is_blacked, start_time, max_seconds)
            return true
        end

        return false  # Aucun des deux choix n’a fonctionné
    end

    return false
end


function has_duplicates(grid::Matrix{Int64}, is_blacked::Matrix{Bool})::Bool
    n, m = size(grid)
    for i in 1:n
        seen = Dict{Int64, Bool}()
        for j in 1:m
            if !is_blacked[i, j]
                if haskey(seen, grid[i, j])
                    return true
                end
                seen[grid[i, j]] = true
            end
        end
    end
    for j in 1:m
        seen = Dict{Int64, Bool}()
        for i in 1:n
            if !is_blacked[i, j]
                if haskey(seen, grid[i, j])
                    return true
                end
                seen[grid[i, j]] = true
            end
        end
    end
    return false
end

function has_duplicate_at(grid::Matrix{Int64}, is_blacked::Matrix{Bool}, i::Int, j::Int)::Bool
    val = grid[i, j]
    n, m = size(grid)
    for k in 1:m
        if k != j && grid[i, k] == val && !is_blacked[i, k]
            return true
        end
    end
    for k in 1:n
        if k != i && grid[k, j] == val && !is_blacked[k, j]
            return true
        end
    end
    return false
end

function is_valid_black(is_blacked::Matrix{Bool})::Bool
    n, m = size(is_blacked)
    for i in 1:n, j in 1:m
        if is_blacked[i, j]
            for (di, dj) in ((1,0), (-1,0), (0,1), (0,-1))
                ni, nj = i + di, j + dj
                if 1 ≤ ni ≤ n && 1 ≤ nj ≤ m && is_blacked[ni, nj]
                    return false
                end
            end
        end
    end
    return true
end

function is_white_connected(is_blacked::Matrix{Bool})::Bool
    n, m = size(is_blacked)
    visited = falses(n, m)

    # Trouver un point blanc de départ
    start = nothing
    for i in 1:n, j in 1:m
        if !is_blacked[i, j]
            start = (i, j)
            break
        end
    end
    if start === nothing
        return false
    end

    # BFS
    queue = [start]
    visited[start...] = true
    while !isempty(queue)
        i, j = pop!(queue)
        for (di, dj) in ((1,0), (-1,0), (0,1), (0,-1))
            ni, nj = i + di, j + dj
            if 1 ≤ ni ≤ n && 1 ≤ nj ≤ m && !is_blacked[ni, nj] && !visited[ni, nj]
                visited[ni, nj] = true
                push!(queue, (ni, nj))
            end
        end
    end

    # Vérifie que toutes les cases blanches sont visitées
    for i in 1:n, j in 1:m
        if !is_blacked[i, j] && !visited[i, j]
            return false
        end
    end
    return true
end


function print_black_only(black::Matrix{Bool})
    h, w = size(black)
    for i in 1:h
        for j in 1:w
            print(black[i, j] ? "X " : ". ")
        end
        println()
    end
end

function solveByGreedyScoredHeuristic(grid::Matrix{Int})::Tuple{Matrix{Bool}, Bool}
    n, m = size(grid)
    is_blacked = fill(false, n, m)

    println("Heuristique gloutonne avec score en cours...")

    while has_duplicates(grid, is_blacked)
        candidates = []

        for i in 1:n, j in 1:m
            if is_blacked[i, j]
                continue
            end

            if has_duplicate_at(grid, is_blacked, i, j)
                score = compute_black_score(grid, is_blacked, i, j)
                push!(candidates, ((i, j), score))
            end
        end

        if isempty(candidates)
            println("Aucun candidat restant.")
            return finish_with_possible_rollback(grid, is_blacked)
        end

        # Trier par score croissant
        sort!(candidates, by = x -> x[2])
        placed = false

        for ((i, j), _) in candidates
            is_blacked[i, j] = true
            if is_valid_black(is_blacked) && is_white_connected(is_blacked)
                println("Case noircie intelligemment en ($i,$j)")
                placed = true
                break
            else
                is_blacked[i, j] = false
            end
        end

        if !placed
            println("Aucun placement possible sans violation.")
            return finish_with_possible_rollback(grid, is_blacked)
        end
    end

    println("Fin de l’heuristique, vérification...")

    return finish_with_possible_rollback(grid, is_blacked)
end

function compute_black_score(grid::Matrix{Int}, is_blacked::Matrix{Bool}, i::Int, j::Int)::Int
    val = grid[i, j]
    n, m = size(grid)

    row_conflicts = count(k -> k != j && grid[i, k] == val && !is_blacked[i, k], 1:m)
    col_conflicts = count(k -> k != i && grid[k, j] == val && !is_blacked[k, j], 1:n)

    white_neighbors = 0
    for (di, dj) in ((-1,0), (1,0), (0,-1), (0,1))
        ni, nj = i + di, j + dj
        if 1 ≤ ni ≤ n && 1 ≤ nj ≤ m && !is_blacked[ni, nj]
            white_neighbors += 1
        end
    end

    return (row_conflicts + col_conflicts) * 10 - white_neighbors
end

function attempt_rollback(grid::Matrix{Int}, is_blacked::Matrix{Bool})::Bool
    black_positions = [(i, j) for i in 1:size(grid,1), j in 1:size(grid,2) if is_blacked[i, j]]

    println("Tentative rollback simple sur $(length(black_positions)) cases...")

    for (i, j) in black_positions
        is_blacked[i, j] = false
        if is_valid_black(is_blacked) && is_white_connected(is_blacked) && !has_duplicates(grid, is_blacked)
            println("Rollback simple réussi en supprimant ($i,$j)")
            return true
        end
        is_blacked[i, j] = true
    end
    println("Rollback simple échoué.")
    return false
end

function attempt_double_rollback(grid::Matrix{Int}, is_blacked::Matrix{Bool})::Bool
    black_positions = [(i, j) for i in 1:size(grid,1), j in 1:size(grid,2) if is_blacked[i, j]]

    println("↩Tentative rollback double sur $(length(black_positions)) cases...")

    for i1 in 1:length(black_positions)
        for i2 in i1+1:length(black_positions)
            (x1, y1) = black_positions[i1]
            (x2, y2) = black_positions[i2]

            is_blacked[x1, y1] = false
            is_blacked[x2, y2] = false

            if is_valid_black(is_blacked) && is_white_connected(is_blacked) && !has_duplicates(grid, is_blacked)
                println("Rollback double réussi en retirant ($x1,$y1) et ($x2,$y2)")
                return true
            end

            is_blacked[x1, y1] = true
            is_blacked[x2, y2] = true
        end
    end

    println("Rollback double échoué.")
    return false
end

function finish_with_possible_rollback(grid::Matrix{Int}, is_blacked::Matrix{Bool})::Tuple{Matrix{Bool}, Bool}
    if is_valid_black(is_blacked) && is_white_connected(is_blacked) && !has_duplicates(grid, is_blacked)
        println("Solution valide.")
        return is_blacked, true
    end

    println("Tentative de rollback simple (final)...")
    if attempt_rollback(grid, is_blacked)
        return is_blacked, true
    end

    println("Tentative de rollback double (final)...")
    if attempt_double_rollback(grid, is_blacked)
        return is_blacked, true
    end

    println("Aucune solution valide trouvée après tous les rollbacks.")
    return is_blacked, false
end
