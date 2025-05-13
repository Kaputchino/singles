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
    @variable(model, blacked[1:n, 1:n], Bin)   # 1 ↔ noir

        # 1) Unicité sur chaque ligne
    for i in 1:n
        for v in unique(grid[i, :])
            idx = findall(j -> grid[i,j] == v, 1:n)
            @constraint(model, sum(1 - blacked[i,j] for j in idx) ≤ 1)
        end
    end

    # 1′) Unicité sur chaque colonne
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

    # 3′) Toute case blanche a au moins un voisin blanc (contrainte tronqué de convexité)
    for i in 1:n, j in 1:n
        neigh = Tuple{Int,Int}[]
        i>1 && push!(neigh, (i-1,j))
        i<n && push!(neigh, (i+1,j))
        j>1 && push!(neigh, (i,j-1))
        j<n && push!(neigh, (i,j+1))
        @constraint(model,
            blacked[i,j] + sum(1 - blacked[p,q] for (p,q) in neigh) ≥ 1)
    end

    @objective(model, Min, 0)           # on cherche juste la faisabilité
    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(model)

    # Return:
    # 1 - true if an optimum is found
    # 2 - the resolution time
    return MOI.get(model, MOI.TerminationStatus()) == MOI.FEASIBLE_POINT || MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL, time() - start, blacked
    
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
"""
function solveDataSet()

    dataFolder = "../data/"
    resFolder = "../res/"
    #TODO know true size


    
    
    # Array which contains the name of the resolution methods
    resolutionMethod = ["cplex"]
    #resolutionMethod = ["cplex", "heuristique"]

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
                    writeSolution(outputFile,mat,solved_blacked)
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
                    
                    # While the grid is not solved and less than 100 seconds are elapsed
                    while !isOptimal && resolutionTime < 100
                        
                        # TODO 
                        println("In file resolution.jl, in method solveDataSet(), TODO: fix heuristicSolve() arguments and returned values")
                        
                        # Solve it and get the results
                        isOptimal, resolutionTime = heuristicSolve()

                        # Stop the chronometer
                        resolutionTime = time() - startingTime
                        
                    end

                    # Write the solution (if any)
                    if isOptimal

                        # TODO
                        println("In file resolution.jl, in method solveDataSet(), TODO: write the heuristic solution in fout")
                        
                    end 
                end

                println(fout, "solveTime = ", resolutionTime) 
                println(fout, "isOptimal = ", isOptimal)
                
                # TODO
                println("In file resolution.jl, in method solveDataSet(), TODO: write the solution in fout") 
                close(fout)
            end


            # Display the results obtained with the method on the current instance
            include(outputFile)
            println(resolutionMethod[methodId], " optimal: ", isOptimal)
            println(resolutionMethod[methodId], " time: " * string(round(solveTime, sigdigits=2)) * "s\n")
        end         
    end 
end
