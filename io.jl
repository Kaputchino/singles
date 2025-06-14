# This file contains functions related to reading, writing and displaying a grid and experimental results

using JuMP
using Plots
import GR

"""
Read an instance from an input file

- Argument:
inputFile: path of the input file
"""
function readInputFile(inputFile::String)::Matrix{Int}
    rows = [parse.(Int, split(strip(l))) for l in readlines(inputFile)]
    return hcat(rows...)'              # ⇒ Matrix{Int64}
end

function parse_result_file(path::String)
    lines = readlines(path)
    time = Inf
    optimal = false
    for line in lines
        if occursin("Temps", line)
            time = parse(Float64, split(line, ":")[2])
        elseif occursin("Statut", line)
            optimal = occursin("true", line)
        end
    end
    return time, optimal
end
function read_solution_grid(path::String)::Matrix{Union{Int, Char}}
    lines = readlines(path)
    grid = Union{Int, Char}[]  # liste plate à transformer après

    rows = 0
    cols = 0

    for line in lines
        if startswith(line, "Temps") || startswith(line, "Statut") || isempty(strip(line))
            break  # stop at metadata
        end
        row = [s == "X" ? 'X' : parse(Int, s) for s in split(line)]
        cols = length(row)
        append!(grid, row)
        rows += 1
    end

    return permutedims(reshape(grid, cols, rows))
end







"""
Create a pdf file which contains a performance diagram associated to the results of the ../res folder
Display one curve for each subfolder of the ../res folder.

Arguments
- outputFile: path of the output file

Prerequisites:
- Each subfolder must contain text files
- Each text file correspond to the resolution of one instance
- Each text file contains a variable "solveTime" and a variable "isOptimal"
"""
function performanceDiagram(outputFile::String)

    resultFolder = "../res/"
    
    # Maximal number of files in a subfolder
    maxSize = 0

    # Number of subfolders
    subfolderCount = 0

    folderName = Array{String, 1}()

    # For each file in the result folder
    for file in readdir(resultFolder)

        path = resultFolder * file
        
        # If it is a subfolder
        if isdir(path)
            
            folderName = vcat(folderName, file)
             
            subfolderCount += 1
            folderSize = size(readdir(path), 1)

            if maxSize < folderSize
                maxSize = folderSize
            end
        end
    end

    # Array that will contain the resolution times (one line for each subfolder)
    results = Array{Float64}(undef, subfolderCount, maxSize)

    for i in 1:subfolderCount
        for j in 1:maxSize
            results[i, j] = Inf
        end
    end

    folderCount = 0
    maxSolveTime = 0

    # For each subfolder
    for file in readdir(resultFolder)
            
        path = resultFolder * file
        
        if isdir(path)

            folderCount += 1
            fileCount = 0

        time, optimal = parse_result_file(path * "/" * resultFile)
        if optimal
            results[folderCount, fileCount] = time
            if time > maxSolveTime
                maxSolveTime = time
            end
        end
        end
    end 

    # Sort each row increasingly
    results = sort(results, dims=2)

    println("Max solve time: ", maxSolveTime)

    # For each line to plot
    for dim in 1: size(results, 1)

        x = Array{Float64, 1}()
        y = Array{Float64, 1}()

        # x coordinate of the previous inflexion point
        previousX = 0
        previousY = 0

        append!(x, previousX)
        append!(y, previousY)
            
        # Current position in the line
        currentId = 1

        # While the end of the line is not reached 
        while currentId != size(results, 2) && results[dim, currentId] != Inf

            # Number of elements which have the value previousX
            identicalValues = 1

             # While the value is the same
            while results[dim, currentId] == previousX && currentId <= size(results, 2)
                currentId += 1
                identicalValues += 1
            end

            # Add the proper points
            append!(x, previousX)
            append!(y, currentId - 1)

            if results[dim, currentId] != Inf
                append!(x, results[dim, currentId])
                append!(y, currentId - 1)
            end
            
            previousX = results[dim, currentId]
            previousY = currentId - 1
            
        end

        append!(x, maxSolveTime)
        append!(y, currentId - 1)

        # If it is the first subfolder
        if dim == 1

            # Draw a new plot
            plot(x, y, label = folderName[dim], legend = :bottomright, xaxis = "Time (s)", yaxis = "Solved instances",linewidth=3)

        # Otherwise 
        else
            # Add the new curve to the created plot
            savefig(plot!(x, y, label = folderName[dim], linewidth=3), outputFile)
        end 
    end
end 

"""
Create a latex file which contains an array with the results of the ../res folder.
Each subfolder of the ../res folder contains the results of a resolution method.

Arguments
- outputFile: path of the output file

Prerequisites:
- Each subfolder must contain text files
- Each text file correspond to the resolution of one instance
- Each text file contains a variable "solveTime" and a variable "isOptimal"
"""
function resultsArray(outputFile::String)
    
    resultFolder = "../res/"
    dataFolder = "../data/"
    
    # Maximal number of files in a subfolder
    maxSize = 0

    # Number of subfolders
    subfolderCount = 0

    # Open the latex output file
    fout = open(outputFile, "w", encoding="UTF-8")
    println("jecris dans le fichier",outputFile)
    # Print the latex file output
    println(fout, raw"""\documentclass{article}

\usepackage[french]{babel}
\usepackage [utf8] {inputenc} % utf-8 / latin1 
\usepackage{multicol}

\setlength{\hoffset}{-18pt}
\setlength{\oddsidemargin}{0pt} % Marge gauche sur pages impaires
\setlength{\evensidemargin}{9pt} % Marge gauche sur pages paires
\setlength{\marginparwidth}{54pt} % Largeur de note dans la marge
\setlength{\textwidth}{481pt} % Largeur de la zone de texte (17cm)
\setlength{\voffset}{-18pt} % Bon pour DOS
\setlength{\marginparsep}{7pt} % Séparation de la marge
\setlength{\topmargin}{0pt} % Pas de marge en haut
\setlength{\headheight}{13pt} % Haut de page
\setlength{\headsep}{10pt} % Entre le haut de page et le texte
\setlength{\footskip}{27pt} % Bas de page + séparation
\setlength{\textheight}{668pt} % Hauteur de la zone de texte (25cm)

\begin{document}""")

    header = raw"""
\begin{center}
\renewcommand{\arraystretch}{1.4} 
 \begin{tabular}{l"""

    # Name of the subfolder of the result folder (i.e, the resolution methods used)
    folderName = Array{String, 1}()

    # List of all the instances solved by at least one resolution method
    solvedInstances = Array{String, 1}()

    # For each file in the result folder
    for file in readdir(resultFolder)

        path = resultFolder * file
        
        # If it is a subfolder
        if isdir(path)

            # Add its name to the folder list
            folderName = vcat(folderName, file)
             
            subfolderCount += 1
            folderSize = size(readdir(path), 1)

            # Add all its files in the solvedInstances array
            for file2 in filter(x->occursin(".txt", x), readdir(path))
                solvedInstances = vcat(solvedInstances, file2)
            end 

            if maxSize < folderSize
                maxSize = folderSize
            end
        end
    end

    # Only keep one string for each instance solved
    unique(solvedInstances)

    # For each resolution method, add two columns in the array
    for folder in folderName
        header *= "rr"
    end

    header *= "}\n\t\\hline\n"

    # Create the header line which contains the methods name
    for folder in folderName
        header *= " & \\multicolumn{2}{c}{\\textbf{" * folder * "}}"
    end

    header *= "\\\\\n\\textbf{Instance} "

    # Create the second header line with the content of the result columns
    for folder in folderName
        header *= " & \\textbf{Temps (s)} & \\textbf{Optimal ?} "
    end

    header *= "\\\\\\hline\n"

    footer = raw"""\hline\end{tabular}
\end{center}

"""
    println(fout, header)

    # On each page an array will contain at most maxInstancePerPage lines with results
    maxInstancePerPage = 30
    id = 1

    # For each solved files
    for solvedInstance in solvedInstances
        println(solvedInstances)
        # If we do not start a new array on a new page
        if rem(id, maxInstancePerPage) == 0
            println(fout, footer, "\\newpage")
            println(fout, header)
        end 

        # Replace the potential underscores '_' in file names
        print(fout, replace(solvedInstance, "_" => "\\_"))

        # For each resolution method
        for method in folderName

            path = resultFolder * method * "/" * solvedInstance

            # If the instance has been solved by this method
            if isfile(path)

                solveTime, isOptimal = parse_result_file(path)
                grid = read_solution_grid(path)

                println(fout, " & ", round(solveTime, digits=2), " & ")

                if isOptimal
                    println(fout, "\$\\times\$")
                    println(fout, "\$\\times\$")
                end 
                
            # If the instance has not been solved by this method
            else
                println(fout, " & - & - ")
            end
        end

        println(fout, "\\\\")

        id += 1
    end

    # Print the end of the latex file
    println(fout, footer)

    println(fout, "\\end{document}")

    close(fout)
    
end 


function displayGrid(matrix::Matrix{Int})
    rows, cols = size(matrix)
    cell_width = maximum(length(string(x)) for x in matrix) + 1  # ajuster largeur dynamique

    # Ligne de séparation horizontale
    function horizontal_line()
        println("+" * join(["-"^cell_width for _ in 1:cols], "+") * "+")
    end

    # Affichage de chaque ligne
    horizontal_line()
    for i in 1:rows
        row_str = "|"
        for j in 1:cols
            val = string(matrix[i, j])
            padding = " " ^ (cell_width - length(val))
            row_str *= padding * val * "|"
        end
        println(row_str)
        horizontal_line()
    end
end

function displaySolution(grid::Matrix{Int}, solution::Matrix{Bool})
    rows, cols = size(grid)
    cell_width = maximum(length(string(x)) for x in grid) + 1  # ajuster largeur dynamique

    # Ligne de séparation horizontale
    function horizontal_line()
        println("+" * join(["-"^cell_width for _ in 1:cols], "+") * "+")
    end

    # Affichage de chaque ligne
    horizontal_line()
    for i in 1:rows
        row_str = "I"
        for j in 1:cols
            if(solution[i,j])
                val = "x"
            else
                val = string(grid[i, j])
            end
            padding = " " ^ (cell_width - length(val))
            row_str *= padding * val * "I"
        end
        println(row_str)
        horizontal_line()
    end
end

function writeSolution(path::AbstractString,grid::Matrix{Int},blacked::Matrix{Bool},isOptimal::Bool,time)
    print_black(blacked)
    @assert size(grid) == size(blacked) "dimensions différentes"
    open(path, "w") do io
        for i in 1:size(grid,1)
            row = [blacked[i,j] ? "X" : string(grid[i,j])
                   for j in 1:size(grid,2)]
            println(io, join(row, " "))
        end
        println(io)  
        println(io, "Temps de résolution : $(time)")
        println(io, "Statut : ", isOptimal)
        #println(io, "solveTime = ", resolutionTime) 
        #println(io, "isOptimal = ", isOptimal)
    end

end

function print_black(black::Matrix{Bool})
    h, w = size(black)
    for i in 1:h
        for j in 1:w
            print(black[i, j] ? "X " : ". ")
        end
        println()
    end
end
