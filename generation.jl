# This file contains methods to generate a data set of instances (i.e., sudoku grids)
include("io.jl")
using Dates

using Random

# un carré latin est techniquement un cas particulier d'un jeux de singles ou il n'y a aucune case à noircir
function generateLatinSquare(n::Int)
    # Créer une permutation de base des nombres 1 à n
    base = collect(1:n)
    
    # Créer la première ligne comme une permutation aléatoire
    first_row = randperm(n)
    
    # Initialiser la matrice avec la première ligne
    square = [first_row]
    
    # Construire les lignes suivantes par décalage circulaire
    for i in 2:n
        new_row = circshift(first_row, i - 1)
        push!(square, new_row)
    end

    # Mélanger les lignes de façon aléatoire (optionnel)
    shuffle!(square)

    # Transposer puis mélanger les colonnes (optionnel)
    square = permutedims(hcat(square...), randperm(n))

    return square
end

function modifyMatrixNTimes(matrix::Matrix{Int}, y::Int)
    rows, cols = size(matrix)

    # Coins à exclure
    corners = Set([(1,1), (1,cols), (rows,1), (rows,cols)])
    modified_positions = Set{Tuple{Int, Int}}()

    for attempt in 1:y
        # Générer la liste des positions valides à chaque étape
        valid_positions = []

        for i in 1:rows, j in 1:cols
            pos = (i, j)

            # Exclure les coins et déjà modifiés
            if pos in corners || pos in modified_positions
                continue
            end

            # Vérifier que les voisins ne sont pas modifiés
            neighbors = [
                (i-1,j), (i+1,j), (i,j-1), (i,j+1)
            ]

            has_modified_neighbor = any(
                (ni, nj) in modified_positions &&
                1 ≤ ni ≤ rows && 1 ≤ nj ≤ cols
                for (ni, nj) in neighbors
            )

            if !has_modified_neighbor
                push!(valid_positions, pos)
            end
        end

        # Si aucune position valide, on arrête
        if isempty(valid_positions)
            println("Modification arrêtée après $(attempt - 1) modifications : plus de positions valides.")
            break
        end

        # Choisir une position valide au hasard
        pos = rand(valid_positions)
        i, j = pos
        current_value = matrix[i, j]

        # Choisir une nouvelle valeur différente
        row_values = setdiff(matrix[i, :], [current_value])
        new_value_candidates = length(row_values) > 0 ? row_values : setdiff(unique(matrix), [current_value])

        if isempty(new_value_candidates)
            continue  # Pas de valeur de remplacement possible, ignorer cette tentative
        end

        new_value = rand(new_value_candidates)
        matrix[i, j] = new_value
        push!(modified_positions, pos)
    end

    return matrix
end






"""
Generate an n*n grid with a given density

Argument
- n: size of the grid
- density: percentage in [0, 1] of initial values in the grid
"""
function generateInstance(n::Int64, density::Float64)
    y = round(Int, density * (n * n - 4) )
    return modifyMatrixNTimes(generateLatinSquare(n),y);
end 

"""
Generate all the instances

Remark: a grid is generated only if the corresponding output file does not already exist
"""

function generateDataSet(
    n::Int = 5,
    density::Float64 = 0.2,
    count::Int = 10,
    folder::String = "./data/"
)
    # Crée le dossier s’il n'existe pas
    isdir(folder) || mkpath(folder)

    for i in 1:count
        matrix = generateInstance(n, density)

        # Timestamp avec millisecondes
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS.sss")

        # Nom du fichier avec n, densité, index et timestamp
        filename = "instance_$(n)x$(n)_$(round(Int, density*100))pct_$i_$timestamp.txt"
        filepath = joinpath(folder, filename)

        # Écriture de la matrice dans le fichier
        open(filepath, "w") do io
            for row in matrix
                println(io, join(row, " "))
            end
        end
    end

    println("Instances generes")
end




