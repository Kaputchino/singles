include("generation.jl")
include("resolution.jl")
include("io.jl")

function launch()
    println("Lancement complet de la chaîne : génération → résolution → export Latex")

    # Étape 1 : génération d'une instance
    println("Génération d'une instance...")
    generateDataSet()  # ou une autre fonction de génération spécifique si tu en as une
    println("Instances générées.\n")

    # Étape 2 : résolution
    println("Résolution des instances...")
    solveDataSet()
    println("Résolution terminée.\n")

    println("Génération du tableau Latex...")
    resultsArray("./res/resultats.tex")
    println("Fichier Latex généré dans ./res/resultats.tex\n")

    println("Processus complet terminé.")
end
