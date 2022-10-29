# MEC8211 Projet Devoirs 1 et 2 Vérification de code SOTTAS-GOUREAU 

## Informations générales
Ce projet a été réalisé avec le langage de programmation Python. Le thème de ce projet porte sur le processus de diffusion du sel de mer dans un pilier cylindrique de béton.
Les principaux objectifs de ce projet sont:
* Étudier l'évolution de la concentration de sel à l'intérieur du pilier en béton en fonction du temps
* Établir un code de différences finies en fonction des données du problème
* Appliquer les principes de vérification de code

### Contexte du problème
La responsable de projet de la compagnie DUPONT & Associés Inc., Mme d'AVIGNON nous charge d'étudier l'évolution de la concentration de sel dans une structure poreuse telle qu'un pilier de béton d'un pont.

### Fichiers

Deux fichiers sont disponibles pour le devoir 1:
* CodeDevoirCartésien.py traite le problème en coordonnées cylindriques
* CodeDevoirCylindrique.py traite le problème en coordonnées cartésiennes

Quatre nouveaux fichiers sont disponibles pour le devoir 2:
* [ComparaisonCOMSOL.py](ComparaisonCOMSOL.py) effectue une comparaison code à code entre la solution proposée par le code Python et celle obtenue par le logiciel COMSOL
* [ComparaisonCOMSOL2.py](ComparaisonCOMSOL2.py) effectue également une comparaison code à code entre la solution proposée par le code Python et celle obtenue par le logiciel COMSOL en effectuant un changement de pas temporel à 1e6s
* [CodeDevoir2_Espace.py](CodeDevoir2_Espace.py) réalise une analyse de convergence spatiale pour un pas de temps fixe
* [CodeDevoir2_Temps.py](CodeDevoir2_Temps.py) réalise une analyse de convergence temporelle pour un pas de maillage fixe

Les fichiers de comparaisons code à code sont accompagnés de 6 fichiers textes de données issues des simulations réalisées sur COMSOL. 

### Résultats

Les résultats sont disponibles dans un fichier Powerpoint.


