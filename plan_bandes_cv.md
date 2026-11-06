# Sensibilité : bandes choisies sur l'apprentissage

Plan fixé le 30 octobre avant le calcul. Pour les 600 paires d'apprentissage de l'application, choisir les deux bandes d'un estimateur conditionnel de densité par validation croisée leave-one-day-out : exclusion de toutes les observations du même jour, puis maximisation de la somme des log-densités conditionnelles. Grille fixée a priori : hx en {0.35, 0.70, 1.20}, ht en {0.15, 0.30, 0.60}, en radians. Il ne s'agit pas d'une sélection optimisée pour le test de contraste, ni d'une correction complète de la dépendance temporelle entre jours.

Utiliser les 276 paires d'évaluation et les 50 mêmes tirages archivés pour recalculer le champ sur les bandes sélectionnées et sur (0.70,0.30). Comparer leurs contrastes en L2 discret, maxima et supports, sans sélectionner a posteriori les résultats à rapporter. Un échec ou une densité prédictive nulle est conservé et documenté.
