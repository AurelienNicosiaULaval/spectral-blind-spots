# Plan de comparaison fixé avant simulation

9 octobre 2026. Étude de diagnostics, sans étude humaine et sans classement universel d'imputeurs. Les contrôles de géométrie ont été exécutés avant ce plan. Aucune fréquence de rejet du benchmark n'a encore été observée lors de sa rédaction.

## Données et hypothèse

Paires indépendantes (X,Y,Z). La référence Y et la candidate Z sont tirées indépendamment conditionnellement à X. Deux tailles : n=240 et n=720. Le prédicteur est uniforme sauf les deux alternatives locales et la nullité conditionnelle, où sa densité vaut (1+0.75 cos x)/(2 pi). L'hypothèse de fidélité complète impose la même loi conditionnelle à Y et Z. Cette hypothèse sera testée avec deux nullités : réponse uniforme ; réponse dépendante avec perturbation 0.8 sin(x-theta).

Alternatives, candidate uniforme : perturbation sin[k(x-theta)] avec k=1, amplitude 0.4 ; k=2, amplitudes 0.4 et 0.8 ; k=4 et k=6, amplitude 0.8 ; k=7, amplitude 0.8, hors de la bande de conception des angles. Ajouter une perturbation 0.8 sin(2x)sin(2theta), puis deux bosses locales d'amplitude 0.8, de demi-largeur 0.9, centrées respectivement en pi et en 0, multipliant sin(2theta). Toutes les densités sont non négatives.

500 jeux indépendants par condition nulle et taille ; 200 par alternative et taille. Total : 5600 jeux. Graines déterministes : 2026101000 + 10000*indice_condition + indice_replication, avec conditions ordonnées par scénario puis taille. Le coût a été testé uniquement sur une multiplication matricielle synthétique ; aucun scénario du benchmark n'a été utilisé pour sélectionner les angles, bandes ou amplitudes.

## Comparateurs

Les neuf diagnostics appariés reçoivent exactement les mêmes X,Y,Z et 199 vecteurs indépendants de signes. Ils utilisent une statistique maximale standardisée sauf le MMD quadratique. Le niveau demandé est 5 %, avec p=(1+nombre de dépassements)/(B+1). Les différences dans le nombre de composantes sont consignées ; un budget d'information identique n'est pas un coût de calcul identique.

1. Quatre vues cosinus croisées aux angles 0 et pi/2, équivalentes aux quatre vues sinus/cosinus canoniques.
2. Ces mêmes vues après quatre rotations communes : 0, pi/8, pi/3 et pi/4, soit seize vues.
3. Quatre vues croisées aux angles 0 et pi/4.
4. Quatre vues croisées aux angles 0 et pi/7, choisies pour la bande de conception 1 à 6. La fréquence 7 est un échec prévu.
5. Quatre vues croisées aux angles 0 et pi/sqrt(2), avec absence de résonance exacte mais possibilité de quasi-résonance.
6. Contrastes de Fourier complets : fréquences prédicteur -8 à 8 et réponse 1 à 8, parties réelle et imaginaire.
7. Fourier marginal de réponse, ordres 1 à 8.
8. Champ périodique : grille 16 par 32, bandes prédicteur 0.35 et 0.70, bandes réponse 0.15 et 0.30, support effectif minimal 20 ; maximum commun sur les quatre bandes et le support retenu. Noyaux wrapped-normal normalisés sur la grille de réponse.
9. MMD apparié quadratique avec noyau produit exp[-(1-cos(dx))/0.5²] exp[-(1-cos(dtheta))/0.5²]. Le noyau est fixé, sans optimisation sur les résultats.

Les diagnostics de vues utilisent les moments des polynômes de Chebyshev des coordonnées projetées jusqu'à l'ordre 8, soit cos[k(x-a)]cos[l(theta-b)], k=0,...,8 et l=1,...,8. Ils interrogent des moments de la loi projetée, pas sa densité complète. La fréquence 7 est donc dans les moments analysés, tout en étant hors de la bande utilisée pour concevoir l'angle pi/7. La tronquature commune est explicite.

Un dixième diagnostic, supplémentaire, utilise la CDF exacte de la candidate : couverture locale à 19 seuils de probabilité, mêmes poids périodiques du prédicteur, standardisation par variance Bernoulli connue et 199 rééchantillonnages uniformes de PIT sous le modèle. Il suit le principe de couverture conditionnelle, mais ne sera pas présenté comme la reproduction de l'algorithme LCT de Zhao et al. Il dispose de la loi candidate connue et génère une référence sous cette loi, alors que les neuf autres comparent deux tirages. Ses résultats seront séparés dans l'interprétation.

## Mesures et incertitude

Rapporter taux de rejet et intervalles de Wilson à 95 %, avec le nombre de jeux. Rapporter aussi les différences appariées des décisions, leur erreur Monte Carlo, les temps de calcul et le nombre de composantes. Sous les alternatives projetées invisibles, ne pas appeler le taux observé un niveau exact : la validité par échanges est garantie sous la fidélité complète, et l'invariance simultanée de toutes les vues n'est pas déduite de leurs lois marginales.

Pour le champ, consigner le nombre de positions actives, l'identité avec le témoin normalisé et la grille. Les rejets ne valident pas une loi imputée ; une absence de rejet ne prouve pas sa fidélité. Les mêmes graines servent à tous les comparateurs d'un jeu, sans choisir le diagnostic après observation.

## Décision et conservation

Conserver tous les résultats, y compris les échecs du choix pi/7, les faibles puissances du champ et les écarts de niveau. La comparaison doit établir ce que changent les représentations, bandes et hypothèses dans ce plan. Elle ne démontre ni compréhension humaine ni originalité à elle seule.

Le calcul sera interrompu pour une erreur d'implémentation ou de validité, pas pour un résultat défavorable. Toute modification de plan nécessitée par une erreur sera datée et identifiée. Les résultats préliminaires corrigés seront distingués du benchmark final.
