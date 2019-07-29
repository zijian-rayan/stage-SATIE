# stage-SATIE
Traitement de signaux de vibrométrie Laser  par méthodes statistiques/probabiliste

elliptic_plot.m：Ce fichier est une fonction qui correspond à l'ellipse

fmax2.m：Une fonction qui renvoie toutes les valeurs maximales dans un ensemble de données (matrices). La valeur maximale est calculée par interpolation quadratique.

direction_3.m: Ce fichier est utilisé pour sélectionner le seuil optimal pour le filtrage du domaine fréquentiel dans trois directions et pour quantifier l’effet de filtrage.

direction_z.m： Ce fichier est utilisé pour essayer d'extraire automatiquement la limite, la direction z étant utilisée comme référence.

prony_t.m： Ce fichier contient du code permettant de traiter les signaux à l’aide de la méthode Prony, notamment la recherche de pôles et de vitesses de groupe. J'utilise les données dans la direction x.

sans_ref.m： Ce fichier est le fichier d'origine qui correspond à l'ellipse

swam_identification_v2.m： Ce fichier contient le code de test de la méthode prony.

test_prony.m： Ce fichier contient un code de test supplémentaire pour la méthode Prony, ce qui prouve que la méthode Prony peut reconstruire le signal.

trois_diections_bruit.m: Ce fichier contient les effets de filtrage des signaux dans trois directions dans les domaines temps et fréquence, et s'insérant dans une ellipse
