# Description générale du projet

## Contexte large

Ce projet s'intéresse à l'application de techniques modernes d’inférence statistique (Approximate Bayesian Computing) pour une question centrale en biologie : la dynamique de génération de diversité génétique (mutagénèse) dans des populations bactériennes (par exemple la génération d'un variant conférant une résistance accrue à un antibiotique).

Les approches habituellement utilisées pour caractériser cette génération de diversité génétique s'appuient sur l'estimation du taux de mutations à partir d'expériences simples réalisées *in vitro* en laboratoire. Cette estimation se fait par maximum likelihood (vraissemblance maximale) grace à un modèle analytique simple.

Ce modèle est cependant trop simpliste (c'est à dire qu'il fait des hypothèse trop restrictives) pour répondre aux questions ouvertes en évolution microbienne. Il lui manque notamment la possibilité de prendre en compte :

- l'écologie complexe caratérisant l'environnement dans lequel vivent habituellement les bactéries (e.g. au sein d'un microbiome complexe et dans un hôte doté d'un système immunitaire plutôt que dans un tube à essais;)

- un éventuel traitement antibiotique

- l'effet potentiel des mutations sur la dynamique de croissance


Pour dépasser les limitations de ce modèle historique, nous proposons l'utilisation de simulations incorporant ces différentes extensions du modèle standard. L'inférence des différents paramètres à partir de données expérimentales peut alors être faite par des méthodes de calcul bayésien approché (Approximate Bayesian Computing -- ABC)). Il s'agit en bref de méthodes permettant de déterminer quels paramètres d'une simulation permettent de se rapprocher au maximum des données expérimentales, avec des contraintes de parcimonie (ne pas appeler le simulateur plus que nécessaire), de convergence et d'intervalles de confiance.

Nous avons pour cela développé un simulateur super-méga-rapide (3 ordres de grandeurs plus rapides que l'état de l'art), capable de s'affranchir des limites du modèle historique détaillées ci-dessus, et de simuler une dynamique écologique arbitrairement complexe. Ce projet s'appuiera sur ce simulateur, avec donc comme objectif de développer une méthode d'inférence rapide et efficace qui poutrera l'état de l'art.

Pour plus de détails, voir un résumé scientifique plus poussé dans `atreyu_short.pdf`


# Description des fichiers fournis

- `gillespie.cpp` (et `forward_simulator.h`) : un simulateur (obtenu suivant la méthode standard -- algorithme de Gillespie) pour le problème biologique qui nous intéresse (émergence de mutants pendant la croissance d'une population); la syntaxe est indiquée en appelant avec l'option `-h`
- `atreyu_forward_simulator.cpp`, `forward_simulator.h` : un simulateur nettement plus rapide au prix de quelques approximations biologiquement réalistes
- `atreyu.py` : des fonctions d'analyse en vrac, pour l'instant je vous fourni le strict minimum, l'idée étant que vous trouverez des meilleurs solutions que moi
- `run_*_tests.py` (scripts python) et `*.test` (fichiers de paramètres): un système de tests basiques des simulateurs (cf contenu des scripts et du Makefile)
- `requirements.txt` : pas fait de façon très fine (ne pas hésiter à corriger), mais indique au moins les versions que j'ai utilisées
- `Makefile` non seulement pour la compilation mais aussi pour lancer les tests, ne pas hésiter à y rajouter des choses

# Description du modèle biologique
![description du modèle biologique de croissance avec mutation](model.png)

# Propositions générales de taches à réaliser

Tout ceci peut être adapté à vos envies, et vous pouvez vouloir proposer d'autres choses. Tenez moi au courant des taches que vous choisissez d'aborder (par exemple en ouvrant une issue) et de vos propositions techniques.

- Développement d’une méthode d'ABC (calcul bayésien approché) pour estimer le taux de mutation lorsque les autres paramètres (paramètres écologiques) sont connus à partir de données expérimentales, en s'appuyant sur le simulateur déjà développé.

- Extension de cette méthode à l'inférence simultanée du taux de mutation et de paramètres écologiques inconnus. Le nombre de populations expérimentales nécessaires à l'inférence simultanée de deux paramètres ou plus et l'identifiabilité (ou non) de chacun des paramètres seront d'abord étudiés par des simulations.

- (exploratoire) Extension du simulateur à l'étude de dynamiques plus complexes (eg plusieurs catégories de mutants)

- (très exploratoire) Amélioration de l’estimateur ABC en remplaçant le simulateur direct par un modèle d’apprentissage profond (modèle de substitution – surrogate model – qui sont utilisés de manière croissante en physique mais pas encore popularisés en biologie) entraîné sur un grand jeu de données de simulation.

# Premier data challenge pour prise en main du projet

Dans `datachallengeA`, 384 simulations ont été réalisées avec les paramètres suivants : 
- taille de population initiale 10
- taille de population finale 1e7
- pas de mort (valeur du paramètre 0)
- effet sur la fitness (croissance relative du mutant) 0.7
- pas d'echantillonage (valeur du paramètre 1)
- taux de mutation inconnu : **pouvez-vous le retrouver ?**


