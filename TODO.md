# Meeting 12-12-2023:
## TODO :
- [x] Ajouter Column Generation
- [x] Faire des restarts de la subgradient method
- [x] Pourquoi Bundle Proximal Level Method ne marche pas ? en étude.
  - x0 est le best point trop longtemps donc on ne cesse de le projeter sur des level sets différents et il n'y a pas de progression. Je relance BPLM en ne changeant que la façon dont les level sets sont construits
- [x] Pourquoi Column-and-Row Generation ne marche pas ? 
  - Problème dans la génération du modèle qui prend un temps fou, il faut aller voir dans les JuMP-related issues.
- [x] Explications textuelles de ce que je calcule et de ce que je plot sur Overleaf.
- [x] Pour corriger FGM, tester GD sur le cas smooth, tracker la distance entre f tilde et f après un run donné.
  - 