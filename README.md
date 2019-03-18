
## Bugs

* Edge struct is a mess, I guess I should give it some member functions
* Need to think about traverse of neighborhood: I was using boundary_v == BOTH as
  the check of on boundary but it not true in fact, problem?
* Then do i still need to use the constraint plane if I have INVALID?
* This program currently never simplifiy non-boundary edge whose endpoints are on boundary
* What should be the weight of constraint plane quadric? result looks bad!!

## To-dos

* with color and texture: Shiwei says it does not help much
* cleanup wrong topology: instead of skipping ecol which changes topo, collapse possibly
  more than one edges and change topo to collapse the target edge finally
* optional fix boundary: implement the paper method
* openmp: for loop in pre processing
* aspect ratio: while scanning neighborhood, examine aspect ratio and lower its error if
  too large, so that this long face will become the target in next ecol
* ecol coding

## Observations

1. Face area weighting produces remarkably less bad-shaped triangles
