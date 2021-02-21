# Sigmoid-Function-
For more information please read the paper below:
https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b01171

Sigmoid function as the Reaction Coordinate of the Umbrella Sampling method on protein Trp-Cage

the reaction coordinate is based on the fraction of native
contacts. The set of native contacts is defined from the native
structure. Specifically, a pair of heavy atoms (i, j) in residues Ri
and Rj is counted as a native contact if |Ri − Rj| > 3 and the
interatomic distance rij_0 in the native structure is smaller than 4.5
Å. In our case, the number of native contacts identified from
the crystal structure is N = 156 and N = 279 for Trp-Cage. Assuming that the atom pair (i, j) is one of
the native contacts, we use rij(X) to denote the distance between the two atoms in a given protein conformation X. The
reaction coordinate Q for any conformation X is then determined by the distances for the N pairs of atoms in this
conformation
q(X) = sum( 1 / (1 + exp(-β (rij(x) - λ rij_0)))
Q(x) = q(x) / N 
with λ = 1.8 and a smoothing parameter β = 5.0 1/Å. The
summand in the equation above is effectively a pairwise contact
strength that approaches 1 when the distance rij is small and
approaches 0 when rij is large, thus quantifying the degree of
contact between the two atoms. The reaction coordinate (Q) is
the average over all pairwise contact strengths, thus
representing the collective fraction of the native contacts
present in a given conformation. A value of Q close to 1
indicates that the protein is in the native state because all of the
native contacts are intact. In contrast, Q ∼ 0 corresponds to
completely non-native structures with all the native contacts
broken.
