function Flux = GodunovFlux(A,absA,UL,UR)

Flux = 0.5*A*(UL+UR) - 0.5*absA*(UR-UL);

end