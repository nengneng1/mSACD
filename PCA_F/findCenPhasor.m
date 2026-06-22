

function [G_mean, S_mean] = findCenPhasor(sph_ref)
G_mean = mean(sph_ref.G(abs(sph_ref.G)>=1e-05));
S_mean = mean(sph_ref.S(abs(sph_ref.S)>=1e-05));
end