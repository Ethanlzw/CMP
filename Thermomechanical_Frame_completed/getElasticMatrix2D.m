function D = getElasticMatrix2D(E, nu, mechType)

switch lower(mechType)
%% Complete the following session
    case 'planestress'
        D = E / (1 - nu^2) * [1  nu 0;
                              nu 1  0;
                              0  0  (1 - nu)/2];
    case 'planestrain'
        D = E / ((1 + nu) * (1 - 2*nu)) * ...
            [1 - nu, nu,      0;
             nu,     1 - nu,  0;
             0,      0,       (1 - 2*nu)/2];
    otherwise
        error('Unknown mechanical type. Use planeStress or planeStrain.');
end

end