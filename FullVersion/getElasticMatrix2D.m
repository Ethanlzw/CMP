function D = getElasticMatrix2D(E, nu, mechType)
% Return constitutive matrix for 2D linear elasticity
%
% strain vector = [eps_xx; eps_yy; gamma_xy]

if nargin < 3 || isempty(mechType)
    mechType = 'planeStress';
end

switch lower(mechType)
    case lower('planeStress')
        D = E / (1 - nu^2) * ...
            [1,   nu,            0;
             nu,  1,             0;
             0,   0,   (1 - nu)/2];

    case lower('planeStrain')
        D = E / ((1 + nu) * (1 - 2*nu)) * ...
            [1 - nu,  nu,               0;
             nu,      1 - nu,           0;
             0,       0,       (1 - 2*nu)/2];

    otherwise
        error('Unsupported mechType: %s', mechType);
end

end