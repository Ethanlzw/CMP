function D = getElasticMatrix2D(E, nu, mechType)
% getElasticMatrix2D  Return 3×3 constitutive matrix for 2D elasticity.
%
%   mechType: 'planeStress' | 'planeStrain'  (case-insensitive)
%
%   Stress-strain relation:  [sxx; syy; sxy] = D * [exx; eyy; gxy]

switch lower(mechType)
    case 'planestress'
        c = E / (1 - nu^2);
        D = c * [1,  nu,          0;
                 nu, 1,           0;
                 0,  0,  (1-nu)/2];

    case 'planestrain'
        c = E / ((1 + nu) * (1 - 2*nu));
        D = c * [1-nu, nu,       0;
                 nu,   1-nu,     0;
                 0,    0,  (1-2*nu)/2];

    otherwise
        error('getElasticMatrix2D: unknown mechType ''%s''. Use planeStress or planeStrain.', mechType);
end

end
