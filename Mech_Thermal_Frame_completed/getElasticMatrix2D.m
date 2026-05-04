function D = getElasticMatrix2D(E, nu, mechType)
% Elastic matrix for 2D linear isotropic elasticity
% strain = [exx eyy gamma_xy]'

if nargin < 3 || isempty(mechType)
    mechType = 'planeStress';
end

switch lower(mechType)
    case {'planestress','plane_stress'}
        D = E/(1 - nu^2) * [1   nu  0;
                            nu  1   0;
                            0   0   (1 - nu)/2];

    case {'planestrain','plane_strain'}
        D = E/((1 + nu)*(1 - 2*nu)) * [1 - nu   nu       0;
                                       nu       1 - nu   0;
                                       0        0        (1 - 2*nu)/2];

    otherwise
        error('Unsupported mechType: %s. Use planeStress or planeStrain.', mechType);
end

end
