function [ m_e ] = elementMassMatrix( rho, A, l, type )
%elementMassMatrix: Returns either the consistent or lumped element mass
%matrix for bending and axial loads.
%   Determine if the consistent or lumped matrix is needed
    if strcmp(type,'consistent')
        C1 = (rho*A*l/420);
        C2 = (rho*A*l/6);
        m_e = [   2*C2,         0,          0,          C2,         0,          0;
                     0,    156*C1,    22*l*C1,           0,     54*C1,   -13*l*C1;
                     0,   22*l*C1,   4*l^2*C1,           0,   13*l*C1,  -3*l^2*C1;
                    C2,         0,          0,        2*C2,         0,          0;
                     0,     54*C1,    13*l*C1,           0,    156*C1,   -22*l*C1;
                     0,  -13*l*C1,  -3*l^2*C1,           0,  -22*l*C1,   4*l^2*C1   ];
    elseif strcmp(type,'lumped')
        C1 = (rho*A*l/420);
        C2 = (rho*A*l/6);
        m_e = [   2*C2,         0,          0,          C2,         0,          0;
                     0,        C1,          0,           0,         0,          0;
                     0,         0,          0,           0,         0,          0;
                    C2,         0,          0,        2*C2,         0,          0;
                     0,         0,          0,           0,        C1,          0;
                     0,         0,          0,           0,         0,          0   ];
    else
        disp('ERROR: Invalid mass matrix type. Valid types are ''consistent'' or ''lumped.''')
    end
end


