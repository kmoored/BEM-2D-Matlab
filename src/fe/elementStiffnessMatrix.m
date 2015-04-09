function [ k_e ] = elementStiffnessMatrix( E, I, A, l )
%elementStiffnessMatrix: Returns the element stiffness matrix for bending
%and axial loads.
    C1 = (E*A / l);
    C2 = (E*I/l^3);
    k_e =  [  C1,          0,           0,         -C1,           0,           0;
               0,      12*C2,      6*l*C2,           0,      -12*C2,      6*l*C2;
               0,     6*l*C2,    4*l^2*C2,           0,     -6*l*C2,    2*l^2*C2;
             -C1,          0,           0,          C1,           0,           0;
               0,     -12*C2,     -6*l*C2,           0,       12*C2,     -6*l*C2;
               0,     6*l*C2,    2*l^2*C2,           0,     -6*l*C2,     4*l^2*C2     ];
end

