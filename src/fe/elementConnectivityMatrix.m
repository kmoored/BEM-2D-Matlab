function [ l_e ] = elementConnectivityMatrix( Nelements, element, theta )
%elementConnectivityMatrix: Returns the element connetivity matrix for a
%beam element. The numbering convention is from left to right.
    l_e = zeros(6,3*(Nelements+1));
    C = cos(theta);
    S = sin(theta);
    temp = [   C,  S,  0,  0,  0,  0;
              -S,  C,  0,  0,  0,  0;
               0,  0,  1,  0,  0,  0;
               0,  0,  0,  C,  S,  0;
               0,  0,  0, -S,  C,  0;
               0,  0,  0,  0,  0,  1    ];
    l_e(:,3*element-2:5+3*element-2) = temp;
end

