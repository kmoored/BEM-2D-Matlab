% Flight trajectory: Translation
if free == 1
    Q0 = kron(-[Uinf 0]',ones(1,Ncyc*Nstep+1));     % Initializing the body velocity
else
    Q0 = kron(-[Uinf Winf]',ones(1,Ncyc*Nstep+1));  % Initializing the body velocity
end
x_b = zeros(1,3);
z_b = D*ones(1,3);
a_b = 0;