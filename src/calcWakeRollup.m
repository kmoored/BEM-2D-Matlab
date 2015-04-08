%% Wake rollup calculations
if (Rollup == 1 && i_t > 1)

    if i_t > Nlump*Nstep + 1
        lump = 1;
    else
        lump = [];
    end

    [u,w] = WakeRollupVelocity2D([xw(2:wakeInd+1) lump*xl(1)],[zw(2:wakeInd+1) lump*zl(1)],mu(:,1)',muTE(1),[muW(1:wakeInd) lump*muLump(1)],muLEt,muLEb,sigma',xp',zp',vt',vn',xTE',zTE',vtTE',vnTE',[xw(1:wakeInd+1) lump*xl(1)],[zw(1:wakeInd+1) lump*zl(1)],[vtw(1:wakeInd,:); lump*vtlw(1) lump*vtlw(2)]',[vnw(1:wakeInd,:); lump*vnlw(1) lump*vnlw(2)]',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);

    u_2 = 0*u;
    w_2 = 0*w;

    if grd == 1
        [u_2,w_2] = WakeRollupVelocity2D([xw(2:wakeInd+1) lump*xl(1)],[zw(2:wakeInd+1) lump*zl(1)],mu(:,1)',muTE(1),[muW(1:wakeInd) lump*muLump(1)],muLEt,muLEb,sigma',xp',-zp',[vt(:,1) -vt(:,2)]',[vn(:,1) -vn(:,2)]',xTE',-zTE',[vtTE(1,1) -vtTE(1,2)]',[vnTE(1,1) -vnTE(1,2)]',[xw(1:wakeInd+1) lump*xl(1)],[-zw(1:wakeInd+1) -lump*zl(1)],[[vtw(1:wakeInd,1) -vtw(1:wakeInd,2)]; lump*vtlw(1) -lump*vtlw(2)]',[[vnw(1:wakeInd,1) -vnw(1:wakeInd,2)]; lump*vnlw(1) -lump*vnlw(2)]',xpt_LES,zpt_LES,vtLEt,vnLEt,xpb_LES,zpb_LES,vtLEb,vnLEb,i_t,Ncyc,Nstep,LES,ep,epBod,SC);
    end

    % Applying fencing scheme
    [xstar_w,zstar_w,fence] = fencing(c,Vp,[xw(2:wakeInd+1) lump*xl(1)],[zw(2:wakeInd+1) lump*zl(1)],u + u_2,w + w_2,xp,zp,vn,x_b(2),z_b(2),delT,1/100*epBod);
    xstar_w = xstar_w.*fence + (xstar_w + (u + u_2)*delT).*(~fence);
    zstar_w = zstar_w.*fence + (zstar_w + (w + w_2)*delT).*(~fence);

    % Wake rollup
    if i_t > Nlump*Nstep + 1
        xw(2:wakeInd+1) = xstar_w(1:end-1);
        zw(2:wakeInd+1) = zstar_w(1:end-1);

        vttemp = [(xw(1,2:end) - xw(1,1:end-1))' (zw(1,2:end) - zw(1,1:end-1))'];
        vtw = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
        vnw = [-vtw(:,2) vtw(:,1)]; 

        xl(1) = xstar_w(end);
        zl(1) = zstar_w(end);

        vttemp = [(xl(1) - xw(Nlump*Nstep+1))' (zl(1) - zw(Nlump*Nstep+1))'];
        vtlw(1,:) = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
        vnlw(1,:) = [-vtlw(1,2) vtlw(1,1)]; 
    else
        xw(2:wakeInd+1) = xstar_w;
        zw(2:wakeInd+1) = zstar_w;

        vttemp = [(xw(1,2:end) - xw(1,1:end-1))' (zw(1,2:end) - zw(1,1:end-1))'];
        vtw = diag(1./sqrt(vttemp(:,1).^2 + vttemp(:,2).^2))*vttemp;
        vnw = [-vtw(:,2) vtw(:,1)];
    end    

end
