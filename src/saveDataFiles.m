if SaveData == 1           
    Data = [wakeInd,Q0(1,i_t),Q0(2,i_t),x_b(2),z_b(2),a_b,Fx,Fz,Pow,L,T,D_visc,Gamma,Cf,Cl,Ct,Cpow,...
        Fx_s,Fx_us,Fz_s,Fz_us,Pow_s,Pow_us,L_s,L_us,T_s,T_us,Cl_s,Cl_us,...
        Ct_s,Ct_us,Cpow_s,Cpow_us];
    fprintf(fid_Data,'%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n',Data');

    PanelProp = [xp_0(1:end-1)'; zp_0(1:end-1)';xp(1:end-1)';zp(1:end-1)';...
        Vc(:,1)';Vc(:,2)';xc';zc';vt(:,1)';vt(:,2)';vn(:,1)';vn(:,2)';...
        dL';Xc';Zc';sigma';mu(:,1)';Qp';Qt';Cp_s';Cp_us';Cp';dFshear']';
    fprintf(fid_PanelProp,'%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n',PanelProp');

    WakeProp = [xw xl(1); zw zl(1);vtTE' vtw' vtlw(1,:)';...
        vnTE' vnw' vnlw(1,:)'; muTE(1) muW muLump(1); 
        GammaW -muLump(1)]';
    fprintf(fid_WakeProp,'%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n',WakeProp');
end
