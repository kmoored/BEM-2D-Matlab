if SaveData == 1
    savefilename = ['_Scratch_f',num2str(f),...
            '_A_c',num2str(A_c),...
            '_d_c',num2str(d_c)];
        
    save(['scratch/Parameters',savefilename,'.mat'],'-v7.3')

    fid_Data = fopen(['scratch/Data',savefilename,'.txt'],'w');
    fid_PanelProp = fopen(['scratch/PanelProp',savefilename,'.txt'],'w');
    fid_WakeProp = fopen(['scratch/WakeProp',savefilename,'.txt'],'w');
end
