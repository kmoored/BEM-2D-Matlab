function [  ] = progTitle( version )
    fprintf('/*---------------------------------------------------------------------------*\\\n')
    fprintf('| .   )\\      B oundary     |                                                 |\n')
    fprintf('| \\`.-'' `-oo  E lement      | Written by:  Dr. Keith Moored                   |\n')
    fprintf('|  ) _  __,0) F luid        | Version:     %s                              |\n',version)
    fprintf('| /.'' )/      S tructure    | Web:         http://www.lehigh.edu/             |\n')  
    fprintf('| ''           I nteraction  |                                                 |\n')
    fprintf('\\*---------------------------------------------------------------------------*/\n')
    [~, name] = system('hostname');
    date = datestr(datetime(clock,'Format','MMMM, d y'),'mmmm dd, yyyy');
    clockTime = datestr(datetime(clock,'Format','HH:mm:ss'),'HH:MM:SS.FFF AM');
    fprintf('Date     : %s\n', date)
    fprintf('Time     : %s\n', clockTime)
    fprintf('Host     : %s\n', name(1:end-1))
    fprintf('PID      : %i\n', feature('getpid'))
    fprintf('Case     : %s\n', pwd)
    fprintf('\n')
    fprintf('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n')
    fprintf('Calculated inputs:\n')
end

