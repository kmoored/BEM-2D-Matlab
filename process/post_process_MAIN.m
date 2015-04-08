clear
close all
clc

%%%% Source paths on local machine
addpath('../src/geom')
addpath('../src/kine')
addpath('../src/util')
addpath('./util')

% Parameters for identifying data file.
f = 1;
A_c = 0.25;
d_c = 0;

savefilename = ['_scratch_f',num2str(f),...
            '_A_c',num2str(A_c),...
            '_d_c',num2str(d_c)];
folder = 'scratch';

%% Process data to a single .mat file
process_data(folder,savefilename)


%% Plot Flowfield data
process_flowfield(folder,savefilename)


%% Plot Time-Averaged Flowfield data
process_TAflowfield(folder,savefilename)


