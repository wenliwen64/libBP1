 function openmfiles(project)
% clear all
% project='tmp'
load(['~/wk/matlab/projectlist/' project '.mat']);% Load the previously saved mat file
Files{:}
matlab.desktop.editor.openDocument(Files) % % Open all the file previously saved in the MATLAB editor