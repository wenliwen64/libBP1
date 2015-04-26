function savemfiles(project)
allDocs = matlab.desktop.editor.getAll; %Save all opened documents to allDocs
Files = {allDocs.Filename};	%Save all the path of opened documents to Files
['~/wk/matlab/projectlist/' project '.mat']
save(['~/wk/matlab/projectlist/' project '.mat'],'Files');	 % Save the variable Files into filename.mat file