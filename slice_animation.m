function slice_animation(casename,run,nworkers)
parpool(nworkers);
hfcase = DNS_case(casename, run);
%hfcase = DNS_case('cwl90_window_turb_clean',3);
%runs =[5];
%hfcase.readKSlices([],61:85);
imgfolder = fullfile(hfcase.casepath,sprintf('run%d',run),'animation_images','k_Mach');
%%
if ~exist(imgfolder, 'dir')
       mkdir(imgfolder);
end

%%
parfor i=1:hfcase.nSlices
    fprintf('Plotting slice %d/%d\n',i,hfcase.nSlices)
    slice = hfcase.readSingleKSlice(i);
    slice2kPlot(slice, 'M', fullfile(imgfolder,sprintf('img_%03d.png',i)), [0 1.4], 'M');
end

end
%%

% clear
% close all
%
% 
% hfcase = DNS_case('cwl90_tripped2',1);
% %runs =[5];
% %hfcase.readKSlices([],61:85);
% imgfolder = fullfile(hfcase.casepath,'run1','animation_images','k_Mach');
% %%
% if ~exist(imgfolder, 'dir')
%        mkdir(imgfolder);
% end
% 
% h = figure;
% h.Visible = 'off';
% ax = gca;
% 
% %%
% for i=11:hfcase.nSlices
%     fprintf('Plotting slice %d/%d\n',i,hfcase.nSlices)
%     slice = hfcase.readSingleKSlice(i);
%     hfcase.kPlot(slice,'M',ax,[0 1.4],'M')
%     %hfcase.jPlot(hfcase.jSlices(i),'tau_w',[],[0 800],'\tau_w')
%     exportgraphics(gcf,fullfile(imgfolder,sprintf('img_%03d.png',i)),'Resolution',600);
%     clear slice
% end
% %%
