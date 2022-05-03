clear
close all

hfcase = DNS_case('cwl90_window_turb_clean',3);
%runs =[5];
hfcase.readKSlices([],61:85);
imgfolder = fullfile(hfcase.casepath,['run3','animation_images','k_Mach');
%%
if ~exist(imgfolder, 'dir')
       mkdir(imgfolder);
end

h = figure;
h.Visible = 'off';
ax = gca;

%%
for i=1:hfcase.nSlices
    fprintf('Plotting slice %d/%d\n',i,hfcase.nSlices)\
    slice = hfcase.readSingleKSlice(i);
    hfcase.kPlot(slice,'M',ax,[0 1.4],'M')
    %hfcase.jPlot(hfcase.jSlices(i),'tau_w',[],[0 800],'\tau_w')
    exportgraphics(gcf,fullfile(imgfolder,sprintf('img_%03d.png',i)),'Resolution',600);
end
%%
