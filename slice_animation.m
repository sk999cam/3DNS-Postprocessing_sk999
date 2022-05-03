clear
close all

hfcase = DNS_case('r150_cwl90_hf2',5);
runs =[5];
hfcase.readKSlices([],61:85);
imgfolder = fullfile(hfcase.casepath,['run' num2str(runs(end))],'animation_images','k_Mach');
%%
if ~exist(imgfolder, 'dir')
       mkdir(imgfolder);
end

h = figure;
h.Visible = 'off';
ax = gca;

%%
for i=61:85%:hfcase.nSlices]
    fprintf('Plotting slice %d/%d\n',i,hfcase.nSlices)
    hfcase.kPlot(hfcase.kSlices(i),'M',ax,[0 1.4],'M')
    %hfcase.jPlot(hfcase.jSlices(i),'tau_w',[],[0 800],'\tau_w')
    exportgraphics(gcf,fullfile(imgfolder,sprintf('img_%03d.png',i)),'Resolution',400);
end
%%
