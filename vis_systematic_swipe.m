%%% visualization of systematic validation (heatmaps) %%%

clear all
close all
clc

%%

mode = 1;
Kexp = 10;
filename = strcat('syst_swipe_mode',num2str(mode),'K',num2str(Kexp),'_250.mat');
load(filename,'D_com','R0max_swipe_ref','Tin_swipe_ref')

[I,J] = size(D_com);

D_av = nan(I,J);
for i = 1:I
    for j = 1:J
        D_av(i,j) = nanmean(D_com{i,j});
    end
end

h = figure('units','centimeters','position',[3,3,8,8]);
aux_D = D_av';
aux_D1  = flip(aux_D,1);
heatmap(round(R0max_swipe_ref,2),round(flip(Tin_swipe_ref),1),aux_D1)
xlabel('Rel. rep number second variant')
ylabel('Int. time second variant (relative to peak) [weeks]')
caxis([1, 10]);
set(gca,'FontSize',8)

%% calculo de diferencia con el constant protocol

aux_D1_mod1 = aux_D1;
mode = 0;
filename = strcat('syst_swipe_mode',num2str(mode),'K',num2str(Kexp),'_250.mat');
load(filename,'D_com','R0max_swipe_ref','Tin_swipe_ref')

[I,J] = size(D_com);

D_av = nan(I,J);
for i = 1:I
    for j = 1:J
        D_av(i,j) = nanmean(D_com{i,j});
    end
end
aux_D = D_av';
aux_D2  = flip(aux_D,1);
h = figure('units','centimeters','position',[3,3,8,8]);

white = [1 1 1];
red = [1 0 0];
n_col = 64;
phi = linspace(0,1,n_col)';
map = phi.*repmat(white,n_col,1) + (1-phi).*repmat(red,n_col,1);

heatmap(round(R0max_swipe_ref,2),round(flip(Tin_swipe_ref),1),aux_D1-aux_D2,'Colormap',map)
xlabel('Rel. rep number second variant')
ylabel('Int. time second variant (relative to peak) [weeks]')
caxis([-3, 0]);
set(gca,'FontSize',8)


