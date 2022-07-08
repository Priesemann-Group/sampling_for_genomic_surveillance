%%% main plot detection time %%%

clear all
close all
clc

%% Preload other parameters

load('Fig_1_Experiment_1_K_25_mode_1.mat')
load('DefColors.mat')
load('Colores.mat')
labels_resorted = {'v. 2 (const)','v. 2 (adapt)',...
    'v. 3 (const)','v. 3 (adapt)',...
    'v. 4 (const)','v. 4 (adapt)',...
    'v. 5 (const)','v. 5 (adapt)'};

id_imp = 1:2:2*(k-1);
id_imp = [1+id_imp; id_imp];
Def = zeros(2*(k-1),3);
Def(id_imp(1,:),:) = Default(2:5,:);
Def(id_imp(2,:),:) = Default(2:5,:);

%% Defining empty matrices and cells for extracting results

Matrix = cell(4,3);
Tin_matrix = cell(4,3);
Kindice = [25 50 100];
Mode = [0 1];

for expno = 1:4
    for i = 1:3
        P.K = Kindice(i);
        Matrix{expno,i} = nan(m_max,2*(k-1));
        Tin_matrix{expno,i} = nan(1,k-1);
        for j = 1:2
            mode = Mode(j);
            load(strcat('Fig_1_Experiment_', num2str(expno),'_K_',num2str(P.K),'_mode_',num2str(mode),'.mat'),'Tdet_var','Tin')
            Matrix{expno,i}(:,id_imp(j,:)) = Tdet_var(:,2:end);
        end
        Tin_matrix{expno,i} = Tin;
    end
end

h = figure('units','centimeters','position',[3,3,18,16]);

for expno = 1:4
    for Kind = 1:3
        sp_id = Kind + 3*(expno-1);
        ax = subplot(4,3,sp_id);
        Tin = Tin_matrix{expno,1};
        M = Matrix{expno,Kind};
        for i = 1:length(Tin)
            plot([Tin(i) Tin(i)],[2*i 2*(i+1)],'--','Color',Default(i+1,:),'LineWidth',1.5)
            hold on
        end
        a = boxplot(M,labels_resorted,'PlotStyle','compact','LabelOrientation','horizontal','Orientation','horizontal','Color',Def);
        set(gca,'YTick',[])
        set(gca,'YTickLabels',[])
        set(gca,'XTick',0:10:60)
        set(a,{'linew'},{1.5})
        set(gca,'FontSize',8,'XColor','k','YColor','k','TickLength',[0.025 0.025])
        xlim([0 55])
    end
end

% print(h,'Figure_1_results','-dpdf')
