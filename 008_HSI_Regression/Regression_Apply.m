function [CartoY, CartoYlabel,Rv_s, Ss, Yc]=Regression_Apply(M,RGB,d,wl,Model,figi,indx_Reg,Yref,dref,pas,pas2)
% Function to apply the model created with the function CreateModel.
%
% INPUT:
%       M: Hyperspectral datacube (n*m*p)
%       RGB: Associated RGB image (n*m*3)
%       d: Associated depth (1*m)
%       wl: Associated wavelengths (1*p)
%       Model: Structure containing the model parameters
%       figi: display figure with value different from 0
%       Optionnal:
%           Yref: Reference data to compare with the prediction map
%           dref: Depth of the reference sampling
%           pas: Width of the reference sampling
%           pas2: Depth limits of the reference sampling (Max, Min,
%                   Center)
% OUTPUT:
%       CartoY: Predicted abundance map
%       CartoYlabel: Name of the variable predicted
%       Rv_s (only if ref available): Correlation between predicted and
%           reference data
%       Ss (only if ref available): Standard deviation of the index
%       Yc (only if ref available): Subsampling predicted values at the
%           reference resolution
%
% This is the software from the papers:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Jacq, K., Perrette, Y., Fanget, B., Sabatier, P., Coquin, D.,
% Martinez-Lamas, R., Debret, M., Arnaud, F., 2019. High-resolution
% prediction of organic matter concentration with hyperspectral imaging
% on a sediment core. Sci. Total Environ. 663, 236–244.
% https://doi.org/10.1016/j.scitotenv.2019.01.320

% Please cite our papers if you use our code for your research.

% Verification of inputs
if nargin<6
    figi=0;
end
if nargin<7
    indx_Reg=0;
end

% Wavelength reduce
if isfield(Model,'idx_wl1')
    M=M(:,:,Model.idx_wl1);
    wl=wl(Model.idx_wl1);
end

m=median(M(:));
if m>1000
    M=M/100;
end

% Conversion in absorbance if necessary
if strcmp(Model.Spectral,'Absorbance')
    S=reshape(M,size(M,1)*size(M,2),size(M,3));
    [a,b]=find(S==0);
    S(a,b)=eps;
    S=log10(1./S);
else
    if strcmp(Model.Spectral,'Absorbance logN')
        S=reshape(M,size(M,1)*size(M,2),size(M,3));
        [a,b]=find(S==0);
        S(a,b)=eps;
        S=log(1./S);
    else
        S=reshape(M,size(M,1)*size(M,2),size(M,3));
    end
end

% Reduced wavelengths
if isfield(Model,'idx_wl2')
    idx_wl2=Model.idx_wl2;
else
    idx_wl2=zeros(1,length(Model.idx_wl1));
end

% Remove abberant pixels
% mask=AbberantPixels(M,RGB,d,wl);
mask=ones(size(M,1),size(M,2));

% First preprocessing
if strcmp(Model.Preprocessing,'Center')
    [Sp,~,~]=pretreat(S,'center',Model.PreprocessingPara.para1,Model.PreprocessingPara.para2);
else
    if strcmp(Model.Preprocessing,'Autoscaling')
        [Sp,~,~]=pretreat(S,'autoscaling',Model.PreprocessingPara.para1,Model.PreprocessingPara.para2);
    else
        Sp=S;
    end
end

% Second preprocessing
if strcmp(Model.SpectralProcessing,'Raw')==0
    pret={'Detrend','SNV','SNVD','MSC','D1','D2','SNV+D1','SNVD+D1','MSC+D1','SNV+D2','SNVD+D2','MSC+D2','CR'};
    [~,b]=find(double(strcmp(pret,Model.SpectralProcessing))==1);
    Sp = AllPret(Sp,wl,b);
    Pretidx=b+1;
else
    Pretidx=1;
end

% Regression methods
Reg={'ANN','PLS','PLS-BPNN','RT','RTE','MARS','SVM'};
Regidx=[];
RegPerf=[];
if strcmp(Model.Type,'Complete')
    if size(Model.ANN_Perf,2)==2
        if isfield(Model,'ANN_Perf')
            Regidx=[Regidx 1];
            RegPerf=[RegPerf Model.ANN_Perf{6, 2}];
        end
        if isfield(Model,'PLSR_Perf')
            Regidx=[Regidx 2];
            RegPerf=[RegPerf Model.PLSR_Perf{7, 2}];
        end
        if isfield(Model,'BPNN_Perf')
            Regidx=[Regidx 3];
            RegPerf=[RegPerf Model.BPNN_Perf{6, 2}];
        end
        if isfield(Model,'RT_Perf')
            Regidx=[Regidx 4];
            RegPerf=[RegPerf Model.RT_Perf{6, 2}];
        end
        if isfield(Model,'RTE_Perf')
            Regidx=[Regidx 5];
            RegPerf=[RegPerf Model.RTE_Perf{6, 2}];
        end
        if isfield(Model,'MARS_Perf')
            Regidx=[Regidx 6];
            RegPerf=[RegPerf Model.MARS_Perf{6, 2}];
        end
        if isfield(Model,'SVM_Perf')
            Regidx=[Regidx 7];
            RegPerf=[RegPerf Model.SVM_Perf{6, 2}];
        end
    else
        if isfield(Model,'ANN_Perf')
            Regidx=[Regidx 1];
            RegPerf=[RegPerf Model.ANN_Perf{1, 6}];
        end
        if isfield(Model,'PLSR_Perf')
            Regidx=[Regidx 2];
            RegPerf=[RegPerf Model.PLSR_Perf{1, 8}];
        end
        if isfield(Model,'BPNN_Perf')
            Regidx=[Regidx 3];
            RegPerf=[RegPerf Model.BPNN_Perf{1, 6}];
        end
        if isfield(Model,'RT_Perf')
            Regidx=[Regidx 4];
            RegPerf=[RegPerf Model.RT_Perf{1, 6}];
        end
        if isfield(Model,'RTE_Perf')
            Regidx=[Regidx 5];
            RegPerf=[RegPerf Model.RTE_Perf{1,6}];
        end
        if isfield(Model,'MARS_Perf')
            Regidx=[Regidx 6];
            RegPerf=[RegPerf Model.MARS_Perf{1,6}];
        end
        if isfield(Model,'SVM_Perf')
            Regidx=[Regidx 7];
            RegPerf=[RegPerf Model.SVM_Perf{1,6}];
        end
    end
else
    if isfield(Model,'ANN_Perf')
        Regidx=[Regidx 1];
        RegPerf=[RegPerf Model.ANN_Perf.R2val];
    end
    if isfield(Model,'PLSR_Perf')
        Regidx=[Regidx 2];
        RegPerf=[RegPerf Model.PLSR_Perf.R2val];
    end
    if isfield(Model,'BPNN_Perf')
        Regidx=[Regidx 3];
        RegPerf=[RegPerf Model.BPNN_Perf.R2val];
    end
    if isfield(Model,'RT_Perf')
        Regidx=[Regidx 4];
        RegPerf=[RegPerf Model.RT_Perf.R2val];
    end
    if isfield(Model,'RTE_Perf')
        Regidx=[Regidx 5];
        RegPerf=[RegPerf Model.RTE_Perf.R2val];
    end
    if isfield(Model,'MARS_Perf')
        Regidx=[Regidx 6];
        RegPerf=[RegPerf Model.MARS_Perf.R2val];
    end
    if isfield(Model,'SVM_Perf')
        Regidx=[Regidx 7];
        RegPerf=[RegPerf Model.SVM_Perf.R2val];
    end
end

iter=0;
for i=1:7
    if sum(Regidx==i)>0
        iter=iter+1;
        list{iter} = strcat(Reg{i},' : R²val=',num2str(round(RegPerf(iter),2)));
    end
end
if indx_Reg==0
    [indx_Reg,~] = listdlg('ListString',list);
end

%% Prediction
Pret={'Raw','Detrend','SNV','SNVD','MSC','D1','D2','SNV+D1','SNVD+D1','MSC+D1','SNV+D2','SNVD+D2','MSC+D2','CR'};
iter=0;
if sum(Regidx(indx_Reg)==1)>0 % ANN
    iter=iter+1;
    if strcmp(Model.Type,'Complete')
        if size(Model.ANN_Perf,2)==2
            Yp{iter} = Model.ANN(Sp(:,idx_wl2==0)')';
        else
            Yp{iter} = Model.ANN{1,1}(Sp(:,idx_wl2==0)')';
        end
    else
        Spr=Sp(:,idx_wl2==0);
        Yp{iter} = Model.ANN(Spr(:,Model.ANN_Perf.idx_vs)')';
    end
end
clear Spret prettmp

if sum(Regidx(indx_Reg)==2)>0 % PLSR
    iter=iter+1;
    if strcmp(Model.Type,'Complete')
        if size(Model.PLSR_Perf,2)==2
            Yp{iter}=Uncenterval(Centerval(Sp(:,idx_wl2==0),Model.PLSR.Pred.Scal)*Model.PLSR.para.B,Model.PLSR.Pred.Ycal);
        else
            Yp{iter}=Uncenterval(Centerval(Sp(:,idx_wl2==0),Model.PLSR{1,1}.Pred.Scal)*Model.PLSR{1,1}.para.B,Model.PLSR{1,1}.Pred.Ycal);
        end
    else
        Spr=Sp(:,idx_wl2==0);
        Yp{iter}=Uncenterval(Centerval(Spr(:,Model.PLSR_Perf.idx_vs),Model.PLSR.Pred.Scal)*Model.PLSR.para.B,Model.PLSR.Pred.Ycal);
    end
end
clear Spret prettmp

if sum(Regidx(indx_Reg)==3)>0 % BPNN
    iter=iter+1;
    if strcmp(Model.Type,'Complete')
        if size(Model.BPNN_Perf,2)==2
            Yp{iter} = Model.BPNN((Sp(:,idx_wl2==0)*Model.BPNN_Perf.PLSRloadings)')';
        else
            Yp{iter} = Model.BPNN{1,1}((Sp(:,idx_wl2==0)*Model.BPNNplsloadings)')';
        end
    else
        Spr=Sp(:,idx_wl2==0);
        Yp{iter} = Model.BPNN((Spr(:,Model.BPNN_Perf.idx_vs)*Model.BPNN_Perf.PLSRloadings)')';
    end
end
clear Spret prettmp

if sum(Regidx(indx_Reg)==4)>0 % RT
    iter=iter+1;
    if strcmp(Model.Type,'Complete')
        if size(Model.RT_Perf,2)==2
            Yp{iter} = predict(Model.RT,Sp(:,idx_wl2==0));
        else
            Yp{iter} = predict(Model.RT{1, 1},Sp(:,idx_wl2==0));
        end
    else
        Spr=Sp(:,idx_wl2==0);
        Yp{iter} = predict(Model.RT,Spr(:,Model.RT_Perf.idx_vs));
    end
end
clear Spret prettmp

if sum(Regidx(indx_Reg)==5)>0 % RTE
    iter=iter+1;
    if strcmp(Model.Type,'Complete')
        if size(Model.RTE_Perf,2)==2
            Yp{iter} = predict(Model.RTE,Sp(:,idx_wl2==0));
        else
            Yp{iter} = predict(Model.RTE{1, 1},Sp(:,idx_wl2==0));
        end
    else
        Spr=Sp(:,idx_wl2==0);
        Yp{iter} = predict(Model.RTE,Spr(:,Model.RTE_Perf.idx_vs));
    end
end
clear Spret prettmp

if sum(Regidx(indx_Reg)==6)>0 % MARS
    iter=iter+1;
    if strcmp(Model.Type,'Complete')
        if size(Model.MARS_Perf,2)==2
            [Yp{iter}, ~] = arespredict(Model.MARS, Sp(:,idx_wl2==0));
        else
            [Yp{iter}, ~] = arespredict(Model.MARS{1, 1}, Sp(:,idx_wl2==0));
        end
    else
        Spr=Sp(:,idx_wl2==0);
        [Yp{iter}, ~] = arespredict(Model.MARS, Spr(:,Model.MARS_Perf.idx_vs));
    end
end
clear Spret prettmp

if sum(Regidx(indx_Reg)==7)>0 % SVM
    iter=iter+1;
    if strcmp(Model.Type,'Complete')
        if size(Model.SVM_Perf,2)==2
            Yp{iter} = predict(Model.SVM,Sp(:,idx_wl2==0));
        else
            Yp{iter} = predict(Model.SVM{1, 1},Sp(:,idx_wl2==0));
        end
    else
        Spr=Sp(:,idx_wl2==0);
        Yp{iter} = predict(Model.SVM,Spr(:,Model.SVM_Perf.idx_vs));
    end
end
clear Spret prettmp

% Folding of the prediction
for i=1:length(Yp)
    CartoY{i}=reshape(Yp{i},size(M,1),size(M,2),size(Yp{i},2)).*mask;
    [a,b]=find(CartoY{i}==inf);
    for j=1:length(a)
        CartoY{i}(a(j),b(j))=NaN;
    end
    [a,b]=find(CartoY{i}==-inf);
    for j=1:length(a)
        CartoY{i}(a(j),b(j))=NaN;
    end
    if size(Yp{i},2)==1
        if isstruct(Model.Proxy)
            CartoYlabel{i}=strcat(Model.Proxy{1,1},{' '},Reg{Regidx(indx_Reg(i))},{' '},Pret{Pretidx(i)});
        else
            CartoYlabel{i}=strcat(Model.Proxy,{' '},Reg{Regidx(indx_Reg(i))},{' '},Pret{Pretidx});
        end
    else
        CartoYlabel{i}=Model.Proxy(i);
    end
end

% Saving
assignin('base', 'CartoY', CartoY);
assignin('base', 'CartoYmlabel', CartoYlabel);

% Display
if figi>0
    for k=1:length(CartoY)
        if nargin>7
            figure
            CartoYi=CartoY{k};
            ha(1)=subplot(3,1,1);
            imagesc(d,d(1:size(CartoYi,1)),CartoYi);
            colormap([0 0 0;jet]);
            colorbar,
            title('Carte d''abondance')
            ylabel('Largeur (cm)')
            caxis([nanmean(CartoYi(:))-3*nanstd(CartoYi(:)) nanmean(CartoYi(:))+3*nanstd(CartoYi(:))])
            set(findobj('type','axes'),'fontsize',14)
            xlim([min(d) max(d)])
            title(CartoYlabel{k})
            
            ha(2)=subplot(3,1,2);
            imagesc(d,d(1:size(CartoYi,1)),RGB*(0.5/mean(RGB(:))))
            ylabel('Largeur (cm)')
            colorbar
            xlim([min(d) max(d)])
            set(findobj('type','axes'),'fontsize',14)
            
            ha(3)=subplot(3,1,3);
            if median(CartoYi(:))+3*std(CartoYi(:))>median(Yref)&&median(CartoYi(:))-3*std(CartoYi(:))<median(Yref)
            else
                yyaxis left
            end
            Cartomin=zeros(1,size(CartoY,2));
            Cartomax=zeros(1,size(CartoY,2));
            for i=1:size(CartoYi,2)
                if i==874
                    a=1;
                end
                CartoYii= CartoYi(:,i);
                if sum(isnan(CartoYii))<size(CartoYi,1)-10
                    Cartomax(i)=nanmax(CartoYii(CartoYii<nanmean(CartoYii)+2*nanstd(CartoYii)));
                    Cartomin(i)=nanmin(CartoYii(CartoYii>nanmean(CartoYii)-2*nanstd(CartoYii)));
                else
                    Cartomax(i)=NaN;
                    Cartomin(i)=NaN;
                end
            end
            h=area(d,[Cartomin; Cartomax-Cartomin]','LineStyle','none');
            h(1).FaceColor = [1 1 1];
            h(2).FaceColor = [0 1 0];
            if nargin>7
                if median(CartoYi(:))+3*std(CartoYi(:))>median(Yref)&&median(CartoYi(:))-3*std(CartoYi(:))<median(Yref)
                    hold on
                else
                    yyaxis right
                end
                if length(pas)==1
                    pas=repmat(pas,1,length(dref));
                end
                for i=1:length(dref)
                    hold on
                    if strcmp(pas2,'Max')
                        h1=plot([dref(i)-pas(i) dref(i)],[Yref(i) Yref(i)],'-','LineWidth',3,'color',[1 0.66 0]);
                    else
                        if strcmp(pas2,'Min')
                            h1=plot([dref(i) dref(i)+pas(i)],[Yref(i) Yref(i)],'-','LineWidth',3,'color',[1 0.66 0]);
                        else
                            if strcmp(pas2,'Center')
                                h1=plot([dref(i)-0.5*pas(i) dref(i)+0.5*pas(i)],[Yref(i) Yref(i)],'-','LineWidth',3,'color',[1 0.66 0]);
                            end
                        end
                    end
                end
                if median(CartoYi(:))+3*std(CartoYi(:))>median(Yref)&&median(CartoYi(:))-3*std(CartoYi(:))<median(Yref)
                    ylim([nanmin([nanmean(Cartomin(:))-3*nanstd(Cartomin(:)) Yref']) max([nanmean(Cartomax(:))+3*nanstd(Cartomax(:)) Yref'])])
                else
                    yyaxis left
                    ylim([nanmean(Cartomin(:))-3*nanstd(Cartomin(:)) nanmean(Cartomax(:))+3*nanstd(Cartomax(:))])
                end
            else
                ylim([nanmin([nanmean(Cartomin(:))-3*nanstd(Cartomin(:)) Yref']) max([nanmean(Cartomax(:))+3*nanstd(Cartomax(:)) Yref'])])
            end
            set(gca,'Layer','top')
            if isstruct(Model.Proxy)
                ylabel(Model.Proxy{1,1})
            else
                ylabel(Model.Proxy)
            end
            xlim([min(d) max(d)])
            legend([h(2) h1],'Prédit','Observé')
            grid on
            set(findobj('type','axes'),'fontsize',14)
            
            linkaxes(ha,'x')
        else
            figure
            CartoYi=CartoY{k};
            ha(1)=subplot(3,1,1);
            imagesc(d,d(1:size(RGB,1)),CartoYi);
            colormap(jet);
            colorbar,
            xlim([min(d) max(d)])
            title('Carte d''abondance')
            xlabel('Profondeur (cm)')
            ylabel('Largeur (cm)')
            [a,b]=find(CartoYi==inf);
            for j=1:length(a)
                CartoYi(a(j),b(j))=NaN;
            end
            [a,b]=find(CartoYi==-inf);
            for j=1:length(a)
                CartoYi(a(j),b(j))=NaN;
            end
            caxis([nanmean(CartoYi(:))-3*nanstd(CartoYi(:)) nanmean(CartoYi(:))+3*nanstd(CartoYi(:))])
            set(findobj('type','axes'),'fontsize',14)
            title(CartoYlabel{k})
            
            ha(2)=subplot(3,1,2);
            imagesc(d,d(1:size(RGB,1)),RGB*(0.5/mean(RGB(:))))
            xlabel('Profondeur (cm)')
            ylabel('Largeur (cm)')
            colorbar
            xlim([min(d) max(d)])
            set(findobj('type','axes'),'fontsize',14)
            
            ha(3)=subplot(3,1,3);
            Cartomin=zeros(1,size(CartoY{1,1},2));Cartomax=zeros(1,size(CartoY{1,1},2));
            for i=1:size(CartoYi,2)
                CartoYii= CartoYi(:,i);
                if isnan(nanmax(CartoYii(CartoYii<nanmean(CartoYii)+2*nanstd(CartoYii))))==0
                    Cartomax(i)=nanmax(CartoYii(CartoYii<nanmean(CartoYii)+2*nanstd(CartoYii)));
                    Cartomin(i)=nanmin(CartoYii(CartoYii>nanmean(CartoYii)-2*nanstd(CartoYii)));
                end
            end
            h=area(d,[Cartomin; Cartomax-Cartomin]','LineStyle','none');
            h(1).FaceColor = [1 1 1];
            h(2).FaceColor = [0 1 0];
            set(gca,'Layer','top')
            xlabel('Profondeur (cm)')
            %                 ylabel(Model.Proxy{1,1})
            xlim([min(d) max(d)])
            ylim([nanmin(nanmean(Cartomin(:))-3*nanstd(Cartomin(:))) nanmax(nanmean(Cartomax(:))+3*nanstd(Cartomax(:)))])
            legend(h(2),'Prédit')
            grid on
            colorbar
            set(findobj('type','axes'),'fontsize',14)
            
            linkaxes(ha,'x')
        end
    end
    
    figure
    ha(1)=subplot(211);
    imagesc(d,d(1:size(RGB,1)),RGB*(0.5/mean(RGB(:))))
    xlabel('Profondeur (cm)')
    ylabel('Largeur (cm)')
    xlim([min(d) max(d)])
    set(gca,'fontsize',14)
    
    for k=1:length(CartoY)
        ha(2)=subplot(212);
        if nargin>7&&median(CartoYi(:))+3*std(CartoYi(:))>median(Yref)&&median(CartoYi(:))-3*std(CartoYi(:))<median(Yref)
        else
            yyaxis left
        end
        plot(d,nanmedian(CartoY{k}),'linewidth',2)
        %             hold on
        if isstruct(Model.Proxy)
            ylabel(strcat(Model.Proxy{1,1},' predicted'))
        else
            ylabel(strcat(Model.Proxy,' predicted'))
        end
    end
    if nargin>7
        if median(CartoYi(:))+3*std(CartoYi(:))>median(Yref)&&median(CartoYi(:))-3*std(CartoYi(:))<median(Yref)
            hold on
        else
            yyaxis right
        end
        plot(dref,Yref,'.','markersize',20)
        if isstruct(Model.Proxy)
            ylabel(strcat(Model.Proxy{1,1},' observed'))
        else
            ylabel(strcat(Model.Proxy,' observed'))
        end
    end
    grid on
    xlabel('Profondeur (cm)')
    
    xlim([min(d) max(d)])
    legend(Reg{indx_Reg})
    linkaxes(ha,'x')
    set(gca,'fontsize',14)
end

% Correlation with reference values
if nargin>7
    for k=1:length(CartoY)
        [Rv_s(k),Ss(k), Yc{k}] = CorrelCartoY(CartoY{k},d,Yref,dref,figi,pas,pas2);
    end
else
    Rv_s=[];
    Yc=[];
end

end