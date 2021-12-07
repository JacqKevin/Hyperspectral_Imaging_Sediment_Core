function Regression_Create(M,d,wl,Y,dy,Yn,svm)
% Function to create a PLSR model to predict a destructive variable with a
% hyperpectral image.
% INPUT :
%           M: Hyperspectral datacube (n*m*p)
%           d: Associated depth (1*m)
%           wl: Associated wavelengths (1*p)
%           Y: Reference values
%           dy: Reference depth
%           Yn: Name of the variable
%           Optionnal:
%           	svm: Variable selection method (0: Combination of all methods,
%                   1: VIP, 2: UVE, 3: CARS, 4: RF, 5: VCN)
% This is the software from the papers:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Jacq, K., Perrette, Y., Fanget, B., Sabatier, P., Coquin, D.,
% Martinez-Lamas, R., Debret, M., Arnaud, F., 2019. High-resolution
% prediction of organic matter concentration with hyperspectral imaging
% on a sediment core. Sci. Total Environ. 663, 236–244.
% https://doi.org/10.1016/j.scitotenv.2019.01.320

% Please cite our papers if you use our code for your research.

%% Input verification
if size(M,2)~=length(d)
    error('M and dm do not have the same size')
end
if size(M,3)~=length(wl)
    error('M and wl do not have the same size')
end

if size(Y,1)~=length(dy)
    error('Y and dy do not have the same dimension')
end
if iscell(Yn)||isnumeric(Yn)
    if size(Y,2)~=length(Yn)
        error('Y and Yn do not have the same dimension')
    end
else
    if size(Y,2)>1
        error('Y and Yn do not have the same dimension')
    end
end
if nargin<7
    svm=0;
end

m=median(M(:));
if m>1000
    M=M/100; % % reflectance
end

Mdl.Method='Regression';
Mdlvs.Method='Regression';
Mdl.Type='Complete';
Mdlvs.Type='Reduce';
% if size(M,3)==98||size(M,3)==144
Mdl.idx_wl1=15:size(M,3)-10;
Mdlvs.idx_wl1=15:size(M,3)-10;
M=M(:,:,15:end-10);
wl=wl(15:end-10);
% end

if mean(wl)>1000
    water = questdlg('Do you want to remove moisture bands ?','Moisture','Yes','No','No');
    if strcmp(water,'Yes')
        z1=[1094 1176];
        z2=[1339 1466];
        z3=[1768 2005];
        idxz1=wl>z1(1)&wl<z1(2);
        idxz2=wl>z2(1)&wl<z2(2);
        idxz3=wl>z3(1)&wl<z3(2);
        idxr=idxz1+idxz2+idxz3;
        idxr=idxr>0;% bandes enlevées
    else
        idxr=zeros(1,length(wl));
    end
else
    idxr=zeros(1,length(wl));
end
Mdl.idx_wl2=idxr;
Mdlvs.idx_wl2=idxr;

pseudoabs = questdlg('Do you want to convert your data in pseudo-absorbance ?','Pseudo-Absorbance','Yes','No','No');
datapretreat2 = questdlg('Do you want to preprocess the spectra ?','Preprocess','Yes','No','No');
if strcmp(datapretreat2,'Yes')
    pret={'Detrend','SNV','SNVD','MSC','D1','D2','SNV+D1','SNVD+D1','MSC+D1','SNV+D2','SNVD+D2','MSC+D2','CR'};
    datapretreat2b = listdlg('ListString',pret);
end
datapretreat = questdlg('Do you want to normalize your data?','Normalization','No','Center','Autoscaling','No');
dataref = questdlg('How your thickness are calculated?','Reference thickness','Min','Center','Max','Max');
if size(dy,2)==1
    datadref = inputdlg('What is the thickness:');
    datadref=repmat(str2double(datadref),1,length(d));
else
    datadref=dy(:,2);
    dy=dy(:,1);
end
SubSelect={'Random selection','Mean','Median','Multi-Mean','Multi-Median'};
[idx,~] = listdlg('PromptString','How would you like to subsample?','SelectionMode','single','ListString',SubSelect);
ech=SubSelect{idx};
Mdl.Proxy=Yn;
Mdlvs.Proxy=Yn;
if strcmp(pseudoabs,'Yes')
    Mdl.Spectral='Absorbance';
    Mdlvs.Spectral='Absorbance';
else
    Mdl.Spectral='Reflectance';
    Mdlvs.Spectral='Reflectance';
end
Mdl.Preprocessing=datapretreat;
Mdlvs.Preprocessing=datapretreat;

svch=questdlg('How do you want to select the variable?','Variable selection','Manual','Automatic','Manual');

%% Conversion in pseudo-absorbance
if strcmp(pseudoabs,'Yes')
    S=reshape(M,size(M,1)*size(M,2),size(M,3));
    [a,b]=find(S==0);
    S(a,b)=eps;
    S=log10(1./S);
    M=reshape(S,size(M,1),size(M,2),size(M,3));
end

%% Spectral preprocessing
if strcmp(datapretreat2,'Yes')
    M = reshape(AllPret(reshape(M,[],size(M,3)),wl,datapretreat2b),size(M,1),size(M,2),size(M,3));
    Mdl.SpectralProcessing=pret{datapretreat2b};
    Mdlvs.SpectralProcessing=pret{datapretreat2b};
else
    Mdl.SpectralProcessing='Raw';
    Mdlvs.SpectralProcessing='Raw';
end

%% Subsampling of the hyperspectral data at the reference resolution
dyb=zeros(2,length(dy));
if strcmp(dataref,'Max')
    for i=1:length(dy)
        dyb(1,i)=dy(i)-2/3*datadref(i);
        dyb(2,i)=dy(i)-1/3*datadref(i);
    end
else if strcmp(dataref,'Min')
        for i=1:length(dy)
            dyb(1,i)=dy(i)+1/3*datadref(i);
            dyb(2,i)=dy(i)+2/3*datadref(i);
        end
    else if strcmp(dataref,'Center')
            for i=1:length(dy)
                dyb(1,i)=dy(i)-1/6*datadref(i);
                dyb(2,i)=dy(i)+1/6*datadref(i);
            end
        end
    end
end

% Depth registration
iter=1;
diy=zeros(2,length(dy));
for i=1:length(dy)
    for j=1:2
        [~, a]=find(abs(d-dyb(j,i))==min(abs(d-dyb(j,i))));
        if abs(d(a(1))-dyb(j,i))<0.5
            dis(j,iter)=a(1);
            Ye(iter,:)=Y(i,:);
        else
            diy(j)=1;
        end
    end
    if abs(d(a(1))-dyb(j,i))<0.5
        iter=iter+1;
    end
end
if dis(2,end)==0
    dis=dis(:,1:end-1);
    Ye=Ye(1:end-1,:);
end

% Keep the central third
Mi=M;

% Preprocessing
Si=reshape(Mi,size(Mi,1)*size(Mi,2),size(Mi,3));
para1=mean(Si,1);para2=std(Si);

disx=dis(2,:)-dis(1,:);
if sum(disx==0)>0
    a=disx==0;
    disx=disx(1,a==0);
    Ye=Ye(a==0,:);
end
if sum(disx<0)>0
    [~,b]=find(disx<0);
    for i=1:length(b)
        disxi= dis(1,b(i));
        dis(1,b(i))=dis(2,b(i));
        dis(2,b(i))=disxi;
    end
end

% Subsampling methods
Sw=[];
if strcmp(ech,'Random selection')
    nbiter=100;
    for i = 1 : nbiter
        for j = 1:size(Ye,1)
            Mii=Mi(:,dis(1,j):dis(2,j),:);
            Si=reshape(Mii,[size(Mii,1)*size(Mii,2) size(Mii,3)]);
            % Remove outliers
            [~,outliers] = ACPauto(Si, 0, wl);
            Si=Si(outliers(:,7)==0,:);
            % Remove outliers
            [~,outliers] = ACPauto(Si, 0, wl);
            Si=Si(outliers(:,7)==0,:);
            
            rdm=rand(1,size(Si,1));
            [~,idx_rdm]=max(rdm);
            
            if strcmp(datapretreat,'Center')
                [Sw(i,j,:),~,~]=pretreat(Si(idx_rdm,:)','center',para1,para2);
                Mdl.PreprocessingPara.para1=para1;
                Mdl.PreprocessingPara.para2=para2;
                Mdlvs.PreprocessingPara.para1=para1;
                Mdlvs.PreprocessingPara.para2=para2;
            else if strcmp(datapretreat,'Autoscaling')
                    [Sw(i,j,:),~,~]=pretreat(Si(idx_rdm,:)','autoscaling',para1,para2);
                    Mdl.PreprocessingPara.para1=para1;
                    Mdl.PreprocessingPara.para2=para2;
                    Mdlvs.PreprocessingPara.para1=para1;
                    Mdlvs.PreprocessingPara.para2=para2;
                else
                    Sw(i,j,:)=Si(idx_rdm,:);
                end
            end
        end
    end
else if strcmp(ech,'Mean')
        for j = 1:size(Ye,1)
            Mii=Mi(:,dis(1,j):dis(2,j),:);
            Si=reshape(Mii,[size(Mii,1)*size(Mii,2) size(Mii,3)]);
            % Remove outliers
            [~,outliers] = ACPauto(Si, 0, wl);
            Si=Si(outliers(:,7)==0,:);
            % Remove outliers
            [~,outliers] = ACPauto(Si, 0, wl);
            Si=Si(outliers(:,7)==0,:);
            if strcmp(datapretreat,'Center')
                [Si,~,~]=pretreat(Si,'center',para1,para2);
                Mdl.PreprocessingPara.para1=para1;
                Mdl.PreprocessingPara.para2=para2;
                Mdlvs.PreprocessingPara.para1=para1;
                Mdlvs.PreprocessingPara.para2=para2;
            else if strcmp(datapretreat,'Autoscaling')
                    [Si,~,~]=pretreat(Si,'autoscaling',para1,para2);
                    Mdl.PreprocessingPara.para1=para1;
                    Mdl.PreprocessingPara.para2=para2;
                    Mdlvs.PreprocessingPara.para1=para1;
                    Mdlvs.PreprocessingPara.para2=para2;
                end
            end
            Sw(1,j,:)=mean(Si);
            Yee(j,:)=Ye(j,:);
        end
        nbiter=1;
    else if strcmp(ech,'Median')
            for j = 1:size(Ye,1)
                Mii=Mi(:,dis(1,j):dis(2,j),:);
                Si=reshape(Mii,[size(Mii,1)*size(Mii,2) size(Mii,3)]);
                % Remove outliers
                [~,outliers] = ACPauto(Si, 0, wl);
                Si=Si(outliers(:,7)==0,:);
                % Remove outliers
                [~,outliers] = ACPauto(Si, 0, wl);
                Si=Si(outliers(:,7)==0,:);
                if strcmp(datapretreat,'Center')
                    [Si,~,~]=pretreat(Si,'center',para1,para2);
                    Mdl.PreprocessingPara.para1=para1;
                    Mdl.PreprocessingPara.para2=para2;
                    Mdlvs.PreprocessingPara.para1=para1;
                    Mdlvs.PreprocessingPara.para2=para2;
                else if strcmp(datapretreat,'Autoscaling')
                        [Si,~,~]=pretreat(Si,'autoscaling',para1,para2);
                        Mdl.PreprocessingPara.para1=para1;
                        Mdl.PreprocessingPara.para2=para2;
                        Mdlvs.PreprocessingPara.para1=para1;
                        Mdlvs.PreprocessingPara.para2=para2;
                    end
                end
                Sw(1,j,:)=median(Si);
                Yee(j,:)=Ye(j,:);
            end
            nbiter=1;
        else if strcmp(ech,'Multi-Mean')
                Subset=floor(size(Mi,1)/50);
                if Subset==0
                    Subset=1;
                end
                Subsetpos=1:size(Mi,1)/Subset-1:size(Mi,1);
                Subsetpos(end)=size(Mi,1);
                for k=1:Subset
                    for j = 1:size(Ye,1)
                        Mii=Mi(Subsetpos(k):Subsetpos(k+1),dis(1,j):dis(2,j),:);
                        Si=reshape(Mii,[size(Mii,1)*size(Mii,2) size(Mii,3)]);
                        % Remove outliers
                        [~,outliers] = ACPauto(Si, 0, wl);
                        Si=Si(outliers(:,7)==0,:);
                        % Remove outliers
                        [~,outliers] = ACPauto(Si, 0, wl);
                        Si=Si(outliers(:,7)==0,:);
                        if strcmp(datapretreat,'Center')
                            [Si,~,~]=pretreat(Si,'center',para1,para2);
                            Mdl.PreprocessingPara.para1=para1;
                            Mdl.PreprocessingPara.para2=para2;
                            Mdlvs.PreprocessingPara.para1=para1;
                            Mdlvs.PreprocessingPara.para2=para2;
                        else if strcmp(datapretreat,'Autoscaling')
                                [Si,~,~]=pretreat(Si,'autoscaling',para1,para2);
                                Mdl.PreprocessingPara.para1=para1;
                                Mdl.PreprocessingPara.para2=para2;
                                Mdlvs.PreprocessingPara.para1=para1;
                                Mdlvs.PreprocessingPara.para2=para2;
                            end
                        end
                        Sw(1,j,:)=mean(Si);
                        Yee(j,:)=Ye(j,:);
                    end
                end
                nbiter=1;
            else if strcmp(ech,'Multi-Median')
                    Subset=floor(size(Mi,1)/50);
                    if Subset==0
                        Subset=1;
                    end
                    Subsetpos=1:size(Mi,1)/Subset-1:size(Mi,1);
                    Subsetpos(end)=size(Mi,1);
                    for k=1:Subset
                        for j = 1:size(Ye,1)
                            Mii=Mi(Subsetpos(k):Subsetpos(k+1),dis(1,j):dis(2,j),:);
                            Si=reshape(Mii,[size(Mii,1)*size(Mii,2) size(Mii,3)]);
                            % Remove outliers
                            [~,outliers] = ACPauto(Si, 0, wl);
                            Si=Si(outliers(:,7)==0,:);
                            % Remove outliers
                            [~,outliers] = ACPauto(Si, 0, wl);
                            Si=Si(outliers(:,7)==0,:);
                            if strcmp(datapretreat,'Center')
                                [Si,~,~]=pretreat(Si,'center',para1,para2);
                                Mdl.PreprocessingPara.para1=para1;
                                Mdl.PreprocessingPara.para2=para2;
                                Mdlvs.PreprocessingPara.para1=para1;
                                Mdlvs.PreprocessingPara.para2=para2;
                            else if strcmp(datapretreat,'Autoscaling')
                                    [Si,~,~]=pretreat(Si,'autoscaling',para1,para2);
                                    Mdl.PreprocessingPara.para1=para1;
                                    Mdl.PreprocessingPara.para2=para2;
                                    Mdlvs.PreprocessingPara.para1=para1;
                                    Mdlvs.PreprocessingPara.para2=para2;
                                end
                            end
                            Sw(1,j,:)=median(Si);
                            Yee(j,:)=Ye(j,:);
                        end
                    end
                    nbiter=1;
                end
            end
        end
    end
end
Swk=Sw;
Num=(1:size(Ye,1))';

% Saving subsampling data
DataPLS.S=Swk;
DataPLS.Y=Ye;
DataPLS.Num=Num;
DataPLS.Pret=datapretreat;
if strcmp(datapretreat,'center')||strcmp(datapretreat,'autoscaling')
    DataPLS.Pret1=para1;
    DataPLS.Pret2=para2;
else
    DataPLS.Pret1=[1 1];
    DataPLS.Pret2=[1 1];
end
assignin('base', 'DataPLS', DataPLS);

%% Histogram
% if size(Ye,2)<12
%     fig = figure;
%     fig.PaperPositionMode = 'auto';
%     fig.InvertHardcopy = 'off';
%     set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name', 'Histogram');
%     na=floor(sqrt(size(Ye,2)));
%     nb=ceil(sqrt(size(Ye,2)));
%     if na*nb<size(Ye,2)
%         na=nb;
%     end
%     for i=1:size(Ye,2)
%         subplot(na,nb,i);hist(Ye(:,i),100);grid on
%     end
% else
%     nby=1;
%     nbiter=floor(size(Ye,2)/12);
%     for i=1:nbiter
%         fig = figure;
%         fig.PaperPositionMode = 'auto';
%         fig.InvertHardcopy = 'off';
%         set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name', 'Histogram');
%         for j=1:12
%             subplot(3,4,j);hist(Ye(:,nby),100);grid on;title(Yn(nby))
%             nby=nby+1;
%         end
%     end
%     na=floor(sqrt(size(Ye,2)-nby+1));
%     nb=ceil(sqrt(size(Ye,2)-nby+1));
%     if na*nb<size(Ye,2)-nby+1
%         na=nb;
%     end
%     fig = figure;
%     fig.PaperPositionMode = 'auto';
%     fig.InvertHardcopy = 'off';
%     set(fig,'units','characters', 'position',[0 0 a_characters(1) a_characters(2)], 'name', 'Histogram');
%     for i=1:size(Ye,2)-nby+1
%         subplot(na,nb,i);hist(Ye(:,nby),100);grid on;title(Yn(nby))
%         nby=nby+1;
%     end
% end

clear M dm Mii Mi Mia Si S

R2cal_PLSR=zeros(size(Ye,2),nbiter,1);
R2val_PLSR=zeros(size(Ye,2),nbiter,1);
Nbvl_PLSR=zeros(nbiter,1);
h = waitbar(0,'PLS-CV');
iter=1;
for j=1:nbiter
    Spret=squeeze(Sw(j,:,:));
    
    % Subset function of Y distribution
    Ymax=max(Y);
    Ymin=min(Y);
    Yet=(Ymax-Ymin)/10;
    Ysep=[Ymin:Yet:Ymax];
    Ysep(1)=floor(Ysep(1));
    Ysep(end)=ceil(Ysep(end));
    for l=1:10
        [a,~]=find(Y>Ysep(l)&Y<Ysep(l+1));
        nb(l)=length(a);
    end
    ical=[];
    ival=[];
    itest=[];
    for l=1:10
        a=[];
        if nb(l)<0.1*max(nb)
            [a,~]=find(Y>Ysep(l)&Y<Ysep(l+1));
            ical=[ical; a];
        else if nb(l)<0.7*max(nb)
                [a,~]=find(Y>Ysep(l)&Y<Ysep(l+1));
                [trainInd,valInd,testInd] = dividerand(nb(l),70/100,15/100,15/100);
                ical=[ical; a(trainInd)];
                ival=[ival; a(valInd)];
                itest=[itest; a(testInd)];
            else
                [a,~]=find(Y>Ysep(l)&Y<Ysep(l+1));
                [trainInd,valInd,testInd] = dividerand(nb(l),50/100,25/100,25/100);
                ical=[ical; a(trainInd)];
                ival=[ival; a(valInd)];
                itest=[itest; a(testInd)];
            end
        end
    end
    trainInd=ical;
    valInd=ival;
    testInd=itest;
    ical=trainInd;
    ival=[valInd; testInd];
    
    for l=2:10
        %% ANN
        tic
        % Create a Fitting Network
        hiddenLayerSize = l;
        Mdl_ANNit{l,j} = fitnet(hiddenLayerSize);
        % Set up Division of Data for Training, Validation, Testing
        %         Mdl_ANN{1,j}.divideParam.trainRatio = 70/100;
        %         Mdl_ANN{1,j}.divideParam.valRatio = 15/100;
        %         Mdl_ANN{1,j}.divideParam.testRatio = 15/100;
        Mdl_ANNit{l,j}.divideFcn='divideind';
        Mdl_ANNit{l,j}.divideParam.trainInd = trainInd;
        Mdl_ANNit{l,j}.divideParam.valInd  = valInd;
        Mdl_ANNit{l,j}.divideParam.testInd = testInd;
        Mdl_ANNit{l,j}.trainParam.showWindow = false;
        % Train the Network
        [Mdl_ANNit{l,j},tr] = train(Mdl_ANNit{l,j},Spret(:,idxr==0)',Ye');
        Mdl_ANNtrit{l,j}=tr;
        % Test the Network
        for k=1:size(Y,2)
            Yp = Mdl_ANNit{l,j}(Spret(ical,idxr==0)');
            R2cal_ANNit(k,l,j)=corr2(Yp(k,:)',Ye(tr.trainInd,k))^2;
            SEC_ANNit(k,l,j)=sqrt(sum((Yp(k,:)'-Ye(ical,k)).^2)/length(ical));
            Yp = Mdl_ANNit{l,j}(Spret(ival,idxr==0)');
            R2val_ANNit(k,l,j)=corr2(Yp(k,:)',Ye([tr.valInd tr.testInd],k))^2;
            SEP_ANNit(k,l,j)=sqrt(sum((Yp(k,:)'-Ye(ival,k)).^2)/length(ival));
        end
        
        tann=toc;
    end
    R2m=mean([R2cal_ANNit(:,:,j); R2val_ANNit(:,:,j)]);
    R2p=(abs(SEC_ANNit(:,:,j)-SEP_ANNit(:,:,j))-min(abs(SEC_ANNit(:,:,j)-SEP_ANNit(:,:,j))))/(max(abs(SEC_ANNit(:,:,j)-SEP_ANNit(:,:,j)))-min(abs(SEC_ANNit(:,:,j)-SEP_ANNit(:,:,j))))+1;
    R2mp=R2m./R2p;
    [~,idxann]=max(R2mp);
    Mdl_ANN{1,j}=Mdl_ANNit{idxann,j};
    Mdl_ANNtr{1,j}=Mdl_ANNtrit{idxann,j};
    R2cal_ANN(k,j)=squeeze(R2cal_ANNit(k,idxann,j));
    SEC_ANN(k,j)=squeeze(SEC_ANNit(k,idxann,j));
    R2val_ANN(k,j)=squeeze(R2val_ANNit(k,idxann,j));
    SEP_ANN(k,j)=squeeze(SEP_ANNit(k,idxann,j));
    
    %% PLS
    tic
    Mdl_PLSR{1,j}=PLS(Spret(:,idxr==0),Ye,0, ical, ival);
    t=toc;
    Mdl_PLSR{1,j}.time=t;
    
    for k=1:size(Y,2)
        R2cal_PLSR(k,j)=Mdl_PLSR{1,j}.coeff.Recap{20,2}(k);
        R2val_PLSR(k,j)=Mdl_PLSR{1,j}.coeff.Recap{24,2}(k);
        SEC_PLSR(k,j)=Mdl_PLSR{1,j}.coeff.Recap{22,2}(k);
        SEP_PLSR(k,j)=Mdl_PLSR{1,j}.coeff.Recap{26,2}(k);
    end
    Nbvl_PLSR(j,1)=Mdl_PLSR{1,j}.coeff.Recap{19,2}(1);
    
    for l=2:10
        %% PLSBPNN
        tic
        % Create a Fitting Network
        hiddenLayerSize = l;
        Mdl_PLSBPNNit{l,j} = fitnet(hiddenLayerSize);
        % Set up Division of Data for Training, Validation, Testing
        %         Mdl_ANN{1,j}.divideParam.trainRatio = 70/100;
        %         Mdl_ANN{1,j}.divideParam.valRatio = 15/100;
        %         Mdl_ANN{1,j}.divideParam.testRatio = 15/100;
        Mdl_PLSBPNNit{l,j}.divideFcn='divideind';
        Mdl_PLSBPNNit{l,j}.divideParam.trainInd = trainInd;
        Mdl_PLSBPNNit{l,j}.divideParam.valInd  = valInd;
        Mdl_PLSBPNNit{l,j}.divideParam.testInd = testInd;
        Mdl_PLSBPNNit{l,j}.trainParam.showWindow = false;
        % Train the Network
        [Mdl_PLSBPNNit{l,j},tr] = train(Mdl_PLSBPNNit{l,j},(Spret(:,idxr==0)*Mdl_PLSR{1, j}.para.P)',Ye');
        Mdl_PLSBPNNtrit{l,j}=tr;
        % Test the Network
        for k=1:size(Y,2)
            Yp = Mdl_PLSBPNNit{l,j}((Spret(ical,idxr==0)*Mdl_PLSR{1, j}.para.P)');
            R2cal_PLSBPNNit(k,l,j)=corr2(Yp(k,:)',Ye(tr.trainInd,k))^2;
            SEC_PLSBPNNit(k,l,j)=sqrt(sum((Yp(k,:)'-Ye(ical,k)).^2)/length(ical));
            Yp = Mdl_PLSBPNNit{l,j}((Spret(ival,idxr==0)*Mdl_PLSR{1, j}.para.P)');
            R2val_PLSBPNNit(k,l,j)=corr2(Yp(k,:)',Ye([tr.valInd tr.testInd],k))^2;
            SEP_PLSBPNNit(k,l,j)=sqrt(sum((Yp(k,:)'-Ye(ival,k)).^2)/length(ival));
        end
        
        tann=toc;
    end
    R2m=mean([R2cal_PLSBPNNit(:,:,j); R2val_PLSBPNNit(:,:,j)]);
    R2p=(abs(SEC_PLSBPNNit(:,:,j)-SEP_PLSBPNNit(:,:,j))-min(abs(SEC_PLSBPNNit(:,:,j)-SEP_PLSBPNNit(:,:,j))))/(max(abs(SEC_PLSBPNNit(:,:,j)-SEP_PLSBPNNit(:,:,j)))-min(abs(SEC_PLSBPNNit(:,:,j)-SEP_PLSBPNNit(:,:,j))))+1;
    R2mp=R2m./R2p;
    [~,idxbpnn]=max(R2mp);
    Mdl_BPNN{1,j}=Mdl_PLSBPNNit{idxbpnn,j};
    Mdl_BPNNtr{1,j}=Mdl_PLSBPNNtrit{idxbpnn,j};
    Mdl_BPNNplsloadings{1,j}=Mdl_PLSR{1, j}.para.P;
    R2cal_BPNN(k,j)=squeeze(R2cal_PLSBPNNit(k,idxbpnn,j));
    SEC_BPNN(k,j)=squeeze(SEC_PLSBPNNit(k,idxbpnn,j));
    R2val_BPNN(k,j)=squeeze(R2val_PLSBPNNit(k,idxbpnn,j));
    SEP_BPNN(k,j)=squeeze(SEP_PLSBPNNit(k,idxbpnn,j));
    
    if size(Ye,2)==1
        %% MLR
        %             if size(Ye,2)==1
        %                 tic
        %                 [beta,Sigma,E,CovB,logL] = mvregress(Spret,Ye);
        %                 tmlr=toc;
        %             end
        
        %% Regression Tree
        tic
        Mdl_RT{1,j} = fitrtree(Spret(ical,idxr==0),Ye(ical,:));
        Yp = predict(Mdl_RT{1,j},Spret(ical,idxr==0));
        R2cal_RT(1,j)=corr(Yp,Ye(ical,:))^2;
        SEC_RT(1,j)=sqrt(sum((Yp-Ye(ical,:)).^2)/length(ical));
        Yp=[];
        Yp = predict(Mdl_RT{1,j},Spret(ival,idxr==0));
        R2val_RT(1,j)=corr(Yp,Ye(ival,:))^2;
        SEP_RT(1,j)=sqrt(sum((Yp-Ye(ival,:)).^2)/length(ival));
        Yp=[];
        ttree=toc;
        
        %% Regression Tree Ensembles
        tic
        Mdl_RTE{1,j} = fitrensemble(Spret(ical,idxr==0),Ye(ical,:),'Method','Bag');
        Yp = predict(Mdl_RTE{1,j},Spret(ical,idxr==0));
        R2cal_RTE(1,j)=corr(Yp,Ye(ical,:))^2;
        SEC_RTE(1,j)=sqrt(sum((Yp-Ye(ical,:)).^2)/length(ical));
        Yp=[];
        Yp = predict(Mdl_RTE{1,j},Spret(ival,idxr==0));
        R2val_RTE(1,j)=corr(Yp,Ye(ival,:))^2;
        SEP_RTE(1,j)=sqrt(sum((Yp-Ye(ival,:)).^2)/length(ival));
        Yp=[];
        trtree=toc;
        
        %% MARS / AresLab
        [Mdl_MARS{1,j}, Mdl_MARS{1,j}.time, Mdl_MARS{1,j}.resultsEval] = aresbuild(Spret(ical,idxr==0),Ye(ical,:));
        [Yp, BX] = arespredict(Mdl_MARS{1,j}, Spret(ical,idxr==0));
        R2cal_MARS(1,j)=corr(Yp,Ye(ical,:))^2;
        SEC_MARS(1,j)=sqrt(sum((Yp-Ye(ical,:)).^2)/length(ical));
        Yp=[];
        [Yp, BX] = arespredict(Mdl_MARS{1,j}, Spret(ival,idxr==0));
        R2val_MARS(1,j)=corr(Yp,Ye(ival,:))^2;
        SEP_MARS(1,j)=sqrt(sum((Yp-Ye(ival,:)).^2)/length(ival));
        Yp=[];
        
        %% SVM
        tic
        Mdl_SVM{1,j} = fitrsvm(Spret(ical,idxr==0),Ye(ical,:));
        Yp = predict(Mdl_SVM{1,j},Spret(ical,idxr==0));
        R2cal_SVM(1,j)=corr(Yp,Ye(ical,:))^2;
        SEC_SVM(1,j)=sqrt(sum((Yp-Ye(ical,:)).^2)/length(ical));
        Yp=[];
        Yp = predict(Mdl_SVM{1,j},Spret(ival,idxr==0));
        R2val_SVM(1,j)=corr(Yp,Ye(ival,:))^2;
        SEP_SVM(1,j)=sqrt(sum((Yp-Ye(ival,:)).^2)/length(ival));
        Yp=[];
        tsvm=toc;
    end
    
    h = waitbar(iter/(nbiter));
    iter=iter+1;
end
close(h)

Pret={'Raw' 'Detrend' 'MSC' 'SNV' 'SNVD' 'D1' 'D2' 'MSCD1' 'SNVD1' 'SNVDD1' 'MSCD2' 'SNVD2' 'SNVDD2' 'CR'};
if nbiter==1
    Mdl.ANN=Mdl_ANN;
    Mdl.ANNtr=Mdl_ANNtr;
    Mdl.ANN_Perf={'R2cal',median(R2cal_ANN),'SEC',median(SEC_ANN),'R2val',median(R2val_ANN),'SEP',median(SEP_ANN),'Nb neurons',idxann};
    Mdl.PLSR=Mdl_PLSR;
    Mdl.PLSR_Perf={'VL',Nbvl_PLSR,'R2cal',median(R2cal_PLSR),'SEC',median(SEC_PLSR),'R2val',median(R2val_PLSR),'SEP',median(SEP_PLSR)};
    Mdl.BPNN=Mdl_BPNN;
    Mdl.BPNNtr=Mdl_BPNNtr;
    Mdl.BPNNplsloadings=Mdl_PLSR{1,1}.para.P;
    Mdl.BPNN_Perf={'R2cal',median(R2cal_BPNN),'SEC',median(SEC_BPNN),'R2val',median(R2val_BPNN),'SEP',median(SEP_BPNN),'Nb neurons',idxbpnn};
    if size(Ye,2)==1
        Mdl.RT=Mdl_RT;
        Mdl.RT_Perf={'R2cal',R2cal_RT,'SEC',SEC_RT,'R2val',R2val_RT,'SEP',SEP_RT};
        Mdl.RTE=Mdl_RTE;
        Mdl.RTE_Perf={'R2cal',R2cal_RTE,'SEC',SEC_RTE,'R2val',R2val_RTE,'SEP',SEP_RTE};
        Mdl.MARS=Mdl_MARS;
        Mdl.MARS_Perf={'R2cal',R2cal_MARS,'SEC',SEC_MARS,'R2val',R2val_MARS,'SEP',SEP_MARS};
        Mdl.SVM=Mdl_SVM;
        Mdl.SVM_Perf={'R2cal',R2cal_SVM,'SEC',SEC_SVM,'R2val',R2val_SVM,'SEP',SEP_SVM};
    end
else % Optimum set
    R=R2cal_ANN+2*R2val_ANN;
    [~,p]=max(mean(squeeze(median(R,1)).*double(squeeze(median(R,1)>0))));
    tmp=squeeze(R2val_ANN(:,:,p));
    a=median(tmp(tmp>0));
    b=std(tmp(tmp>0));
    %     [~,c]=find(abs(a-tmp)==min(abs(a-tmp)));
    tmp=squeeze(R(:,:,p));
    [~,c]=max(tmp);
    Mdl.ANN=Mdl_ANN{p,c(1)};
    set_ANN=c(1);
    set_ANNr2=R(c(1));
    Mdl.ANN_Perf={'R2cal',R2cal_ANN(:,c(1),p),'SEC',SEC_ANN(:,c(1),p),'R2val',R2val_ANN(:,c(1),p),'SEP',SEP_ANN(:,c(1),p),'Random median',a,'Random std',b};
    
    R=R2cal_PLSR+2*R2val_PLSR;
    [~,p]=max(mean(squeeze(median(R,1)).*double(squeeze(median(R,1)>0))));
    tmp=squeeze(R2val_PLSR(:,:,p));
    a=median(tmp(tmp>0));
    b=std(tmp(tmp>0));
    %     [~,c]=find(abs(a-tmp)==min(abs(a-tmp)));
    tmp=squeeze(R(:,:,p));
    [~,c]=max(tmp);
    set_PLSR=c(1);
    set_PLSRr2=R(c(1));
    Mdl.PLSR=Mdl_PLSR{p,c(1)};
    Mdl.PLSR_Perf={'VL',Nbvl_PLSR(c(1),p),'R2cal',R2cal_PLSR(:,c(1),p),'SEC',SEC_PLSR(:,c(1),p),'R2val',R2val_PLSR(:,c(1),p),'SEP',SEP_PLSR(:,c(1),p),'Random median',a,'Random std',b};
    
    R=R2cal_BPNN+2*R2val_BPNN;
    [~,p]=max(mean(squeeze(median(R,1)).*double(squeeze(median(R,1)>0))));
    tmp=squeeze(R2val_BPNN(:,:,p));
    a=median(tmp(tmp>0));
    b=std(tmp(tmp>0));
    %     [~,c]=find(abs(a-tmp)==min(abs(a-tmp)));
    tmp=squeeze(R(:,:,p));
    [~,c]=max(tmp);
    Mdl.BPNN=Mdl_BPNN{p,c(1)};
    Mdl.BPNNplsloadings=Mdl_PLSR{p,c(1)}.para.P;
    set_BPNN=c(1);
    set_BPNNr2=R(c(1));
    Mdl.BPNN_Perf={'R2cal',R2cal_BPNN(:,c(1),p),'SEC',SEC_BPNN(:,c(1),p),'R2val',R2val_BPNN(:,c(1),p),'SEP',SEP_BPNN(:,c(1),p),'Random median',a,'Random std',b};
    
    if size(Ye,2)==1
        R=R2cal_RT+2*R2val_RT;
        [~,p]=max(mean(squeeze(median(R,1)).*double(squeeze(median(R,1)>0))));
        tmp=squeeze(R2val_RT(:,:,p));
        a=median(tmp(tmp>0));
        b=std(tmp(tmp>0));
        %         [~,c]=find(abs(a-tmp)==min(abs(a-tmp)));
        tmp=squeeze(R(:,:,p));
        [~,c]=max(tmp);
        set_RT=c(1);
        set_RTr2=R(c(1));
        Mdl.RT=Mdl_RT{p,c(1)};
        Mdl.RT_Perf={'R2cal',R2cal_RT(:,c(1),p),'SEC',SEC_RT(:,c(1),p),'R2val',R2val_RT(:,c(1),p),'SEP',SEP_RT(:,c(1),p),'Random median',a,'Random std',b};
        
        R=R2cal_RTE+2*R2val_RTE;
        [~,p]=max(mean(squeeze(median(R,1)).*double(squeeze(median(R,1)>0))));
        tmp=squeeze(R2val_RTE(:,:,p));
        a=median(tmp(tmp>0));
        b=std(tmp(tmp>0));
        %         [~,c]=find(abs(a-tmp)==min(abs(a-tmp)));
        tmp=squeeze(R(:,:,p));
        [~,c]=max(tmp);
        set_RTE=c(1);
        set_RTEr2=R(c(1));
        Mdl.RTE=Mdl_RTE{p,c(1)};
        Mdl.RTE_Perf={'R2cal',R2cal_RTE(:,c(1),p),'SEC',SEC_RTE(:,c(1),p),'R2val',R2val_RTE(:,c(1),p),'SEP',SEP_RTE(:,c(1),p),'Random median',a,'Random std',b};
        
        R=R2cal_MARS+2*R2val_MARS;
        [~,p]=max(mean(squeeze(median(R,1)).*double(squeeze(median(R,1)>0))));
        tmp=squeeze(R2val_MARS(:,:,p));
        a=median(tmp(tmp>0));
        b=std(tmp(tmp>0));
        %         [~,c]=find(abs(a-tmp)==min(abs(a-tmp)));
        tmp=squeeze(R(:,:,p));
        [~,c]=max(tmp);
        set_MARS=c(1);
        set_MARSr2=R(c(1));
        Mdl.MARS=Mdl_MARS{p,c(1)};
        Mdl.MARS_Perf={'R2cal',R2cal_MARS(:,c(1),p),'SEC',SEC_MARS(:,c(1),p),'R2val',R2val_MARS(:,c(1),p),'SEP',SEP_MARS(:,c(1),p),'Random median',a,'Random std',b};
        
        R=R2cal_SVM+2*R2val_SVM;
        [~,p]=max(mean(squeeze(median(R,1)).*double(squeeze(median(R,1)>0))));
        tmp=squeeze(R2val_SVM(:,:,p));
        a=median(tmp(tmp>0));
        b=std(tmp(tmp>0));
        %         [~,c]=find(abs(a-tmp)==min(abs(a-tmp)));
        tmp=squeeze(R(:,:,p));
        [~,c]=max(tmp);
        set_SVM=c(1);
        set_SVMr2=R(c(1));
        Mdl.SVM=Mdl_SVM{p,c(1)};
        Mdl.SVM_Perf={'R2cal',R2cal_SVM(:,c(1),p),'SEC',SEC_SVM(:,c(1),p),'R2val',R2val_SVM(:,c(1),p),'SEP',SEP_SVM(:,c(1),p),'Random median',a,'Random std',b};
    end
end

assignin('base', 'Modelwvs', Mdl);

%% Discrimination areas for each method are displayed
if nbiter==1
    if size(Ye,2)==1
        figure;
        plot(wl(idxr==0),sum(cell2mat(Mdl.ANN{1,1}.IW))/max(sum(cell2mat(Mdl.ANN{1,1}.IW))),'linewidth',2)
        hold on
        plot(wl(idxr==0),Mdl.PLSR{1,1}.para.B/max(Mdl.PLSR{1,1}.para.B),'linewidth',2)
        hold on
        plot(wl(idxr==0),Mdl.SVM{1,1}.Beta/max(Mdl.SVM{1,1}.Beta),'linewidth',2)
        legend('ANN','PLSR','SVM')
    else
        figure;
        [~,hLine1,hLine2] = plotyy(wl(idxr==0),sum(cell2mat(Mdl.ANN{1,1}.IW)),wl(idxr==0),sum(Mdl.PLSR{1,1}.para.B,2));
        set(hLine1,'LineWidth',2);
        set(hLine2,'LineWidth',2);
        legend('ANN','PLSR')
    end
else
    if size(Ye,2)==1
        figure;
        plot(wl(idxr==0),sum(cell2mat(Mdl.ANN.IW))/max(sum(cell2mat(Mdl.ANN.IW))),'linewidth',2)
        hold on
        plot(wl(idxr==0),Mdl.PLSR.para.B/max(Mdl.PLSR.para.B),'linewidth',2)
        hold on
        plot(wl(idxr==0),Mdl.SVM.Beta/max(Mdl.SVM.Beta),'linewidth',2)
        legend('ANN','PLSR','SVM')
    else
        figure;
        [~,hLine1,hLine2] = plotyy(wl(idxr==0),sum(cell2mat(Mdl.ANN.IW)),wl(idxr==0),sum(Mdl.PLSR.para.B,2));
        set(hLine1,'LineWidth',2);
        set(hLine2,'LineWidth',2);
        legend('ANN','PLSR')
    end
end
grid on
set(gca,'fontsize',14)

%% Variable selection
if nbiter>1
    if size(Ye,2)==1
        setr2=[set_ANNr2, set_PLSRr2, set_BPNNr2, set_RTr2, set_RTEr2, set_MARSr2, set_SVMr2];
        set_opt=[set_ANN, set_PLSR, set_BPNN, set_RT, set_RTE, set_MARS, set_SVM];
    else
        setr2=[set_ANNr2 set_PLSRr2 set_BPNNr2];
        set_opt=[set_ANN set_PLSR set_BPNN];
    end
    [~,idx] = max(setr2);
    idx_set=set_opt(idx);
    Spret=squeeze(Sw(idx_set,:,:));
else
    idx_set=1;
    Spret=squeeze(Sw(idx_set,:,:));
    set_PLSR=1;
end
SV=AllSelectVar(squeeze(Sw(idx_set,:,idxr==0)),Y,wl(idxr==0),Nbvl_PLSR(set_PLSR),1);

if svm==0
    list = {'Combined','VIP','UVE','CARS','RF','VCN'};
    [svm,~] = listdlg('PromptString','Select a variable selection method:',...
        'SelectionMode','single',...
        'ListString',list);
end

if svm==1
    [B,I] = sort(SV.tot,'descend');
    h=waitbar(0,'Test variable selection');
    for i=2:sum(B>median(B))
        Swk=Spret(:,idxr==0);
        waitbar(i/sum(B>median(B)))
        
        %% ANN
        tic
        % Create a Fitting Network
        hiddenLayerSize = idxann;
        Mdl_ANNvsi{i} = fitnet(hiddenLayerSize);
        % Set up Division of Data for Training, Validation, Testing
        Mdl_ANNvsi{i}.divideFcn='divideind';
        Mdl_ANNvsi{i}.divideParam.trainInd = trainInd;
        Mdl_ANNvsi{i}.divideParam.valInd  = valInd;
        Mdl_ANNvsi{i}.divideParam.testInd = testInd;
        Mdl_ANNvsi{i}.trainParam.showWindow = false;
        % Train the Network
        [Mdl_ANNvsi{i},tr] = train(Mdl_ANNvsi{i},Swk(:,I(1:i))',Ye');
        % Test the Network
        ical=tr.trainInd;
        ival=[tr.valInd tr.testInd];
        
        for k=1:size(Y,2)
            Yp = Mdl_ANNvsi{i}(Swk(ical,I(1:i))');
            Mdl_ANNvsiPerf.R2cal_ANN(k,i)=corr2(Yp(k,:)',Ye(tr.trainInd,k))^2;
            Mdl_ANNvsiPerf.SEC_ANN(k,i)=sqrt(sum((Yp(k,:)'-Ye(ical,k)).^2)/length(ical));
            Yp = Mdl_ANNvsi{i}(Swk(ival,I(1:i))');
            Mdl_ANNvsiPerf.R2val_ANN(k,i)=corr2(Yp(k,:)',Ye([tr.valInd tr.testInd],k))^2;
            Mdl_ANNvsiPerf.SEP_ANN(k,i)=sqrt(sum((Yp(k,:)'-Ye(ival,k)).^2)/length(ival));
        end
        tann=toc;
        
        %% PLS
        tic
        Mdl_PLSRvsi{i}=PLS(Swk(:,I(1:i)),Ye,0, ical, ival);
        t=toc;
        Mdl_PLSRvsi{i}.time=t;
        
        for k=1:size(Y,2)
            Mdl_PLSRvsiPerf.R2cal_PLSR(k,i)=Mdl_PLSRvsi{i}.coeff.Recap{20,2}(k);
            Mdl_PLSRvsiPerf.R2val_PLSR(k,i)=Mdl_PLSRvsi{i}.coeff.Recap{24,2}(k);
            Mdl_PLSRvsiPerf.SEC_PLSR(k,i)=Mdl_PLSRvsi{i}.coeff.Recap{22,2}(k);
            Mdl_PLSRvsiPerf.SEP_PLSR(k,i)=Mdl_PLSRvsi{i}.coeff.Recap{26,2}(k);
        end
        Nbvl_PLSR(i)=Mdl_PLSRvsi{i}.coeff.Recap{19,2}(1);
        
        %% BPNN
        tic
        % Create a Fitting Network
        hiddenLayerSize = idxbpnn;
        Mdl_BPNNvsi{i} = fitnet(hiddenLayerSize);
        % Set up Division of Data for Training, Validation, Testing
        Mdl_BPNNvsi{i}.divideFcn='divideind';
        Mdl_BPNNvsi{i}.divideParam.trainInd = trainInd;
        Mdl_BPNNvsi{i}.divideParam.valInd  = valInd;
        Mdl_BPNNvsi{i}.divideParam.testInd = testInd;
        Mdl_BPNNvsi{i}.trainParam.showWindow = false;
        % Train the Network
        [Mdl_BPNNvsi{i},tr] = train(Mdl_BPNNvsi{i},(Swk(:,I(1:i))*Mdl_PLSRvsi{i}.para.P)',Ye');
        % Test the Network
        ical=tr.trainInd;
        ival=[tr.valInd tr.testInd];
        
        for k=1:size(Y,2)
            Yp = Mdl_BPNNvsi{i}((Swk(ical,I(1:i))*Mdl_PLSRvsi{i}.para.P)');
            Mdl_BPNNvsiPerf.R2cal_BPNN(k,i)=corr2(Yp(k,:)',Ye(tr.trainInd,k))^2;
            Mdl_BPNNvsiPerf.SEC_BPNN(k,i)=sqrt(sum((Yp(k,:)'-Ye(ical,k)).^2)/length(ical));
            Yp = Mdl_BPNNvsi{i}((Swk(ival,I(1:i))*Mdl_PLSRvsi{i}.para.P)');
            Mdl_BPNNvsiPerf.R2val_BPNN(k,i)=corr2(Yp(k,:)',Ye([tr.valInd tr.testInd],k))^2;
            Mdl_BPNNvsiPerf.SEP_BPNN(k,i)=sqrt(sum((Yp(k,:)'-Ye(ival,k)).^2)/length(ival));
        end
        tbpnn=toc;
        
        if size(Ye,2)==1
            %% MLR
            %             if size(Ye,2)==1
            %                 tic
            %                 [beta,Sigma,E,CovB,logL] = mvregress(Spret,Ye);
            %                 tmlr=toc;
            %             end
            
            %% Regression Tree
            tic
            Mdl_RTvsi{i} = fitrtree(Swk(ical,I(1:i)),Ye(ical,:));
            Yp = predict(Mdl_RTvsi{i},Swk(ical,I(1:i)));
            Mdl_RTvsiPerf.R2cal_RT(1,i)=corr(Yp,Ye(ical,:))^2;
            Mdl_RTvsiPerf.SEC_RT(1,i)=sqrt(sum((Yp-Ye(ical,:)).^2)/length(ical));
            Yp=[];
            Yp = predict(Mdl_RTvsi{i},Swk(ival,I(1:i)));
            Mdl_RTvsiPerf.R2val_RT(1,i)=corr(Yp,Ye(ival,:))^2;
            Mdl_RTvsiPerf.SEP_RT(1,i)=sqrt(sum((Yp-Ye(ival,:)).^2)/length(ival));
            Yp=[];
            ttree=toc;
            
            %% Regression Tree Ensembles
            tic
            Mdl_RTEvsi{i} = fitrensemble(Swk(ical,I(1:i)),Ye(ical,:),'Method','Bag');
            Yp = predict(Mdl_RTEvsi{i},Swk(ical,I(1:i)));
            Mdl_RTEvsiPerf.R2cal_RTE(1,i)=corr(Yp,Ye(ical,:))^2;
            Mdl_RTEvsiPerf.SEC_RTE(1,i)=sqrt(sum((Yp-Ye(ical,:)).^2)/length(ical));
            Yp=[];
            Yp = predict(Mdl_RTEvsi{i},Swk(ival,I(1:i)));
            Mdl_RTEvsiPerf.R2val_RTE(1,i)=corr(Yp,Ye(ival,:))^2;
            Mdl_RTEvsiPerf.SEP_RTE(1,i)=sqrt(sum((Yp-Ye(ival,:)).^2)/length(ival));
            Yp=[];
            trtree=toc;
            
            %% MARS / AresLab
            [Mdl_MARSvsi{i}, Mdl_MARSvsi{i}.time, Mdl_MARSvsi{i}.resultsEval] = aresbuild(Swk(ical,I(1:i)),Ye(ical,:));
            [Yp, BX] = arespredict(Mdl_MARSvsi{i}, Swk(ical,I(1:i)));
            Mdl_MARSvsiPerf.R2cal_MARS(1,i)=corr(Yp,Ye(ical,:))^2;
            Mdl_MARSvsiPerf.SEC_MARS(1,i)=sqrt(sum((Yp-Ye(ical,:)).^2)/length(ical));
            Yp=[];
            [Yp, BX] = arespredict(Mdl_MARSvsi{i}, Swk(ival,I(1:i)));
            Mdl_MARSvsiPerf.R2val_MARS(1,i)=corr(Yp,Ye(ival,:))^2;
            Mdl_MARSvsiPerf.SEP_MARS(1,i)=sqrt(sum((Yp-Ye(ival,:)).^2)/length(ival));
            Yp=[];
            
            %% SVM
            tic
            Mdl_SVMvsi{i} = fitrsvm(Swk(ical,I(1:i)),Ye(ical,:));
            Yp = predict(Mdl_SVMvsi{i},Swk(ical,I(1:i)));
            Mdl_SVMvsiPerf.R2cal_SVM(1,i)=corr(Yp,Ye(ical,:))^2;
            Mdl_SVMvsiPerf.SEC_SVM(1,i)=sqrt(sum((Yp-Ye(ical,:)).^2)/length(ical));
            Yp=[];
            Yp = predict(Mdl_SVMvsi{i},Swk(ival,I(1:i)));
            Mdl_SVMvsiPerf.R2val_SVM(1,i)=corr(Yp,Ye(ival,:))^2;
            Mdl_SVMvsiPerf.SEP_SVM(1,i)=sqrt(sum((Yp-Ye(ival,:)).^2)/length(ival));
            Yp=[];
            tsvm=toc;
        end
    end
    close(h)
end

%% ANN
if strcmp(svch,'Automatic')
    R2m=mean([Mdl_ANNvsiPerf.R2cal_ANN; Mdl_ANNvsiPerf.R2val_ANN]);
    R2p=(abs(Mdl_ANNvsiPerf.SEC_ANN-Mdl_ANNvsiPerf.SEP_ANN)-min(abs(Mdl_ANNvsiPerf.SEC_ANN-Mdl_ANNvsiPerf.SEP_ANN)))/(max(abs(Mdl_ANNvsiPerf.SEC_ANN-Mdl_ANNvsiPerf.SEP_ANN))-min(abs(Mdl_ANNvsiPerf.SEC_ANN-Mdl_ANNvsiPerf.SEP_ANN)))+1;
    R2mp=R2m./R2p;
    [~,sv]=max(R2mp);
else
    figure
    subplot(211)
    plot(1:sum(B>median(B)),Mdl_ANNvsiPerf.R2cal_ANN,'*',1:sum(B>median(B)),Mdl_ANNvsiPerf.R2val_ANN,'*')
    grid on
    ylim([0 1])
    legend('Calibration','Validation')
    set(gca,'fontsize',14,'xtick',0:5:sum(B>median(B)))
    title('R²')
    subplot(212)
    yyaxis left
    plot(1:sum(B>median(B)),Mdl_ANNvsiPerf.SEC_ANN,'*')
    ylim([0 2*median([Mdl_ANNvsiPerf.SEC_ANN Mdl_ANNvsiPerf.SEP_ANN])])
    yyaxis right
    plot(1:sum(B>median(B)),Mdl_ANNvsiPerf.SEP_ANN,'*')
    grid on
    ylim([0 2*median([Mdl_ANNvsiPerf.SEC_ANN Mdl_ANNvsiPerf.SEP_ANN])])
    legend('Calibration','Validation')
    title('SE')
    set(gca,'fontsize',14,'xtick',0:5:sum(B>median(B)))
    
    sv=ginput(1);
end

Mdlvs.ANN=Mdl_ANNvsi{round(sv(1))};

wl_vs=wl(idxr==0);
Mdlvs.ANN_Perf.wl_vs=wl_vs(I(1:round(sv(1))));
Mdlvs.ANN_Perf.idx_vs=I(1:round(sv(1)));
Mdlvs.ANN_Perf.R2cal=Mdl_ANNvsiPerf.R2cal_ANN(round(sv(1)));
Mdlvs.ANN_Perf.SEC=Mdl_ANNvsiPerf.SEC_ANN(round(sv(1)));
Mdlvs.ANN_Perf.R2val=Mdl_ANNvsiPerf.R2val_ANN(round(sv(1)));
Mdlvs.ANN_Perf.SEP=Mdl_ANNvsiPerf.SEP_ANN(round(sv(1)));

%% PLSR
if strcmp(svch,'Automatic')
    R2m=mean([Mdl_PLSRvsiPerf.R2cal_PLSR; Mdl_PLSRvsiPerf.R2val_PLSR]);
    R2p=(abs(Mdl_PLSRvsiPerf.SEC_PLSR-Mdl_PLSRvsiPerf.SEP_PLSR)-min(abs(Mdl_PLSRvsiPerf.SEC_PLSR-Mdl_PLSRvsiPerf.SEP_PLSR)))/(max(abs(Mdl_PLSRvsiPerf.SEC_PLSR-Mdl_PLSRvsiPerf.SEP_PLSR))-min(abs(Mdl_PLSRvsiPerf.SEC_PLSR-Mdl_PLSRvsiPerf.SEP_PLSR)))+1;
    R2mp=R2m./R2p;
    [~,sv]=max(R2mp);
else
    figure
    subplot(211)
    plot(1:sum(B>median(B)),Mdl_PLSRvsiPerf.R2cal_PLSR,'*',1:sum(B>median(B)),Mdl_PLSRvsiPerf.R2val_PLSR,'*')
    grid on
    ylim([0 1])
    legend('Calibration','Validation')
    set(gca,'fontsize',14,'xtick',0:5:sum(B>median(B)))
    title('R²')
    subplot(212)
    yyaxis left
    plot(1:sum(B>median(B)),Mdl_PLSRvsiPerf.SEC_PLSR,'*')
    ylim([0 2*median([Mdl_PLSRvsiPerf.SEC_PLSR Mdl_PLSRvsiPerf.SEP_PLSR])])
    yyaxis right
    plot(1:sum(B>median(B)),Mdl_PLSRvsiPerf.SEP_PLSR,'*')
    grid on
    ylim([0 2*median([Mdl_PLSRvsiPerf.SEC_PLSR Mdl_PLSRvsiPerf.SEP_PLSR])])
    legend('Calibration','Validation')
    title('SE')
    set(gca,'fontsize',14,'xtick',0:5:sum(B>median(B)))
    
    sv=ginput(1);
end

Mdlvs.PLSR=Mdl_PLSRvsi{round(sv(1))};

wl_vs=wl(idxr==0);
Mdlvs.PLSR_Perf.wl_vs=wl_vs(I(1:round(sv(1))));
Mdlvs.PLSR_Perf.idx_vs=I(1:round(sv(1)));
Mdlvs.PLSR_Perf.R2cal=Mdl_PLSRvsiPerf.R2cal_PLSR(round(sv(1)));
Mdlvs.PLSR_Perf.SEC=Mdl_PLSRvsiPerf.SEC_PLSR(round(sv(1)));
Mdlvs.PLSR_Perf.R2val=Mdl_PLSRvsiPerf.R2val_PLSR(round(sv(1)));
Mdlvs.PLSR_Perf.SEP=Mdl_PLSRvsiPerf.SEP_PLSR(round(sv(1)));

%% BPNN
if strcmp(svch,'Automatic')
    R2m=mean([Mdl_BPNNvsiPerf.R2cal_BPNN; Mdl_BPNNvsiPerf.R2val_BPNN]);
    R2p=(abs(Mdl_BPNNvsiPerf.SEC_BPNN-Mdl_BPNNvsiPerf.SEP_BPNN)-min(abs(Mdl_BPNNvsiPerf.SEC_BPNN-Mdl_BPNNvsiPerf.SEP_BPNN)))/(max(abs(Mdl_BPNNvsiPerf.SEC_BPNN-Mdl_BPNNvsiPerf.SEP_BPNN))-min(abs(Mdl_BPNNvsiPerf.SEC_BPNN-Mdl_BPNNvsiPerf.SEP_BPNN)))+1;
    R2mp=R2m./R2p;
    [~,sv]=max(R2mp);
else
    figure
    subplot(211)
    plot(1:sum(B>median(B)),Mdl_BPNNvsiPerf.R2cal_BPNN,'*',1:sum(B>median(B)),Mdl_BPNNvsiPerf.R2val_BPNN,'*')
    grid on
    ylim([0 1])
    legend('Calibration','Validation')
    set(gca,'fontsize',14,'xtick',0:5:sum(B>median(B)))
    title('R²')
    subplot(212)
    yyaxis left
    plot(1:sum(B>median(B)),Mdl_BPNNvsiPerf.SEC_BPNN,'*')
    ylim([0 2*median([Mdl_BPNNvsiPerf.SEC_BPNN Mdl_BPNNvsiPerf.SEP_BPNN])])
    yyaxis right
    plot(1:sum(B>median(B)),Mdl_BPNNvsiPerf.SEP_BPNN,'*')
    grid on
    ylim([0 2*median([Mdl_BPNNvsiPerf.SEC_BPNN Mdl_BPNNvsiPerf.SEP_BPNN])])
    legend('Calibration','Validation')
    title('SE')
    set(gca,'fontsize',14,'xtick',0:5:sum(B>median(B)))
    
    sv=ginput(1);
end

Mdlvs.BPNN=Mdl_BPNNvsi{round(sv(1))};

wl_vs=wl(idxr==0);
Mdlvs.BPNN_Perf.PLSRloadings=Mdl_PLSRvsi{round(sv(1))}.para.P;
Mdlvs.BPNN_Perf.wl_vs=wl_vs(I(1:round(sv(1))));
Mdlvs.BPNN_Perf.idx_vs=I(1:round(sv(1)));
Mdlvs.BPNN_Perf.R2cal=Mdl_BPNNvsiPerf.R2cal_BPNN(round(sv(1)));
Mdlvs.BPNN_Perf.SEC=Mdl_BPNNvsiPerf.SEC_BPNN(round(sv(1)));
Mdlvs.BPNN_Perf.R2val=Mdl_BPNNvsiPerf.R2val_BPNN(round(sv(1)));
Mdlvs.BPNN_Perf.SEP=Mdl_BPNNvsiPerf.SEP_BPNN(round(sv(1)));

if size(Ye,2)==1
    %% Regression Tree
    if strcmp(svch,'Automatic')
        R2m=mean([Mdl_RTvsiPerf.R2cal_RT; Mdl_RTvsiPerf.R2val_RT]);
        R2p=(abs(Mdl_RTvsiPerf.SEC_RT-Mdl_RTvsiPerf.SEP_RT)-min(abs(Mdl_RTvsiPerf.SEC_RT-Mdl_RTvsiPerf.SEP_RT)))/(max(abs(Mdl_RTvsiPerf.SEC_RT-Mdl_RTvsiPerf.SEP_RT))-min(abs(Mdl_RTvsiPerf.SEC_RT-Mdl_RTvsiPerf.SEP_RT)))+1;
        R2mp=R2m./R2p;
        [~,sv]=max(R2mp);
    else
        figure
        subplot(211)
        plot(1:sum(B>median(B)),Mdl_RTvsiPerf.R2cal_RT,'*',1:sum(B>median(B)),Mdl_RTvsiPerf.R2val_RT,'*')
        grid on
        ylim([0 1])
        legend('Calibration','Validation')
        set(gca,'fontsize',14,'xtick',0:5:sum(B>median(B)))
        title('R²')
        subplot(212)
        yyaxis left
        plot(1:sum(B>median(B)),Mdl_RTvsiPerf.SEC_RT,'*')
        ylim([0 2*median([Mdl_RTvsiPerf.SEC_RT Mdl_RTvsiPerf.SEP_RT])])
        yyaxis right
        plot(1:sum(B>median(B)),Mdl_RTvsiPerf.SEP_RT,'*')
        grid on
        ylim([0 2*median([Mdl_RTvsiPerf.SEC_RT Mdl_RTvsiPerf.SEP_RT])])
        legend('Calibration','Validation')
        title('SE')
        set(gca,'fontsize',14,'xtick',0:5:sum(B>median(B)))
        
        sv=ginput(1);
    end
    
    Mdlvs.RT=Mdl_RTvsi{round(sv(1))};
    
    wl_vs=wl(idxr==0);
    Mdlvs.RT_Perf.wl_vs=wl_vs(I(1:round(sv(1))));
    Mdlvs.RT_Perf.idx_vs=I(1:round(sv(1)));
    Mdlvs.RT_Perf.R2cal=Mdl_RTvsiPerf.R2cal_RT(round(sv(1)));
    Mdlvs.RT_Perf.SEC=Mdl_RTvsiPerf.SEC_RT(round(sv(1)));
    Mdlvs.RT_Perf.R2val=Mdl_RTvsiPerf.R2val_RT(round(sv(1)));
    Mdlvs.RT_Perf.SEP=Mdl_RTvsiPerf.SEP_RT(round(sv(1)));
    
    %% Regression Tree Ensembles
    if strcmp(svch,'Automatic')
        R2m=mean([Mdl_RTEvsiPerf.R2cal_RTE; Mdl_RTEvsiPerf.R2val_RTE]);
        R2p=(abs(Mdl_RTEvsiPerf.SEC_RTE-Mdl_RTEvsiPerf.SEP_RTE)-min(abs(Mdl_RTEvsiPerf.SEC_RTE-Mdl_RTEvsiPerf.SEP_RTE)))/(max(abs(Mdl_RTEvsiPerf.SEC_RTE-Mdl_RTEvsiPerf.SEP_RTE))-min(abs(Mdl_RTEvsiPerf.SEC_RTE-Mdl_RTEvsiPerf.SEP_RTE)))+1;
        R2mp=R2m./R2p;
        [~,sv]=max(R2mp);
    else
        figure
        subplot(211)
        plot(1:sum(B>median(B)),Mdl_RTEvsiPerf.R2cal_RTE,'*',1:sum(B>median(B)),Mdl_RTEvsiPerf.R2val_RTE,'*')
        grid on
        ylim([0 1])
        legend('Calibration','Validation')
        set(gca,'fontsize',14,'xtick',0:5:sum(B>median(B)))
        title('R²')
        subplot(212)
        yyaxis left
        plot(1:sum(B>median(B)),Mdl_RTEvsiPerf.SEC_RTE,'*')
        ylim([0 2*median([Mdl_RTEvsiPerf.SEC_RTE Mdl_RTEvsiPerf.SEP_RTE])])
        yyaxis right
        plot(1:sum(B>median(B)),Mdl_RTEvsiPerf.SEP_RTE,'*')
        grid on
        ylim([0 2*median([Mdl_RTEvsiPerf.SEC_RTE Mdl_RTEvsiPerf.SEP_RTE])])
        legend('Calibration','Validation')
        title('SE')
        set(gca,'fontsize',14,'xtick',0:5:sum(B>median(B)))
        
        sv=ginput(1);
    end
    
    Mdlvs.RTE=Mdl_RTEvsi{round(sv(1))};
    
    wl_vs=wl(idxr==0);
    Mdlvs.RTE_Perf.wl_vs=wl_vs(I(1:round(sv(1))));
    Mdlvs.RTE_Perf.idx_vs=I(1:round(sv(1)));
    Mdlvs.RTE_Perf.R2cal=Mdl_RTEvsiPerf.R2cal_RTE(round(sv(1)));
    Mdlvs.RTE_Perf.SEC=Mdl_RTEvsiPerf.SEC_RTE(round(sv(1)));
    Mdlvs.RTE_Perf.R2val=Mdl_RTEvsiPerf.R2val_RTE(round(sv(1)));
    Mdlvs.RTE_Perf.SEP=Mdl_RTEvsiPerf.SEP_RTE(round(sv(1)));
    
    %% MARS / AresLab
    if strcmp(svch,'Automatic')
        R2m=mean([Mdl_MARSvsiPerf.R2cal_MARS; Mdl_MARSvsiPerf.R2val_MARS]);
        R2p=(abs(Mdl_MARSvsiPerf.SEC_MARS-Mdl_MARSvsiPerf.SEP_MARS)-min(abs(Mdl_MARSvsiPerf.SEC_MARS-Mdl_MARSvsiPerf.SEP_MARS)))/(max(abs(Mdl_MARSvsiPerf.SEC_MARS-Mdl_MARSvsiPerf.SEP_MARS))-min(abs(Mdl_MARSvsiPerf.SEC_MARS-Mdl_MARSvsiPerf.SEP_MARS)))+1;
        R2mp=R2m./R2p;
        [~,sv]=max(R2mp);
    else
        figure
        subplot(211)
        plot(1:sum(B>median(B)),Mdl_MARSvsiPerf.R2cal_MARS,'*',1:sum(B>median(B)),Mdl_MARSvsiPerf.R2val_MARS,'*')
        grid on
        ylim([0 1])
        legend('Calibration','Validation')
        set(gca,'fontsize',14,'xtick',0:5:sum(B>median(B)))
        title('R²')
        subplot(212)
        yyaxis left
        plot(1:sum(B>median(B)),Mdl_MARSvsiPerf.SEC_MARS,'*')
        ylim([0 2*median([Mdl_MARSvsiPerf.SEC_MARS Mdl_MARSvsiPerf.SEP_MARS])])
        yyaxis right
        plot(1:sum(B>median(B)),Mdl_MARSvsiPerf.SEP_MARS,'*')
        grid on
        ylim([0 2*median([Mdl_MARSvsiPerf.SEC_MARS Mdl_MARSvsiPerf.SEP_MARS])])
        legend('Calibration','Validation')
        title('SE')
        set(gca,'fontsize',14,'xtick',0:5:sum(B>median(B)))
        
        sv=ginput(1);
    end
    
    Mdlvs.MARS=Mdl_MARSvsi{round(sv(1))};
    
    wl_vs=wl(idxr==0);
    Mdlvs.MARS_Perf.wl_vs=wl_vs(I(1:round(sv(1))));
    Mdlvs.MARS_Perf.idx_vs=I(1:round(sv(1)));
    Mdlvs.MARS_Perf.R2cal=Mdl_MARSvsiPerf.R2cal_MARS(round(sv(1)));
    Mdlvs.MARS_Perf.SEC=Mdl_MARSvsiPerf.SEC_MARS(round(sv(1)));
    Mdlvs.MARS_Perf.R2val=Mdl_MARSvsiPerf.R2val_MARS(round(sv(1)));
    Mdlvs.MARS_Perf.SEP=Mdl_MARSvsiPerf.SEP_MARS(round(sv(1)));
    
    %% SVM % pb type Y
    if strcmp(svch,'Automatic')
        R2m=mean([Mdl_SVMvsiPerf.R2cal_SVM; Mdl_SVMvsiPerf.R2val_SVM]);
        R2p=(abs(Mdl_SVMvsiPerf.SEC_SVM-Mdl_SVMvsiPerf.SEP_SVM)-min(abs(Mdl_SVMvsiPerf.SEC_SVM-Mdl_SVMvsiPerf.SEP_SVM)))/(max(abs(Mdl_SVMvsiPerf.SEC_SVM-Mdl_SVMvsiPerf.SEP_SVM))-min(abs(Mdl_SVMvsiPerf.SEC_SVM-Mdl_SVMvsiPerf.SEP_SVM)))+1;
        R2mp=R2m./R2p;
        [~,sv]=max(R2mp);
    else
        figure
        subplot(211)
        plot(1:sum(B>median(B)),Mdl_SVMvsiPerf.R2cal_SVM,'*',1:sum(B>median(B)),Mdl_SVMvsiPerf.R2val_SVM,'*')
        grid on
        ylim([0 1])
        legend('Calibration','Validation')
        set(gca,'fontsize',14,'xtick',0:5:sum(B>median(B)))
        title('R²')
        subplot(212)
        yyaxis left
        plot(1:sum(B>median(B)),Mdl_SVMvsiPerf.SEC_SVM,'*')
        ylim([0 2*median([Mdl_SVMvsiPerf.SEC_SVM Mdl_SVMvsiPerf.SEP_SVM])])
        yyaxis right
        plot(1:sum(B>median(B)),Mdl_SVMvsiPerf.SEP_SVM,'*')
        grid on
        ylim([0 2*median([Mdl_SVMvsiPerf.SEC_SVM Mdl_SVMvsiPerf.SEP_SVM])])
        legend('Calibration','Validation')
        title('SE')
        set(gca,'fontsize',14,'xtick',0:5:sum(B>median(B)))
        
        sv=ginput(1);
    end
    
    Mdlvs.SVM=Mdl_SVMvsi{round(sv(1))};
    
    wl_vs=wl(idxr==0);
    Mdlvs.SVM_Perf.wl_vs=wl_vs(I(1:round(sv(1))));
    Mdlvs.SVM_Perf.idx_vs=I(1:round(sv(1)));
    Mdlvs.SVM_Perf.R2cal=Mdl_SVMvsiPerf.R2cal_SVM(round(sv(1)));
    Mdlvs.SVM_Perf.SEC=Mdl_SVMvsiPerf.SEC_SVM(round(sv(1)));
    Mdlvs.SVM_Perf.R2val=Mdl_SVMvsiPerf.R2val_SVM(round(sv(1)));
    Mdlvs.SVM_Perf.SEP=Mdl_SVMvsiPerf.SEP_SVM(round(sv(1)));
end

if size(Ye,2)==1
    Name={'Nb wl','R2cal','SEC','R2val','SEP'}';
    ANNraw=[length(wl_vs) Mdl.ANN_Perf{1, 2}  Mdl.ANN_Perf{1, 4}  Mdl.ANN_Perf{1, 6}  Mdl.ANN_Perf{1, 8}]';
    PLSRraw=[length(wl_vs) Mdl.PLSR_Perf{1, 4}  Mdl.PLSR_Perf{1, 6}  Mdl.PLSR_Perf{1, 8}  Mdl.PLSR_Perf{1, 10}]';
    BPNNraw=[length(wl_vs) Mdl.BPNN_Perf{1, 2}  Mdl.BPNN_Perf{1, 4}  Mdl.BPNN_Perf{1, 6}  Mdl.BPNN_Perf{1, 8}]';
    RTraw=[length(wl_vs) Mdl.RT_Perf{1, 2}  Mdl.RT_Perf{1, 4}  Mdl.RT_Perf{1, 6}  Mdl.RT_Perf{1, 8}]';
    RTEraw=[length(wl_vs) Mdl.RTE_Perf{1, 2}  Mdl.RTE_Perf{1, 4}  Mdl.RTE_Perf{1, 6}  Mdl.RTE_Perf{1, 8}]';
    MARSraw=[length(wl_vs) Mdl.MARS_Perf{1, 2}  Mdl.MARS_Perf{1, 4}  Mdl.MARS_Perf{1, 6}  Mdl.MARS_Perf{1, 8}]';
    SVMraw=[length(wl_vs) Mdl.SVM_Perf{1, 2}  Mdl.SVM_Perf{1, 4}  Mdl.SVM_Perf{1, 6}  Mdl.SVM_Perf{1, 8}]';
    
    ANNsv=[length(Mdlvs.ANN_Perf.wl_vs) Mdlvs.ANN_Perf.R2cal  Mdlvs.ANN_Perf.SEC  Mdlvs.ANN_Perf.R2val  Mdlvs.ANN_Perf.SEP]';
    PLSRsv=[length(Mdlvs.PLSR_Perf.wl_vs) Mdlvs.PLSR_Perf.R2cal  Mdlvs.PLSR_Perf.SEC  Mdlvs.PLSR_Perf.R2val  Mdlvs.PLSR_Perf.SEP]';
    BPNNsv=[length(Mdlvs.BPNN_Perf.wl_vs) Mdlvs.BPNN_Perf.R2cal  Mdlvs.BPNN_Perf.SEC  Mdlvs.BPNN_Perf.R2val  Mdlvs.BPNN_Perf.SEP]';
    RTsv=[length(Mdlvs.RT_Perf.wl_vs) Mdlvs.RT_Perf.R2cal  Mdlvs.RT_Perf.SEC  Mdlvs.RT_Perf.R2val  Mdlvs.RT_Perf.SEP]';
    RTEsv=[length(Mdlvs.RTE_Perf.wl_vs) Mdlvs.RTE_Perf.R2cal  Mdlvs.RTE_Perf.SEC  Mdlvs.RTE_Perf.R2val  Mdlvs.RTE_Perf.SEP]';
    MARSsv=[length(Mdlvs.MARS_Perf.wl_vs) Mdlvs.MARS_Perf.R2cal  Mdlvs.MARS_Perf.SEC  Mdlvs.MARS_Perf.R2val  Mdlvs.MARS_Perf.SEP]';
    SVMsv=[length(Mdlvs.SVM_Perf.wl_vs) Mdlvs.SVM_Perf.R2cal  Mdlvs.SVM_Perf.SEC  Mdlvs.SVM_Perf.R2val  Mdlvs.SVM_Perf.SEP]';
    
    Recap=table(Name,ANNraw,ANNsv,PLSRraw,PLSRsv,BPNNraw,BPNNsv,RTraw,RTsv,RTEraw,RTEsv,MARSraw,MARSsv,SVMraw,SVMsv)
    Mdlvs.recap=Recap;
else
    Name={'Nb wl','R2cal','SEC','R2val','SEP'}';
    ANNraw=[length(Mdlvs.ANN_Perf.wl) Mdl.ANN_Perf{1, 2}  Mdl.ANN_Perf{1, 4}  Mdl.ANN_Perf{1, 6}  Mdl.ANN_Perf{1, 8}]';
    PLSRraw=[length(Mdlvs.PLSR_Perf.wl) Mdl.PLSR_Perf{1, 4}  Mdl.PLSR_Perf{1, 6}  Mdl.PLSR_Perf{1, 8}  Mdl.PLSR_Perf{1, 10}]';
    BPNNraw=[length(wl_vs) Mdl.BPNN_Perf{1, 2}  Mdl.BPNN_Perf{1, 4}  Mdl.BPNN_Perf{1, 6}  Mdl.BPNN_Perf{1, 8}]';
    
    ANNsv=[length(Mdlvs.ANN_Perf.wl_vs) Mdlvs.ANN_Perf.R2cal  Mdlvs.ANN_Perf.SEC  Mdlvs.ANN_Perf.R2val  Mdlvs.ANN_Perf.SEP]';
    PLSRsv=[length(Mdlvs.PLSR_Perf.wl_vs) Mdlvs.PLSR_Perf.R2cal  Mdlvs.PLSR_Perf.SEC  Mdlvs.PLSR_Perf.R2val  Mdlvs.PLSR_Perf.SEP]';
    BPNNsv=[length(Mdlvs.BPNN_Perf.wl_vs) Mdlvs.BPNN_Perf.R2cal  Mdlvs.BPNN_Perf.SEC  Mdlvs.BPNN_Perf.R2val  Mdlvs.BPNN_Perf.SEP]';
    
    Recap=table(Name,ANNraw,ANNsv,PLSRraw,PLSRsv,BPNNraw,BPNNsv)
    Mdlvs.recap=Recap;
end

assignin('base', 'Modelvs', Mdlvs);

end