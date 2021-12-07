function RegressionUnivariate_Create(M,d,wl,Y,dy,Yn)
% Function to estimate simple indices based on 1 or 2 wavelength(s) that
% are related to reference values (Y).
% INPUT :
%           M: Hyperspectral datacube (n*m*p)
%           d: Associated depth (1*m)
%           wl: Associated wavelengths (1*p)
%           Y: Reference values
%           dy: Reference depth
%           Yn: Name of the variable

%% Input verification
if size(M,2)~=length(d)
    error('M and d do not have the same size')
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
if size(Y,2)>1
    error('Y must be a column vector')
end

m=median(M(:));
if m>1000
    M=M/100;
end

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
ech = questdlg('How would you like to subsample?','Subsampling','Random selection','Mean','Median','Random selection');
Mdl.Proxy=Yn;
if strcmp(pseudoabs,'Yes')
    Mdl.Spectral='Absorbance';
else
    Mdl.Spectral='Reflectance';
end
Mdl.Preprocessing=datapretreat;

%% Conversion in pseudo-absorbance
if strcmp(pseudoabs,'Yes')
    S=reshape(M,size(M,1)*size(M,2),size(M,3));
    [a,b]=find(S==0);
    S(a,b)=eps;
    S=log(1./S);
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
Mi=M(round(1/3*size(M,1):2/3*size(M,1)),:,:);

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
            else if strcmp(datapretreat,'Autoscaling')
                    [Sw(i,j,:),~,~]=pretreat(Si(idx_rdm,:)','autoscaling',para1,para2);
                    Mdl.PreprocessingPara.para1=para1;
                    Mdl.PreprocessingPara.para2=para2;
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
            else if strcmp(datapretreat,'Autoscaling')
                    [Si,~,~]=pretreat(Si,'autoscaling',para1,para2);
                    Mdl.PreprocessingPara.para1=para1;
                    Mdl.PreprocessingPara.para2=para2;
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
                else if strcmp(datapretreat,'Autoscaling')
                        [Si,~,~]=pretreat(Si,'autoscaling',para1,para2);
                        Mdl.PreprocessingPara.para1=para1;
                        Mdl.PreprocessingPara.para2=para2;
                    end
                end
                Sw(1,j,:)=median(Si);
                Yee(j,:)=Ye(j,:);
            end
            nbiter=1;
        end
    end
end

clear M dm Mii Mi Mia Si S

h = waitbar(0,'Indice estimations');
iter=1;
Cc=zeros(1,length(wl));
RIc=zeros(1,length(wl)^2-length(wl));
RIc_exp=zeros(1,length(wl)^2-length(wl));
RIwl=zeros(length(wl)^2-length(wl),2);
NDIn=zeros(1,length(wl)^2-length(wl));
NDIn_exp=zeros(1,length(wl)^2-length(wl));
NDInwl=zeros(length(wl)^2-length(wl),2);
NDIp=zeros(1,(length(wl)^2-length(wl))/2);
NDIp_exp=zeros(1,(length(wl)^2-length(wl))/2);
NDIpwl=zeros((length(wl)^2-length(wl))/2,2);
DIc=zeros(1,(length(wl)^2-length(wl))/2);
DIc_exp=zeros(1,(length(wl)^2-length(wl))/2);
DIwl=zeros((length(wl)^2-length(wl))/2,2);
for j=1:nbiter
    Spret=squeeze(Sw(j,:,:));
    
    % Univariate
    for k=1:size(Spret,2)
        Cc(k)=corr(Spret(:,k),Ye);
    end
    
    iter_ri=1;
    iter_di=1;
    for k=1:size(Spret,2)
        for l=1:size(Spret,2)
            if k~=l
                % Ratio Indices
                tmp=Spret(:,k)./Spret(:,l);
                RIc(iter_ri)=corr(tmp,Ye);
                RIc_exp(iter_ri)=corr(real(log(abs(tmp))),Ye);
                RIwl(iter_ri,:)=[wl(k) wl(l)];
		IM_RI(k,l)=RIc_exp(iter_ri);
                
                if l>k
                    % Normalized Difference Indices
                    tmp=(Spret(:,k)-Spret(:,l))./(Spret(:,k)+Spret(:,l));
                    NDIn(iter_di)=corr(tmp,Ye);
                    NDIn_exp(iter_di)=corr(real(log(abs(tmp))),Ye);
                    NDInwl(iter_di,:)=[wl(k) wl(l)];
			IM_NDIn(k,l)=NDIn_exp(iter_di);

                    tmp=(Spret(:,k)+Spret(:,l))./(Spret(:,k)-Spret(:,l));
                    NDIp(iter_di)=corr(tmp,Ye);
                    NDIp_exp(iter_di)=corr(real(log(abs(tmp))),Ye);
                    NDIpwl(iter_di,:)=[wl(k) wl(l)];
			IM_NDIp(k,l)=NDIp_exp(iter_di);
                    
                    % Difference Indices
                    tmp=Spret(:,k)-Spret(:,l);
                    DIc(iter_di)=corr(tmp,Ye);
                    DIc_exp(iter_di)=corr(real(log(abs(tmp))),Ye);
                    DIwl(iter_di,:)=[wl(k) wl(l)];
			IM_DI(k,l)=DIc_exp(iter_di);

                    iter_di=iter_di+1;
                end
                iter_ri=iter_ri+1;
            end
        end
    end
    
    NDIc=[NDIn NDIp];
    NDIc_exp=[NDIn_exp NDIp_exp];
    NDIwl=[NDInwl; NDIpwl];
    
    h = waitbar(iter/(nbiter));
    iter=iter+1;
end
close(h)

% Display
figure;
subplot(221)
imagesc(wl,wl,IM_RI)
colormap(jet)
colorbar
caxis([-1 1])
xlabel('Wavelength (nm)')
ylabel('Wavelength (nm)')
title('Ratio indices')

subplot(222)
imagesc(wl,wl,IM_NDIn)
colormap(jet)
colorbar
caxis([-1 1])
xlabel('Wavelength (nm)')
ylabel('Wavelength (nm)')
title('Normalized Difference Indices')

subplot(223)
imagesc(wl,wl,IM_NDIp)
colormap(jet)
colorbar
caxis([-1 1])
xlabel('Wavelength (nm)')
ylabel('Wavelength (nm)')
title('Normalized Difference Indices')

subplot(224)
imagesc(wl,wl,IM_DI)
colormap(jet)
colorbar
caxis([-1 1])
xlabel('Wavelength (nm)')
ylabel('Wavelength (nm)')
title('Difference Indices')

figure
plot(wl,Cc,'linewidth',2)
hold on
plot([wl(1) wl(end)],[0.9 0.9],'g--',...
    [wl(1) wl(end)],[-0.9 -0.9],'g--',...
    [wl(1) wl(end)],[0.8 0.8],'r--',...
    [wl(1) wl(end)],[-0.8 -0.8],'r--')
xlabel('Wavelength (nm)')
ylabel('Correlation')
title('Univariate correlation')
grid on
ylim([-1 1])
set(gca,'fontsize',14)

figRI=figure;
subplot(211)
plot(RIc,'linewidth',2)
hold on
plot([0 length(RIc)],[0.9 0.9],'g--',...
    [0 length(RIc)],[-0.9 -0.9],'g--',...
    [0 length(RIc)],[0.8 0.8],'r--',...
    [0 length(RIc)],[-0.8 -0.8],'r--')
ylabel('Correlation')
title('Ratio indices linear correlation')
grid on
ylim([-1 1])
xlim([0 length(RIc)])
set(gca,'fontsize',14)

subplot(212)
plot(RIc_exp,'linewidth',2)
hold on
plot([0 length(RIc_exp)],[0.9 0.9],'g--',...
    [0 length(RIc_exp)],[-0.9 -0.9],'g--',...
    [0 length(RIc_exp)],[0.8 0.8],'r--',...
    [0 length(RIc_exp)],[-0.8 -0.8],'r--')
ylabel('Correlation')
title('Ratio indices exponential correlation')
grid on
ylim([-1 1])
xlim([0 length(RIc_exp)])
set(gca,'fontsize',14)
dcm_objRI = datacursormode(figRI);
set(dcm_objRI,'UpdateFcn',{@myupdatefcn,RIwl})

figNDI=figure;
subplot(211)
plot(NDIc,'linewidth',2)
hold on
plot([0 length(NDIc)],[0.9 0.9],'g--',...
    [0 length(NDIc)],[-0.9 -0.9],'g--',...
    [0 length(NDIc)],[0.8 0.8],'r--',...
    [0 length(NDIc)],[-0.8 -0.8],'r--')
ylabel('Correlation')
title('Normalized difference indices linear correlation')
grid on
ylim([-1 1])
xlim([0 length(NDIc)])
set(gca,'fontsize',14)

subplot(212)
plot(NDIc_exp,'linewidth',2)
hold on
plot([0 length(NDIc_exp)],[0.9 0.9],'g--',...
    [0 length(NDIc_exp)],[-0.9 -0.9],'g--',...
    [0 length(NDIc_exp)],[0.8 0.8],'r--',...
    [0 length(NDIc_exp)],[-0.8 -0.8],'r--')
ylabel('Correlation')
title('Normalized difference indices exponential correlation')
grid on
ylim([-1 1])
xlim([0 length(NDIc_exp)])
set(gca,'fontsize',14)
dcm_objNDI = datacursormode(figNDI);
set(dcm_objNDI,'UpdateFcn',{@myupdatefcn,NDIwl})

figDI=figure;
subplot(211)
plot(DIc,'linewidth',2)
hold on
plot([0 length(DIc)],[0.9 0.9],'g--',...
    [0 length(DIc)],[-0.9 -0.9],'g--',...
    [0 length(DIc)],[0.8 0.8],'r--',...
    [0 length(DIc)],[-0.8 -0.8],'r--')
ylabel('Correlation')
title('Difference indices linear correlation')
grid on
ylim([-1 1])
xlim([0 length(DIc)])
set(gca,'fontsize',14)

subplot(212)
plot(DIc_exp,'linewidth',2)
hold on
plot([0 length(DIc_exp)],[0.9 0.9],'g--',...
    [0 length(DIc_exp)],[-0.9 -0.9],'g--',...
    [0 length(DIc_exp)],[0.8 0.8],'r--',...
    [0 length(DIc_exp)],[-0.8 -0.8],'r--')
ylabel('Correlation')
title('Difference indices exponential correlation')
grid on
ylim([-1 1])
xlim([0 length(DIc_exp)])
set(gca,'fontsize',14)
dcm_objDI = datacursormode(figDI);
set(dcm_objDI,'UpdateFcn',{@myupdatefcn,NDIwl})

% Save
UniWl.Corr=Cc;
UniWl.Wl=wl;

RI.Corr=RIc;
RI.Wl=RIwl;

DI.Corr=DIc;
DI.Wl=DIwl;

NDI.Corr=NDIc;
NDI.Wl=NDIwl;

assignin('base', 'IndicesEstimationParamaters', Mdl);
assignin('base', 'RI', RI);
assignin('base', 'DI', DI);
assignin('base', 'NDI', NDI);
assignin('base', 'UnivWl', UniWl);

end

function txt = myupdatefcn(~,event_obj,t)
% Customizes text of data tips
pos = get(event_obj,'Position');
txt = {['X: ',num2str(pos(1))],...
    ['Y: ',num2str(pos(2))],...
    ['Wl1: ',num2str(round(t(pos(1),1)))],...
    ['Wl2: ',num2str(round(t(pos(1),2)))]};
end