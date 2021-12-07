function [S,Y] = Subsampling3_2D(M,d,Y,dy,pas,pas2,typeEchM)
% Function to extract the spectra corresponding to the sampling areas.
% INPUT :
%           M: Hyperspectral datacube (n*m*p)
%           d: Associated depth (1*m)
%           Y: Reference values
%           dy: Reference depth
%           pas: Width of the reference sampling
%           pas2: Depth limits of the reference sampling (Max, Min, Center)
%           typeEchM: Type of spectral subsampling ('Random selection',
%               'Mean','Median','Multi-Mean','Multi-Median'
% OUTPUT:
%           S: Spectra extracted
%           Y: Associated reference values

if length(pas)==1
    pas=ones(length(Y),1)*pas;
end
wl=1:size(M,3);

if nargin<7
    SubSelect={'Random selection','Mean','Median','Multi-Mean','Multi-Median'};
    [idx,~] = listdlg('PromptString','How would you like to subsample?','SelectionMode','single','ListString',SubSelect);
    typeEchM=SubSelect{idx};
end

%% Subsampling of the hyperspectral data at the reference resolution
dyb=zeros(2,length(dy));
if strcmp(pas2,'Max')
    for i=1:length(dy)
        dyb(1,i)=dy(i)-2/3*pas(i);
        dyb(2,i)=dy(i)-1/3*pas(i);
    end
else if strcmp(pas2,'Min')
        for i=1:length(dy)
            dyb(1,i)=dy(i)+1/3*pas(i);
            dyb(2,i)=dy(i)+2/3*pas(i);
        end
    else if strcmp(pas2,'Center')
            for i=1:length(dy)
                dyb(1,i)=dy(i)-1/6*pas(i);
                dyb(2,i)=dy(i)+1/6*pas(i);
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
if strcmp(typeEchM,'Random selection')
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
            
            Sw(i,j,:)=Si(idx_rdm,:);
            
        end
    end
else if strcmp(typeEchM,'Mean')
        for j = 1:size(Ye,1)
            Mii=Mi(:,dis(1,j):dis(2,j),:);
            Si=reshape(Mii,[size(Mii,1)*size(Mii,2) size(Mii,3)]);
            % Remove outliers
            [~,outliers] = ACPauto(Si, 0, wl);
            Si=Si(outliers(:,7)==0,:);
            % Remove outliers
            [~,outliers] = ACPauto(Si, 0, wl);
            Si=Si(outliers(:,7)==0,:);
            Sw(1,j,:)=mean(Si);
            Yee(j,:)=Ye(j,:);
        end
        nbiter=1;
    else if strcmp(typeEchM,'Median')
            for j = 1:size(Ye,1)
                Mii=Mi(:,dis(1,j):dis(2,j),:);
                Si=reshape(Mii,[size(Mii,1)*size(Mii,2) size(Mii,3)]);
                % Remove outliers
                [~,outliers] = ACPauto(Si, 0, wl);
                Si=Si(outliers(:,7)==0,:);
                % Remove outliers
                [~,outliers] = ACPauto(Si, 0, wl);
                Si=Si(outliers(:,7)==0,:);
                Sw(1,j,:)=median(Si);
                Yee(j,:)=Ye(j,:);
            end
            nbiter=1;
        else if strcmp(typeEchM,'Multi-Mean')
                Subset=floor(size(Mi,1)/50);
                if Subset==0
                    Subset=1;
                end
                Subsetpos=1:size(Mi,1)/Subset-1:size(Mi,1);
                Subsetpos(end)=size(Mi,1);
                iter=1;
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
                        Sw(1,iter,:)=mean(Si);
                        Yee(iter,:)=Ye(j,:);
                        iter=iter+1;
                    end
                end
                nbiter=1;
            else if strcmp(typeEchM,'Multi-Median')
                    Subset=floor(size(Mi,1)/50);
                    if Subset==0
                        Subset=1;
                    end
                    Subsetpos=1:size(Mi,1)/Subset-1:size(Mi,1);
                    Subsetpos(end)=size(Mi,1);
                    iter=1;
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
                            Sw(1,iter,:)=median(Si);
                            Yee(iter,:)=Ye(j,:);
                            iter=iter+1;
                        end
                    end
                    nbiter=1;
                end
            end
        end
    end
end
S=Sw;
Y=Yee;
end