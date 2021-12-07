function [Cartlab,lbl] = LabelledMapCreation(IM)
% Function to create a labelled map to estimate a discrimination model.
% INPUT:
%           IM : Hyperspectral datacube or RGB image of the sample
% OUTPUT:
%           Cartlab : Labelled map
%           Lbl : Label of the classes

% This is the Matalb toolbox from the papers:

% Jacq, K., Rapuc, W., Benoit, A., Develle, A.-L., Coquin, D.,
% Fanget, B., Perrette, Y.,  Sabatier, P., Wilhelm, B.,
% Debret, M., Arnaud, F., 2019. Sedimentary structures
% discriminations with hyperspectral imaging on
% sediment cores. Computers & Geosciences

% Please cite our papers if you use our code for your research.

% RGB creation if IM is an hyperspectral image
if size(IM,3)>3
    % ROI
    ROI=IM;
    if size(ROI,3)==144
        RGB=ROI(:,:,[107 111 124])/10000;
    else if size(ROI,3)==98
            RGB=ROI(:,:,[50 26 13])/10000;
        end
    end
else
    ROI=IM;
    RGB=ROI;
end
m=mean(RGB(:));
% RGB=RGB*ceil(0.5/m);
nbpixel=size(ROI,1)*size(ROI,2);

if size(ROI,3)==98||size(ROI,3)==144
    ROI=ROI(:,:,15:end-10);
end

% Initialisation
ROIclass=[];class=[];Mclass=[];
Cartlab=zeros(size(ROI,1),size(ROI,2));
disp('Selection des groupes pour classification supervis�e')
roisearch='Yes';
list = {'1','2','3','4','5'};
listi = {'1','2','3','4','5'};

answer = questdlg('Which ROI selection do you want to use?', ...
    'ROI selection', ...
    'Rectangle','Polygonal','Rectangle');

figure;
imagesc(RGB)

while strcmp(roisearch,'Yes')
    if strcmp(answer,'Rectangle')
        hBox = imrect;
        roiPositioni = wait(hBox);
        roiPositioni=round(roiPositioni);
    else
        BW = roipoly;
    end
    
    [indx,~] = listdlg('ListString',list);
    if strcmp(list{indx},listi{indx})==1
        lab=inputdlg({strcat('Class name ',num2str(indx))},'Labelling');
        lbl{indx}=lab{1,1};
    end
    
    if strcmp(answer,'Rectangle')
        Cartlab(roiPositioni(2):roiPositioni(2)+roiPositioni(4),roiPositioni(1):roiPositioni(1)+roiPositioni(3))=indx;
    else
        Cartlab=Cartlab+double(BW)*indx;
    end
    nbclassi(indx)=sum(Cartlab(:)==indx);
    
    list{indx}=strcat(lbl{indx},' : ',num2str(nbclassi(indx)/nbpixel),'%');
    
    roisearch=questdlg('Do you want to continue ?','ROI selection','Yes','No','Yes');
    
    assignin('base','Cartlab',Cartlab)
end

end