function [C,label] = IndiceRSCalculation(M,RGB,d,wl,indx)
% Function to calculate remote sensing indices
% INPUT :
%           M : Hyperspectral datacube
%           RGB : RGB image of the sample
%           d: Depth vector
%           wl : Wavelength
%           optional:
%               indx : Index of the index to estimate (0 to choose after)
% OUTPUT :
%           C : carte des indices
%           label : indices calculés

if mean(RGB(:))>1000
    RGB=RGB/10000;
end
RGB=RGB*(0.5/mean(RGB(:)));

map_red=cat(2,[0.5:0.001:1;zeros(1,501);zeros(1,501)],[ones(1,1001);0:0.001:1;0:0.001:1])';
map_green=cat(2,[zeros(1,501);0.5:0.001:1;zeros(1,501)],[0:0.001:1;ones(1,1001);0:0.001:1])';
map_orange=[ones(1,1001);1:-0.0005:.5;zeros(1,1001)]';
map_blue=cat(2,[zeros(1,501);zeros(1,501);0.5:0.001:1],[0:0.001:1;0:0.001:1;ones(1,1001)])';
map_gray=[0:0.001:1;0:0.001:1;0:0.001:1]';
p_r=[1 0 0];
p_g=[0 1 0];
p_b=[0 0 1];
p_o=[1 1 0];
p_gr=[0.5 0.5 0.5];

M=M/10000*100; % % Reflectance

%% VNIR
% Check the wavelengths range
iter=0;
if sum(wl<800&wl>600)>0
    
    index_vnir={'Difference Vegetation Index (DVI)',...
        'Enhanced Vegetation Index (EVI)',...
        'Global Environmental Monitoring Index (GEMI)',...
        'Green Atmospherically Resistant Index (GARI)',...
        'Green Chlorophyll Index (GCI)',...
        'Green Difference Vegetation Index (GDVI)',...
        'Green Leaf Index (GLI)',...
        'Green Normalized Difference Vegetation Index (GNDVI)',...
        'Green Optimized Soil Adjusted Vegetation Index (GOSAVI)',...
        'Green Ratio Vegetation Index (GRVI)',...
        'Green Soil Adjusted Vegetation Index (GSAVI)',...
        'Infrared Percentage Vegetation Index (IPVI)',...
        'Leaf Area Index (LAI)',...
        'Modified Non-Linear Index (MNLI)',...
        'Modified Simple Ratio (MSR)',...
        'Modified Soil Adjusted Vegetation Index 2 (MSAVI2)',...
        'Non-Linear Index (NLI)',...
        'Normalized Difference Vegetation Index (NDVI)',...
        'Optimized Soil Adjusted Vegetation Index (OSAVI)',...
        'Renormalized Difference Vegetation Index (RDVI)',...
        'Soil Adjusted Vegetation Index (SAVI)',...
        'Simple Ratio 650 (SR650)',...
        'Sum Green Index (SGI)',...
        'Transformed Difference Vegetation Index (TDVI)',...
        'Triangular Greenness Index (TGI)',...
        'Visible Atmospherically Resistant Index (VARI)',...
        'Wide Dynamic Range Vegetation Index (WDRVI)',...
        'WorldView Improved Vegetation Index (WV-VI)',...
        'Atmospherically Resistant Vegetation Index (ARVI)',...
        'Modified Chlorophyll Absorption Ratio Index (MCARI)',...
        'Modified Chlorophyll Absorption Ratio Index - Improved (MCARI2)',...
        'Modified Red Edge Normalized Difference Vegetation Index (MRENDVI)',...
        'Modified Red Edge Simple Ratio (MRESR)',...
        'Modified Triangular Vegetation Index (MTVI)',...
        'Modified Triangular Vegetation Index - Improved (MTVI2)',...
        'Red Edge Normalized Difference Vegetation Index (RENDVI)',...
        'Transformed Chlorophyll Absorption Reflectance Index (TCARI)',...
        'Triangular Vegetation Index (TVI)',...
        'Vogelmann Red Edge Index 1 (VREI1)',...
        'Vogelmann Red Edge Index 2 (VREI2)',...
        'Simple Ratio 680 (SR680)',...
        'Simple Ratio 705 (SR705)',...
        'Normalized Difference 680 (ND680)',...
        'Normalized Difference 705 (ND705)',...
        'Modified Simple Ratio 705 (mSR705)',...
        'Modified Normalized Difference 705 (mND705)',...
        'Plant Senescence Reflectance Index (PSRI)',...
        'Water Band Index (WBI)',...
        'Anthocyanin Reflectance Index 1 (ARI1)',...
        'Anthocyanin Reflectance Index 2 (ARI2)',...
        'Carotenoid Reflectance Index 1 (CRI1)',...
        'Carotenoid Reflectance Index 2 (CRI2)',...
        'Photochemical Reflectance Index (PRI)',...
        'Structure-Insensitive Pigment Index (SIPI)',...
        'Red Green Ratio Index (RGRI)',...
        'Burn Area Index (BAI)',...
        'Iron Oxide Ratio (IOR)',...
        'WorldView New Iron Index (WV-II)',...
        'WorldView Soil Index (WV-SI)',...
        'Normalized Difference Water Index (NDWI)',...
        'Normalized Difference Mud Index (NDMI)',...
        'WorldView Built-Up Index (WV-BI)',...
        'WorldView Non-Homogeneous Feature Difference (WV-NHFD)',...
        'WorldView Water Index (WV-WI)'};
    index_vnir_acc={'DVI',	'EVI',	'GEMI',	'GARI',	'GCI',	'GDVI',	'GLI',	'GNDVI',	'GOSAVI',	'GRVI',	'GSAVI',	'IPVI',	'LAI',	'MNLI',	'MSR',	'MSAVI2',	'NLI',	'NDVI',	'OSAVI',	'RDVI',	'SAVI',	'SR650',	'SGI',	'TDVI',	'TGI',	'VARI',	'WDRVI',	'WV-VI',	'ARVI',	'MCARI',	'MCARI2',	'MRENDVI',	'MRESR',	'MTVI',	'MTVI2',	'RENDVI',	'TCARI',	'TVI',	'VREI1',	'VREI2',	'SR680',	'SR705',	'ND680',	'ND705',	'mSR705',	'mND705',	'PSRI',	'WBI',	'ARI1',	'ARI2',	'CRI1',	'CRI2',	'PRI',	'SIPI',	'RGRI',	'BAI',	'IOR',	'WV-II',	'WV-SI',	'NDWI',	'NDMI',	'WV-BI',	'WV-NHFD',	'WV-WI'};
    associated_vnir_map={map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_green;	map_blue;	map_orange;	map_orange;	map_orange;	map_orange;	map_green;	map_green;	map_green;	map_gray;	map_red;	map_red;	map_red;	map_blue;	map_green;	map_gray;	map_gray;	map_blue};
    associated_vnir_p={p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_g;	p_b;	p_o;	p_o;	p_o;	p_o;	p_g;	p_g;	p_g;	p_gr;	p_r;	p_r;	p_r;	p_b;	p_g;	p_gr;	p_gr;	p_b};
    
    if nargin<5
        [indx,~] = listdlg('ListString',index_vnir,'ListSize',[300 450]);
    end
    
    for i=1:length(indx)
        index_vnir_selec{i}=index_vnir{indx(i)};
    end
    
    [~,Coastal]=find(abs(wl-425)==min(abs(wl-425)));
    [~,Blue]=find(abs(wl-470)==min(abs(wl-470)));
    [~,Green]=find(abs(wl-550)==min(abs(wl-550)));
    [~,Yellow]=find(abs(wl-605)==min(abs(wl-605)));
    [~,Red]=find(abs(wl-650)==min(abs(wl-650)));
    [~,Red_Edge]=find(abs(wl-725)==min(abs(wl-725)));
    [~,NIR]=find(abs(wl-860)==min(abs(wl-860)));
    [~,NIR1]=find(abs(wl-830)==min(abs(wl-830)));
    [~,NIR2]=find(abs(wl-950)==min(abs(wl-950)));
    
    %% Broadband Greenness
    
    if sum(double(strcmp('Difference Vegetation Index (DVI)',index_vnir_selec)))
        % Tucker, C. "Red and Photographic Infrared Linear Combinations for Monitoring Vegetation. Remote Sensing of Environment 8 (1979): 127–150
        iter=iter+1;
        
        IM=squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Red(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Enhanced Vegetation Index (EVI)',index_vnir_selec)))
        % Huete, A., et al. "Overview of the Radiometric and Biophysical Performance of the MODIS Vegetation Indices." Remote Sensing of Environment 83 (2002):195–213
        iter=iter+1;
        
        IM=2.5*(squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Red(1))))./...
            (squeeze(M(:,:,NIR(1)))+6*squeeze(M(:,:,Red(1)))-...
            7.5*squeeze(M(:,:,Blue(1)))+1);
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Global Environmental Monitoring Index (GEMI)',index_vnir_selec)))
        % Pinty, B., and M. Verstraete. GEMI: a Non-Linear Index to Monitor Global Vegetation From Satellites. Vegetation 101 (1992): 15-20
        iter=iter+1;
        
        eta=(2*(squeeze(M(:,:,NIR(1))).^2-squeeze(M(:,:,Red(1))).^2)+...
            1.5*squeeze(M(:,:,NIR(1)))+0.5*squeeze(M(:,:,Red(1))))./...
            (squeeze(M(:,:,NIR(1)))+squeeze(M(:,:,Red(1)))+0.5);
        IM=eta.*(1-0.25*eta)-((squeeze(M(:,:,Red(1)))-0.125)./...
            (1-squeeze(M(:,:,Red(1)))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        clear eta
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Green Atmospherically Resistant Index (GARI)',index_vnir_selec)))
        % Gitelson, A., Y. Kaufman, and M. Merzylak. "Use of a Green Channel in Remote Sensing of Global Vegetation from EOS-MODIS." Remote Sensing of Environment 58 (1996): 289-298
        iter=iter+1;
        
        IM=(squeeze(M(:,:,NIR(1)))-(squeeze(M(:,:,Green(1)))-...
            1.7*(squeeze(M(:,:,Blue(1)))-squeeze(M(:,:,Red(1))))))./...
            (squeeze(M(:,:,NIR(1)))+(squeeze(M(:,:,Green(1)))-...
            1.7*(squeeze(M(:,:,Blue(1)))-squeeze(M(:,:,Red(1))))));
        % 1.7 = constant define by ENVI
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Green Chlorophyll Index (GCI)',index_vnir_selec)))
        % Gitelson, A., Y. Gritz, and M. Merzlyak. "Relationships Between Leaf Chlorophyll Content and Spectral Reflectance and Algorithms for Non-Destructive Chlorophyll Assessment in Higher Plant Leaves." Journal of Plant Physiology 160 (2003): 271-282
        iter=iter+1;
        
        IM=(squeeze(M(:,:,NIR(1)))./squeeze(M(:,:,Green(1))))-1;
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Green Difference Vegetation Index (GDVI)',index_vnir_selec)))
        % Sripada, R., et al. "Determining In-Season Nitrogen Requirements for Corn Using Aerial Color-Infrared Photography." Ph.D. dissertation, North Carolina State University, 2005
        iter=iter+1;
        
        IM=squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Green(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Green Leaf Index (GLI)',index_vnir_selec)))
        % Louhaichi, M., M. Borman, and D. Johnson. "Spatially Located Platform and Aerial Photography for Documentation of Grazing Impacts on Wheat." Geocarto International 16, No. 1 (2001): 65-70
        iter=iter+1;
        
        IM=((squeeze(M(:,:,Green(1)))-squeeze(M(:,:,Red(1))))+...
            (squeeze(M(:,:,Green(1)))-squeeze(M(:,:,Blue(1)))))./...
            (2*squeeze(M(:,:,Green(1)))+squeeze(M(:,:,Red(1)))+...
            squeeze(M(:,:,Blue(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Green Normalized Difference Vegetation Index (GNDVI)',index_vnir_selec)))
        %Gitelson, A., and M. Merzlyak. "Remote Sensing of Chlorophyll Concentration in Higher Plant Leaves." Advances in Space Research 22 (1998): 689-692
        iter=iter+1;
        
        IM=(squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Green(1))))./...
            (squeeze(M(:,:,NIR(1)))+squeeze(M(:,:,Green(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Green Optimized Soil Adjusted Vegetation Index (GOSAVI)',index_vnir_selec)))
        %Sripada, R., et al. "Determining In-Season Nitrogen Requirements for Corn Using Aerial Color-Infrared Photography." Ph.D. dissertation, North Carolina State University, 2005
        iter=iter+1;
        
        IM=(squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Green(1))))./...
            (squeeze(M(:,:,NIR(1)))+squeeze(M(:,:,Green(1)))+0.16);
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Green Ratio Vegetation Index (GRVI)',index_vnir_selec)))
        %Sripada, R., et al. "Aerial Color Infrared Photography for Determining Early In-season Nitrogen Requirements in Corn." Agronomy Journal 98 (2006): 968-977
        iter=iter+1;
        
        IM=squeeze(M(:,:,NIR(1)))./squeeze(M(:,:,Green(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Green Soil Adjusted Vegetation Index (GSAVI)',index_vnir_selec)))
        %Sripada, R., et al. "Determining In-Season Nitrogen Requirements for Corn Using Aerial Color-Infrared Photography." Ph.D. dissertation, North Carolina State University, 2005
        iter=iter+1;
        
        IM=1.5*((squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Green(1))))./...
            (squeeze(M(:,:,NIR(1)))+squeeze(M(:,:,Green(1)))+0.5));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Infrared Percentage Vegetation Index (IPVI)',index_vnir_selec)))
        %Crippen, R. "Calculating the Vegetation Index Faster." Remote Sensing of Environment 34 (1990): 71-73.
        iter=iter+1;
        
        IM=squeeze(M(:,:,NIR(1)))./(squeeze(M(:,:,NIR(1)))+...
            squeeze(M(:,:,Red(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Leaf Area Index (LAI)',index_vnir_selec)))
        %Boegh, E., H. Soegaard, N. Broge, C. Hasager, N. Jensen, K. Schelde, and A. Thomsen. "Airborne Multi-spectral Data for Quantifying Leaf Area Index, Nitrogen Concentration and Photosynthetic Efficiency in Agriculture." Remote Sensing of Environment 81, no. 2-3 (2002): 179-193
        iter=iter+1;
        
        EVI=2.5*(squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Red(1))))./...
            (squeeze(M(:,:,NIR(1)))+6*squeeze(M(:,:,Red(1)))-...
            7.5*squeeze(M(:,:,Blue(1)))+1);
        IM=3.618*EVI-0.118;
        clear EVI
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Modified Non-Linear Index (MNLI)',index_vnir_selec)))
        %Yang, Z., P. Willis, and R. Mueller. "Impact of Band-Ratio Enhanced AWIFS Image to Crop Classification Accuracy." Proceedings of the Pecora 17 Remote Sensing Symposium (2008), Denver, CO
        iter=iter+1;
        
        IM=(squeeze(M(:,:,NIR(1))).^2*(1+0.5))./...
            (squeeze(M(:,:,NIR(1))).^2+squeeze(M(:,:,Red(1)))+0.5);
        % canopy background adjustment factor (L) value of 0.5
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Modified Simple Ratio (MSR)',index_vnir_selec)))
        %Chen, J. "Evaluation of Vegetation Indices and Modified Simple Ratio for Boreal Applications." Canadian Journal of Remote Sensing 22 (1996): 229-242
        iter=iter+1;
        
        IM=((squeeze(M(:,:,NIR(1)))./squeeze(M(:,:,Red(1))))-1)./...
            (sqrt(squeeze(M(:,:,NIR(1)))./squeeze(M(:,:,Red(1))))+1);
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Modified Soil Adjusted Vegetation Index 2 (MSAVI2)',index_vnir_selec)))
        %Qi, J., A. Chehbouni, A. Huete, Y. Kerr, and S. Sorooshian. "A Modified Soil Adjusted Vegetation Index." Remote Sensing of Environment 48 (1994): 119-126
        iter=iter+1;
        
        IM=(2*squeeze(M(:,:,NIR(1)))+1-...
            sqrt((2*squeeze(M(:,:,NIR(1)))+1).^2-...
            8*(squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Red(1))))))/2;
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Non-Linear Index (NLI)',index_vnir_selec)))
        %Goel, N., and W. Qin. "Influences of Canopy Architecture on Relationships Between Various Vegetation Indices and LAI and Fpar: A Computer Simulation." Remote Sensing Reviews 10 (1994): 309-347
        iter=iter+1;
        
        IM=(squeeze(M(:,:,NIR(1))).^2-squeeze(M(:,:,Red(1))))./...
            (squeeze(M(:,:,NIR(1))).^2+squeeze(M(:,:,Red(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Normalized Difference Vegetation Index (NDVI)',index_vnir_selec)))
        %Rouse, J., R. Haas, J. Schell, and D. Deering. Monitoring Vegetation Systems in the Great Plains with ERTS. Third ERTS Symposium, NASA (1973): 309-317
        iter=iter+1;
        
        IM=(squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Red(1))))./...
            (squeeze(M(:,:,NIR(1)))+squeeze(M(:,:,Red(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Optimized Soil Adjusted Vegetation Index (OSAVI)',index_vnir_selec)))
        %Rondeaux, G., M. Steven, and F. Baret. "Optimization of Soil-Adjusted Vegetation Indices." Remote Sensing of Environment 55 (1996): 95-107
        iter=iter+1;
        
        IM=(squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Red(1))))./...
            (squeeze(M(:,:,NIR(1)))+squeeze(M(:,:,Red(1)))+0.16);
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Renormalized Difference Vegetation Index (RDVI)',index_vnir_selec)))
        %Roujean, J., and F. Breon. "Estimating PAR Absorbed by Vegetation from Bidirectional Reflectance Measurements." Remote Sensing of Environment 51 (1995): 375-384
        iter=iter+1;
        
        IM=(squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Red(1))))./...
            sqrt(squeeze(M(:,:,NIR(1)))+squeeze(M(:,:,Red(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Soil Adjusted Vegetation Index (SAVI)',index_vnir_selec)))
        %Huete, A. "A Soil-Adjusted Vegetation Index (SAVI)." Remote Sensing of Environment 25 (1988): 295-309
        iter=iter+1;
        
        IM=(1.5*squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Red(1))))./...
            (squeeze(M(:,:,NIR(1)))+squeeze(M(:,:,Red(1)))+0.5);
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Simple Ratio 650 (SR650)',index_vnir_selec)))
        %Birth, G., and G. McVey. "Measuring the Color of Growing Turf with a Reflectance Spectrophotometer." Agronomy Journal 60 (1968): 640-643
        iter=iter+1;
        
        IM=squeeze(M(:,:,NIR(1)))./squeeze(M(:,:,Red(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Sum Green Index (SGI)',index_vnir_selec)))
        %Lobell, D., and G. Asner. "Hyperion Studies of Crop Stress in Mexico." Proceedings of the 12th Annual JPL Airborne Earth Science Workshop, Pasadena, CA (2003)
        iter=iter+1;
        
        wlr=[500 600];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        IM=mean(squeeze(M(:,:,ind1(1):ind2(1))),3);
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Transformed Difference Vegetation Index (TDVI)',index_vnir_selec)))
        %Bannari, A., H. Asalhi, and P. Teillet. "Transformed Difference Vegetation Index (TDVI) for Vegetation Cover Mapping" In Proceedings of the Geoscience and Remote Sensing Symposium, IGARSS '02, IEEE International, Volume 5 (2002)
        iter=iter+1;
        
        IM=1.5*((squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Red(1))))./...
            sqrt(squeeze(M(:,:,NIR(1))).^2+squeeze(M(:,:,Red(1)))+0.5));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Triangular Greenness Index (TGI)',index_vnir_selec)))
        % Hunt, E., C. Daughtry, J. Eitel, and D. Long. "Remote Sensing Leaf Chlorophyll Content Using a Visible Band Index." Agronomy Journal 103, No. 4 (2011): 1090-1099
        iter=iter+1;
        
        IM=(100*(squeeze(M(:,:,Red(1)))-squeeze(M(:,:,Blue(1))))-...
            180*(squeeze(M(:,:,Red(1)))-squeeze(M(:,:,Green(1)))))/2;
        % 100 = lred-lgreen
        % 180 = lref-lblue
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Visible Atmospherically Resistant Index (VARI)',index_vnir_selec)))
        % Gitelson, A., et al. "Vegetation and Soil Lines in Visible Spectral Space: A Concept and Technique for Remote Estimation of Vegetation Fraction. International Journal of Remote Sensing 23 (2002): 2537?2562
        iter=iter+1;
        
        IM=(squeeze(M(:,:,Green(1)))-squeeze(M(:,:,Red(1))))./...
            (squeeze(M(:,:,Green(1)))+squeeze(M(:,:,Red(1)))-...
            squeeze(M(:,:,Blue(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Wide Dynamic Range Vegetation Index (WDRVI)',index_vnir_selec)))
        % Gitelson, A. "Wide Dynamic Range Vegetation Index for Remote Quantification of Biophysical Characteristics of Vegetation." Journal of Plant Physiology 161, No. 2 (2004): 165-173.
        % Henebry, G., A. Viña, and A. Gitelson. "The Wide Dynamic Range Vegetation Index and its Potential Utility for Gap Analysis." Gap Analysis Bulletin 12: 50-56
        iter=iter+1;
        
        IM=(0.2*squeeze(M(:,:,NIR(1)))-squeeze(M(:,:,Red(1))))./...
            (0.2*squeeze(M(:,:,NIR(1)))+squeeze(M(:,:,Red(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('WorldView Improved Vegetation Index (WV-VI)',index_vnir_selec)))
        % Wolf, A. Using WorldView 2 Vis-NIR MSI Imagery to Support Land Mapping and Feature Extraction Using Normalized Difference Index Ratios. Unpublished report, Longmont, CO: DigitalGlobe (2010)
        iter=iter+1;
        
        IM=(squeeze(M(:,:,NIR2(1)))-squeeze(M(:,:,Red(1))))./...
            (squeeze(M(:,:,NIR2(1)))+squeeze(M(:,:,Red(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    %% Narrowband Greenness
    if sum(double(strcmp('Atmospherically Resistant Vegetation Index (ARVI)',index_vnir_selec)))
        %Kaufman, Y., and D. Tanre. "Atmospherically Resistant Vegetation Index (ARVI) for EOS-MODIS. IEEE Transactions on Geoscience and Remote Sensing 30, no. 2 (1992): 261-270
        iter=iter+1;
        wlr=[800 680 450];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind1(1)))-(squeeze(M(:,:,ind2(1)))-...
            1*(squeeze(M(:,:,ind3(1)))-squeeze(M(:,:,ind2(1))))))./...
            (squeeze(M(:,:,ind1(1)))+(squeeze(M(:,:,ind2(1)))-...
            1*(squeeze(M(:,:,ind3(1)))-squeeze(M(:,:,ind2(1))))));
        % 1 recommanded by ENVI
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Modified Chlorophyll Absorption Ratio Index (MCARI)',index_vnir_selec)))
        %Daughtry, C., et al. "Estimating Corn Leaf Chlorophyll Concentration from Leaf and Canopy Re?ectance." Remote Sensing Environment 74 (2000): 229–239
        iter=iter+1;
        wlr=[700 670 882];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=((squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))-...
            0.2*(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind3(1))))).*...
            (squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Modified Chlorophyll Absorption Ratio Index - Improved (MCARI2)',index_vnir_selec)))
        %Haboudane, D., et al. "Hyperspectral Vegetation Indices and Novel Algorithms for Predicting Green LAI of Crop Canopies: Modeling and Validation in the Context of Precision Agriculture." Remote Sensing of Environment 90 (2004): 337-352
        iter=iter+1;
        wlr=[800 670 550];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(1.5*(2.5*(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))-...
            1.3*(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind3(1))))))./...
            sqrt((2*squeeze(M(:,:,ind1(1)))+1).^2-...
            (6*squeeze(M(:,:,ind1(1)))-5*sqrt(squeeze(M(:,:,ind2(1)))))-...
            0.5);
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Modified Red Edge Normalized Difference Vegetation Index (MRENDVI)',index_vnir_selec)))
        %Datt, B. "A New Reflectance Index for Remote Sensing of Chlorophyll Content in Higher Plants: Tests Using Eucalyptus Leaves." Journal of Plant Physiology 154 (1999): 30-36.
        %Sims, D. and J. Gamon. "Relationships Between Leaf Pigment Content and Spectral Reflectance Across a Wide Range of Species, Leaf Structures and Developmental Stages." Remote Sensing of Environment 81 (2002): 337-354
        iter=iter+1;
        wlr=[750 705 445];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1)))-...
            2*squeeze(M(:,:,ind3(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Modified Red Edge Simple Ratio (MRESR)',index_vnir_selec)))
        %Sims, D., and J. Gamon. "Relationships Between Leaf Pigment Content and Spectral Reflectance Across a Wide Range of Species, Leaf Structures and Developmental Stages." Remote Sensing of Environment 81 (2002):337-354.
        %Datt, B. "A New Reflectance Index for Remote Sensing of Chlorophyll Content in Higher Plants: Tests Using Eucalyptus Leaves." Journal of Plant Physiology 154 (1999):30-36
        iter=iter+1;
        wlr=[750 445 705];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind3(1)))-squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Modified Triangular Vegetation Index (MTVI)',index_vnir_selec)))
        %Haboudane, D., et al. "Hyperspectral Vegetation Indices and Novel Algorithms for Predicting Green LAI of Crop Canopies: Modeling and Validation in the Context of Precision Agriculture." Remote Sensing of Environment 90 (2004): 337-352
        iter=iter+1;
        wlr=[800 550 670];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=1.2*(1.2*(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))-...
            2.5*(squeeze(M(:,:,ind3(1)))-squeeze(M(:,:,ind2(1)))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Modified Triangular Vegetation Index - Improved (MTVI2)',index_vnir_selec)))
        %Haboudane, D., et al. "Hyperspectral Vegetation Indices and Novel Algorithms for Predicting Green LAI of Crop Canopies: Modeling and Validation in the Context of Precision Agriculture." Remote Sensing of Environment 90 (2004): 337-352
        iter=iter+1;
        wlr=[800 550 670];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(1.5*(1.2*(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))-...
            2.5*(squeeze(M(:,:,ind3(1)))-squeeze(M(:,:,ind2(1))))))./...
            sqrt((2*squeeze(M(:,:,ind1(1)))+1).^2-...
            (6*squeeze(M(:,:,ind1(1)))-...
            5*sqrt(squeeze(M(:,:,ind3(1)))))-0.5);
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Red Edge Normalized Difference Vegetation Index (RENDVI)',index_vnir_selec)))
        %Gitelson, A., and M. Merzlyak. "Spectral Reflectance Changes Associated with Autumn Senescence of Aesculus Hippocastanum L. and Acer Platanoides L. Leaves." Journal of Plant Physiology 143 (1994): 286?292.
        %Sims, D., and J. Gamon. "Relationships Between Leaf Pigment Content and Spectral Reflectance Across a Wide Range of Species, Leaf Structures and Developmental Stages." Remote Sensing of Environment 81 (2002): 337-354
        iter=iter+1;
        wlr=[750 705];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Transformed Chlorophyll Absorption Reflectance Index (TCARI)',index_vnir_selec)))
        %Haboudane, D., et al. "Hyperspectral Vegetation Indices and Novel Algorithms for Predicting Green LAI of Crop Canopies: Modeling and Validation in the Context of Precision Agriculture." Remote Sensing of Environment 90 (2004): 337-352
        iter=iter+1;
        wlr=[700 670 550];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=3*((squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))-...
            0.2*(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind3(1)))).*...
            (squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Triangular Vegetation Index (TVI)',index_vnir_selec)))
        %Broge, N., and E. Leblanc. "Comparing Prediction Power and Stability of Broadband and Hyperspectral Vegetation Indices for Estimation of Green Leaf Area and Canopy Chlorophyll Density." Remote Sensing of Environment 76 (2000): 156-172
        iter=iter+1;
        wlr=[750 550 670];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(120*(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))-...
            200*(squeeze(M(:,:,ind3(1)))-squeeze(M(:,:,ind2(1)))))/2;
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Vogelmann Red Edge Index 1 (VREI1)',index_vnir_selec)))
        %Vogelmann, J., B. Rock, and D. Moss. "Red Edge Spectral Measurements from Sugar Maple Leaves." International Journal of Remote Sensing 14 (1993): 1563-1575
        iter=iter+1;
        wlr=[740 720];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        if ind1==ind2
            ind2=ind2+1;
        end
        
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Vogelmann Red Edge Index 2 (VREI2)',index_vnir_selec)))
        %Vogelmann, J., B. Rock, and D. Moss. "Red Edge Spectral Measurements from Sugar Maple Leaves." International Journal of Remote Sensing 14 (1993): 1563-1575
        iter=iter+1;
        wlr=[734 747 715 726];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        [~,ind4]=find(abs(wl-wlr(1,4))==min(abs(wl-wlr(1,4))));
        if ind4-ind3==3
            IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
                (squeeze(M(:,:,ind3(1)))+squeeze(M(:,:,ind4(1))));
            [a,b]=find(IM==inf);
            for j=1:length(a)
                if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                    IM(a(j),b(j))=NaN;
                else
                    IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
                end
            end
            [a,b]=find(IM==-inf);
            for j=1:length(a)
                if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                    IM(a(j),b(j))=NaN;
                else
                    IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
                end
            end
        else
            IM=zeros(size(M,1),size(M,2));
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Simple Ratio 680 (SR680)',index_vnir_selec))) %
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[800 680];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Simple Ratio 705 (SR705)',index_vnir_selec)))
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[750 705];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Normalized Difference 680 (ND680)',index_vnir_selec)))
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[800 680];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Normalized Difference 705 (ND705)',index_vnir_selec)))
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[750 705];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Modified Simple Ratio 705 (mSR705)',index_vnir_selec)))
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[750 445 705];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind3(1)))-squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Modified Normalized Difference 705 (mND705)',index_vnir_selec)))
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[750 705 445];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1)))-...
            2*squeeze(M(:,:,ind3(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Plant Senescence Reflectance Index (PSRI)',index_vnir_selec)))
        %Merzlyak, J., et al. "Non-destructive Optical Detection of Pigment Changes During Leaf Senescence and Fruit Ripening." Physiologia Plantarum 106 (1999): 135-141
        iter=iter+1;
        wlr=[680 500 750];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind3(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    %% Canopy Water Content
    
    if sum(double(strcmp('Water Band Index (WBI)',index_vnir_selec)))
        %Penuelas, J., et al. "The Reflectance at the 950-970 Region as an Indicator of Plant Water Status." International Journal of Remote Sensing 14 (1993): 1887-1905.
        %Champagne, C., et al. "Mapping Crop Water Status: Issues of Scale in the Detection of Crop Water Stress Using Hyperspectral Indices." Proceedings of the 8th International Symposium on Physical Measurements and Signatures in Remote Sensing, Aussois, France (2001): 79-84
        iter=iter+1;
        wlr=[970 900];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=squeeze(M(:,:,ind1(1)))./squeeze(M(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    %% Leaf Pigments
    
    if sum(double(strcmp('Anthocyanin Reflectance Index 1 (ARI1)',index_vnir_selec)))
        %Gitelson, A., M. Merzlyak, and O. Chivkunova. "Optical Properties and Nondestructive Estimation of Anthocyanin Content in Plant Leaves." Photochemistry and Photobiology 71 (2001): 38-45
        iter=iter+1;
        wlr=[550 700];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=1./squeeze(M(:,:,ind1(1)))-1./squeeze(M(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Anthocyanin Reflectance Index 2 (ARI2)',index_vnir_selec)))
        %Gitelson, A., M. Merzlyak, and O. Chivkunova. "Optical Properties and Nondestructive Estimation of Anthocyanin Content in Plant Leaves." Photochemistry and Photobiology 71 (2001): 38-45
        iter=iter+1;
        wlr=[800 550 700];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=squeeze(M(:,:,ind1(1))).*(1./squeeze(M(:,:,ind2(1)))-...
            1./squeeze(M(:,:,ind3(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Carotenoid Reflectance Index 1 (CRI1)',index_vnir_selec)))
        %Gitelson, A., et al. "Assessing Carotenoid Content in Plant Leaves with Reflectance Spectroscopy." Photochemistry and Photobiology 75 (2002): 272-281
        iter=iter+1;
        wlr=[510 550];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=1./squeeze(M(:,:,ind1(1)))-1./squeeze(M(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Carotenoid Reflectance Index 2 (CRI2)',index_vnir_selec)))
        %Gitelson, A., et al. "Assessing Carotenoid Content in Plant Leaves with Reflectance Spectroscopy." Photochemistry and Photobiology 75 (2002): 272-281
        iter=iter+1;
        wlr=[510 700];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=1./squeeze(M(:,:,ind1(1)))-1./squeeze(M(:,:,ind2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    %% Light Use Efficiency
    
    if sum(double(strcmp('Photochemical Reflectance Index (PRI)',index_vnir_selec)))
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[531 570];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Structure-Insensitive Pigment Index (SIPI)',index_vnir_selec)))
        %Sims, D.A., Gamon, J.A. (2002) Relationships between leaf pigment content and spectral reflectance across a wide range of species, leaf structures and developmental stages. Remote Sensing of Environment 81: 337–354
        iter=iter+1;
        wlr=[800 445 680];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind3(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Red Green Ratio Index (RGRI)',index_vnir_selec)))
        %Gamon, J., and J. Surfus. "Assessing Leaf Pigment Content and Activity With a Reflectometer." New Phytologist 143 (1999): 105-117
        iter=iter+1;
        wlr=[699 600 500];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=sum(squeeze(M(:,:,ind2(1):ind1(1))),3)./...
            sum(squeeze(M(:,:,ind3(1):ind2(1))),3);
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    %% Burn Indices
    
    if sum(double(strcmp('Burn Area Index (BAI)',index_vnir_selec)))
        %Chuvieco, E., M. Pilar Martin, and A. Palacios. “Assessment of Different Spectral Indices in the Red-Near-Infrared Spectral Domain for Burned Land Discrimination.” Remote Sensing of Environment 112 (2002): 2381-2396.
        %Martín, M. Cartografía e inventario de incendios forestales en la Península Iberica a partir de imágenes NOAA AVHRR. Doctoral thesis, Universidad de Alcalá, Alcalá de Henares (1998).
        iter=iter+1;
        
        IM=1./((0.1-squeeze(M(:,:,Red(1)))).^2+...
            (0.06-squeeze(M(:,:,NIR(1)))).^2);
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    %% Geology Indices
    
    if sum(double(strcmp('Iron Oxide Ratio (IOR)',index_vnir_selec)))
        %Segal, D. "Theoretical Basis for Differentiation of Ferric-Iron Bearing Minerals, Using Landsat MSS Data." Proceedings of Symposium for Remote Sensing of Environment, 2nd Thematic Conference on Remote Sensing for Exploratory Geology, Fort Worth, TX (1982): pp. 949-951.
        %Drury, S. Image Interpretation in Geology. London: Allen and Unwin (1987), 243 pp.
        iter=iter+1;
        
        IM=squeeze(M(:,:,Red(1)))./squeeze(M(:,:,Blue(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('WorldView New Iron Index (WV-II)',index_vnir_selec)))
        %Wolf, A. Using WorldView 2 Vis-NIR MSI Imagery to Support Land Mapping and Feature Extraction Using Normalized Difference Index Ratios. Unpublished report, Longmont, CO: DigitalGlobe (2010)
        iter=iter+1;
        
        IM=(squeeze(M(:,:,Green(1))).*squeeze(M(:,:,Yellow(1))))./...
            (squeeze(M(:,:,Green(1)))*1000);
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('WorldView Soil Index (WV-SI)',index_vnir_selec)))
        %Wolf, A. Using WorldView 2 Vis-NIR MSI Imagery to Support Land Mapping and Feature Extraction Using Normalized Difference Index Ratios. Unpublished report, Longmont, CO: DigitalGlobe (2010)
        iter=iter+1;
        
        IM=(squeeze(M(:,:,Green(1)))-squeeze(M(:,:,Yellow(1))))./...
            (squeeze(M(:,:,Green(1)))+squeeze(M(:,:,Yellow(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    %% Miscellaneous Indices
    
    if sum(double(strcmp('Normalized Difference Water Index (NDWI)',index_vnir_selec)))
        %McFeeters, S. "The use of Normalized Difference Water Index (NDWI) in the Delineation of Open Water Features." International Journal of Remote Sensing 17 (1996): 1425-1432
        iter=iter+1;
        
        IM=(squeeze(M(:,:,Green(1)))-squeeze(M(:,:,NIR(1))))./...
            (squeeze(M(:,:,Green(1)))+squeeze(M(:,:,NIR(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('Normalized Difference Mud Index (NDMI)',index_vnir_selec)))
        %Bernstein, L. S., X. Jin, B. Gregor, and S. Adler-Golden. "Quick Atmospheric Correction Code: Algorithm Description and Recent Upgrades." Optical Engineering 51, No. 11 (2012): 111719-1 to 111719-11
        iter=iter+1;
        wlr=[795 990];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(squeeze(M(:,:,ind1(1)))-squeeze(M(:,:,ind2(1))))./...
            (squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('WorldView Built-Up Index (WV-BI)',index_vnir_selec)))
        iter=iter+1;
        
        IM=(squeeze(M(:,:,Coastal(1)))-squeeze(M(:,:,Red_Edge(1))))./...
            (squeeze(M(:,:,Coastal(1)))+squeeze(M(:,:,Red_Edge(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('WorldView Non-Homogeneous Feature Difference (WV-NHFD)',index_vnir_selec)))
        %Wolf, A. Using WorldView 2 Vis-NIR MSI Imagery to Support Land Mapping and Feature Extraction Using Normalized Difference Index Ratios. Unpublished report, Longmont, CO: DigitalGlobe (2010).
        iter=iter+1;
        
        IM=(squeeze(M(:,:,Red_Edge(1)))-squeeze(M(:,:,Coastal(1))))./...
            (squeeze(M(:,:,Red_Edge(1)))+squeeze(M(:,:,Coastal(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    if sum(double(strcmp('WorldView Water Index (WV-WI)',index_vnir_selec)))
        %Wolf, A. Using WorldView 2 Vis-NIR MSI Imagery to Support Land Mapping and Feature Extraction Using Normalized Difference Index Ratios. Unpublished report, Longmont, CO: DigitalGlobe (2010).
        iter=iter+1;
        
        IM=(squeeze(M(:,:,Coastal(1)))-squeeze(M(:,:,NIR2(1))))./...
            (squeeze(M(:,:,Coastal(1)))+squeeze(M(:,:,NIR2(1))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_vnir{indx(iter)};
    end
    
    %% Display
    if sum(wl<2400&wl>1100)==0
    for i=1:size(C,3)
        ha=[];
        IM=squeeze(C(:,:,i));
        figure;
        ha(1)=subplot(3,1,1);
        imagesc(d,d(1:size(C,1)),IM)
        caxis([nanmean(IM(:))-3*nanstd(IM(:)) nanmean(IM(:))+3*nanstd(IM(:))]);
        colorbar
        map=associated_vnir_map{indx(i)};
        colormap(map)
        xlabel('Depth (cm)')
        ylabel('Width (cm)')
        set(gca,'fontsize',14)
        title(index_vnir{indx(i)})
        ha(2)=subplot(3,1,2);
        imagesc(d,d(1:size(IM,1)),RGB)
        colorbar
        xlabel('Depth (cm)')
        ylabel('Width (cm)')
        set(gca,'fontsize',14)
        ha(3)=subplot(313);
        plot(d,nanmean(IM(round(0.4*size(IM,1):0.6*size(IM,1)),:),1),'color',associated_vnir_p{indx(i)})
        grid on
        xlabel('Depth (cm)')
        ylabel('SI')
        set(gca,'fontsize',14)
        colorbar
        xlim([min(d) max(d)])
        linkaxes(ha,'x')
    end
    end
    
    for i=1:length(indx)
        lab_acc{i}=index_vnir_acc{indx(i)};
    end
    
    if sum(wl<2400&wl>1100)==0
        figure;
        imagesc(corr(reshape(C,[],size(C,3)),'rows','complete'))
        colorbar
        colormap(jet)
        xtickangle(45)
        set(gca,'fontsize',14,'xtick',1:length(lab_acc),'xticklabel',lab_acc,...
            'ytick',1:length(lab_acc),'yticklabel',lab_acc)
    end
end

if sum(wl<2400&wl>1100)>0
    
    index_swir={'Normalized Difference Nitrogen Index (NDNI)',...
        'Cellulose Absorption Index (CAI)',...
        'Lignin Cellulose Absorption Index (LCAI)',...
        'Normalized Difference Lignin Index (NDLI)',...
        'Clay Ratio (CR)',...
        'Sand moisture index (SMI)',...
        'SWIR Fine particles Index (FI)'};
    index_swir_acc={'NDNI';	'CAI';	'LCAI';	'NDLI';	'CR';	'SMI';	'FI'};
    associated_swir_map={map_blue;	map_green;	map_green;	map_green;	map_orange;	map_orange;	map_orange};
    associated_swir_p={p_b;	p_g;	p_g;	p_g;	p_o;	p_o;	p_o};
    
    if nargin<5
        [indx,~] = listdlg('ListString',index_swir);
    end
    
    for i=1:length(indx)
        index_swir_selec{i}=index_swir{indx(i)};
    end
    
    [~,SWIR1]=find(abs(wl-1650)==min(abs(wl-1650)));
    [~,SWIR2]=find(abs(wl-2220)==min(abs(wl-2220)));
    
    %% Canopy Nitrogen
    if sum(double(strcmp('Normalized Difference Nitrogen Index (NDNI)',index_swir_selec)))
        % Serrano, L., J. Penuelas, and S. Ustin. "Remote Sensing of Nitrogen and Lignin in Mediterranean Vegetation from AVIRIS Data: Decomposing Biochemical from Structural Signals." Remote Sensing of Environment 81 (2002):355-364.
        %Fourty, T., et al. "Leaf Optical Properties with Explicit Description of Its Biochemical Composition: Direct and Inverse Problems." Remote Sensing of Environment 56 (1996):104-117.
        iter=iter+1;
        wlr=[1510 1680];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(log(squeeze(M(:,:,ind1(1))))-log(squeeze(M(:,:,ind2(1)))))./...
            (log(squeeze(M(:,:,ind1(1))))+log(squeeze(M(:,:,ind2(1)))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    %% Dry or Senescent Carbon
    if sum(double(strcmp('Cellulose Absorption Index (CAI)',index_swir_selec)))
        % Daughtry, C. "Discriminating Crop Residues from Soil by Short-Wave Infrared Reflectance." Agronomy Journal 93 (2001): 125-131.
        %Daughtry, C., E. Hunt Jr., and J. McMurtrey III. "Assessing Crop Residue Cover Using Shortwave Infrared Reflectance." Remote Sensing of Environment 90 (2004): 126-134.
        iter=iter+1;
        wlr=[2000 2200 2100];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=0.5*(squeeze(M(:,:,ind1(1)))+squeeze(M(:,:,ind2(1))))-...
            squeeze(M(:,:,ind3(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Lignin Cellulose Absorption Index (LCAI)',index_swir_selec)))
        % Daughtry, C., E. Hunt, Jr., P. Doraiswamy, and J. McMurtrey III. "Remote Sensing the Spatial Distribution of Crop Residues." Agronomy Journal 97 (2005): 864-871.
        iter=iter+1;
        wlr=[2185 2225 2145 2295 2365];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        [~,ind4]=find(abs(wl-wlr(1,4))==min(abs(wl-wlr(1,4))));
        [~,ind5]=find(abs(wl-wlr(1,5))==min(abs(wl-wlr(1,5))));
        
        IM=100*((sum(squeeze(M(:,:,ind1(1):ind2(1))),3)-...
            sum(squeeze(M(:,:,ind3(1):ind1(1))),3))+...
            (sum(squeeze(M(:,:,ind1(1):ind2(1))),3)-...
            sum(squeeze(M(:,:,ind4(1):ind5(1))),3)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Normalized Difference Lignin Index (NDLI)',index_swir_selec)))
        % Serrano, L., J. Penuelas, and S. Ustin. "Remote Sensing of Nitrogen and Lignin in Mediterranean Vegetation from AVIRIS Data: Decomposing Biochemical from Structural Signals." Remote Sensing of Environment 81 (2002): 355-364.
        %Fourty, T., et al. "Leaf Optical Properties with Explicit Description of Its Biochemical Composition: Direct and Inverse Problems." Remote Sensing of Environment 56 (1996): 104-117.
        %Melillo, J., J. Aber, and J. Muratore. "Nitrogen and Lignin Control of Hardwood Leaf Litter Decomposition Dynamics." Ecology 63 (1982): 621-626.
        iter=iter+1;
        wlr=[1754 1680];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        
        IM=(log(1./squeeze(M(:,:,ind1(1))))-...
            log(1./squeeze(M(:,:,ind2(1)))))./...
            (log(1./squeeze(M(:,:,ind1(1))))+...
            log(1./squeeze(M(:,:,ind2(1)))));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    %% Geology Indices
    
    if sum(double(strcmp('Clay Ratio (CR)',index_swir_selec)))
        % Drury, S. Image Interpretation in Geology. London: Allen and Unwin (1987), 243 pp.
        iter=iter+1;
        
        IM=squeeze(M(:,:,SWIR1(1)))./squeeze(M(:,:,SWIR2(1)));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('Sand moisture index (SMI)',index_swir_selec)))
        iter=iter+1;
        wlr=[1860 1925 2140];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=sqrt((squeeze(M(:,:,ind1(1))).^2+...
            squeeze(M(:,:,ind2(1))).^2+squeeze(M(:,:,ind3(1))).^2)/3);
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    if sum(double(strcmp('SWIR Fine particles Index (FI)',index_swir_selec)))
        iter=iter+1;
        wlr=[2133 2225 2209];
        [~,ind1]=find(abs(wl-wlr(1,1))==min(abs(wl-wlr(1,1))));
        [~,ind2]=find(abs(wl-wlr(1,2))==min(abs(wl-wlr(1,2))));
        [~,ind3]=find(abs(wl-wlr(1,3))==min(abs(wl-wlr(1,3))));
        
        IM=(squeeze(M(:,:,ind2(1))).^2)./...
            (squeeze(M(:,:,ind2(1))).*(squeeze(M(:,:,ind2(1))).^3));
        [a,b]=find(IM==inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        [a,b]=find(IM==-inf);
        for j=1:length(a)
            if a(j)==1||a(j)==size(M,1)||b(j)==1||b(j)==size(M,2)
                IM(a(j),b(j))=NaN;
            else
                IM(a(j),b(j))=median(median(IM(a(j)-1:a(j)+1,b(j)-1:b(j)+1)));
            end
        end
        C(:,:,iter)=IM;
        label{iter}=index_swir{indx(iter)};
    end
    
    for i=1:size(C,3)
        ha=[];
        IM=squeeze(C(:,:,i));
        figure;
        ha(1)=subplot(3,1,1);
        imagesc(d,d(1:size(C,1)),IM)
        caxis([nanmean(IM(:))-3*nanstd(IM(:)) nanmean(IM(:))+3*nanstd(IM(:))])
        colorbar
        map=associated_swir_map{indx(i)};
        colormap(map)
        xlabel('Depth (cm)')
        ylabel('Width (cm)')
        set(gca,'fontsize',14)
        title(index_swir{indx(i)})
        ha(2)=subplot(3,1,2);
        imagesc(d,d(1:size(IM,1)),RGB)
        colorbar
        xlabel('Depth (cm)')
        ylabel('Width (cm)')
        set(gca,'fontsize',14)
        ha(3)=subplot(313);
        plot(d,nanmean(IM(round(0.4*size(IM,1):0.6*size(IM,1)),:),1),'color',associated_swir_p{indx(i)})
        grid on
        xlabel('Depth (cm)')
        ylabel('SI')
        set(gca,'fontsize',14)
        colorbar
        xlim([min(d) max(d)])
        linkaxes(ha,'x')
    end
    
    if sum(wl<800&wl>600)>0
        a=length(lab_acc);
        for i=a+1:a+length(indx)
            lab_acc{i}=index_swir_acc{indx(i-a)};
        end
    else
        for i=1:length(indx)
            lab_acc{i}=index_swir_acc{indx(i)};
        end
    end
    
    figure;
    imagesc(corr(reshape(C,[],size(C,3)), 'rows','complete'))
    colorbar
    colormap(jet)
    xtickangle(45)
    set(gca,'fontsize',14,'xtick',1:length(lab_acc),'xticklabel',lab_acc,'ytick',1:length(lab_acc),'yticklabel',lab_acc)
    
end

end