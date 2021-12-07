function [S_cr,cr,wl]=continuum_removal(wl,S)
% Function to compute Continuum Removal

if length(wl)==size(wl,2)
    wl=wl';
end

if size(S,3)==1
    d=2;
else
    d=3;
end

if size(S,3)>1 % Spectra or Hyperspectral Cube
    M=S;
    S=reshape(M,size(M,1)*size(M,2),size(M,3));
end

if length(wl)==size(S,2)
    S=S';
end

S_cr=zeros(size(S,1),size(S,2));
cr=zeros(size(S,1),size(S,2));

% h=waitbar(0,'Continuum removal in process');
for i=1:size(S,2)
%     waitbar(i/size(S,2))
    
    %Spectrum triangulation
    DT = delaunayTriangulation(wl,S(:,i));
    
    %calculates the convex envelope of the spectrum
    k=convexHull(DT);
    
    %plot( wl,S)
    %hold on, plot(DT.Points(k,1),DT.Points(k,2),'.')
    bl=[DT.Points(k,1) DT.Points(k,2)];
    [~,k]=max(bl(:,1));
    
    cr(:,i)=interp1(bl(k:end,1),bl(k:end,2),wl);
    %hold on, plot(wl,cr)
    S_cr(:,i)=S(:,i)./cr(:,i);
end
% close(h)

S_cr=S_cr';
cr=cr';

if d==3
    S_cr=reshape(S_cr,size(M,1),size(M,2),size(M,3));
end
end