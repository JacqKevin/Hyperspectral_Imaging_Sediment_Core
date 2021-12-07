function [S,para1,para2]=pretreat(S,method,para1,para2)
%+++   data pretreatment
%+++ HD Li, Central South University

if nargin==2
    if size(S,3)>1
        d=1;
        d1=size(S,1);
        d2=size(S,2);
        S=reshape(S,[],size(S,3));
    else
        d=0;
    end
    [~,Nx]=size(S);
    if strcmp(method,'autoscaling')
        para1=mean(S,1);
        para2=std(S);
    elseif strcmp(method,'center')
        para1=mean(S,1);
        para2=ones(1,Nx);
    elseif strcmp(method,'unilength')
        para1=mean(S,1);
        for j=1:size(S,2)
            para2(1,j)=norm(S(:,j)-para1(j));
        end
    elseif strcmp(method,'minmax')
        para1=min(S);
        maxv=max(S);
        para2=maxv-para1;
    elseif strcmp(method,'pareto')
        para1=mean(S,1);
        para2=sqrt(std(S));
    end
    
    for i=1:Nx
        S(:,i)=(S(:,i)-para1(i))/para2(i);
    end
    
    if d==1
        S=reshape(S,d1,d2,size(S,2));
    end
    
elseif nargin==4
    [~,Nx]=size(S);
    for i=1:Nx
        S(:,i)=(S(:,i)-para1(i))/para2(i);
    end
end