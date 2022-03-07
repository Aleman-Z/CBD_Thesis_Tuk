
function [p,q,cont,sig_pq]=convert2fieldtrip(ti,Mono,k,Mx,V,ro)
cont=0;
if ~isempty(Mx{k})
    for j=1:length(Mx{k})
%     ts=find(ti{k}==S{k}(j));
%     tend=find(ti{k}==E{k}(j));
%     sig{j}=Mono{k}(ts:tend);
    
   if nargin>5
    %Ripple-centered window.
    tm=find(ti{k}==Mx{k}(j));
        if  tm+ro<=length(ti{k}) && tm-ro>=1
%             p{j}=[V{k}(tm-ro:tm+ro).';V2{k}(tm-ro:tm+ro).';V3{k}(tm-ro:tm+ro).'];
%             q{j}=[Mono{k}(tm-ro:tm+ro).';Mono2{k}(tm-ro:tm+ro).';Mono3{k}(tm-ro:tm+ro).'];
            p{j}=[V{k}(tm-ro:tm+ro)];
            q{j}=[Mono{k}(tm-ro:tm+ro)];
            
            %sig_pq{j}=Mono{k}(ts:tend);
            
        else
            p{j}=[];
            q{j}=[];
            sig_pq{j}=[];
            cont=cont+1;
        end
   else
    p{j}=[];
    q{j}=[];
    sig_pq{j}=[];
   end
   
    end
else
    sig=[];
    p=[];
    q=[];
    sig_pq=[];
end
end






