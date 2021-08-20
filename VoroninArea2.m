%function for calculating the normalised average Voronoi Cell area
function [MeanVAnorm, SDVAnorm]=VorArea(celldatx,celldaty,s1,s2,ps)
xy(:,1)=celldatx(:,1);
xy(:,2)=celldaty(:,1);
figure; [v,c]=voronoin(xy);

w=0;
for i=1:length(c)
    cells=c{i};
    k=0;
    for j=1:length(cells)
        if v(cells(j),1)>0 && v(cells(j),2)>0 && v(cells(j),1)<s2 && v(cells(j),2)<s1
            k=k+0;
        else
            k=k+1;
        end
    end
        if k==0
            w=w+1;
            c2(w)=i;
            vorcell(w)=cell(c(i));
            
        end
end
vorcell=vorcell';
 
A = zeros(length(vorcell),1) ;
for i = 1:length(vorcell)
    v1 = v(vorcell{i},1) ; 
    v2 = v(vorcell{i},2) ;
    patch(v1,v2,rand(1,3))
    A(i) = polyarea(v1,v2) ;
    A(i)=A(i)*ps*ps;
end
MeanVA=mean(A);
MeanVAnorm=MeanVA/(s1*s2*ps*ps)
SDVA=std(A);
SDVAnorm=SDVA/(s1*s2*ps*ps)