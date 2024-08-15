compatible=zeros(1,3);
count=1;
iter=3;
for i=0:iter:9
    for j=0:iter:9
        for k=0:iter:9
        if(i+j+k==9)
            compatible(count,1)=i;
            compatible(count,2)=j;
            compatible(count,3)=k;
        count=count+1;
        end
            
        end
    end
end
figure
scatter3(compatible(:,1),compatible(:,2),compatible(:,3),'r','filled');
%line(compatible(:,1),compatible(:,2),compatible(:,3));
xlabel('Kupffer Cells');
ylabel('Stellate Cells');
zlabel('Tumor Cells');