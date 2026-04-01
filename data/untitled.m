x_unique = sort(unique(gridPts(:,1))); 
y_unique = sort(unique(gridPts(:,2)));

[~, x_idx] = ismember(gridPts(:,1), x_unique);
[~, y_idx] = ismember(gridPts(:,2), y_unique);

cost_grid = zeros(length(y_unique), length(x_unique));

cost_grid(sub2ind(size(cost_grid), y_idx, x_idx)) = cost;

costMap = figure('Color','white','Position',[100,100,1000,700]);

imagesc(x_unique, y_unique, cost_grid);
axis xy;  
shading flat;

colormap(parula); 
colorbar;          
xlabel('k_h','FontSize',12);
ylabel('k_k','FontSize',12);
title('Cost Map')

axis equal tight;  
grid off;          
set(gca, 'FontSize',10); 

savefig(costMap, 'Cost Map');% save .fig


%%
costMap3D = figure('Color','white');
scatter3(gridPts(:,1), gridPts(:,2), cost, 5, cost, 'filled');
xlabel('k_h'); ylabel('k_k'); zlabel('Cost');
colorbar; grid on; view(30,30);
title('Cost Map');

savefig(costMap3D, 'Cost Map 3D');% save .fig
