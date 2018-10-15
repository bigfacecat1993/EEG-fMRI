threshold = 13;
top = 10;
F = GC;
G_all = digraph(F,'OmitSelfLoops');
G_all.Nodes.Var1(1:123,1) = name;
Z_bin = F > threshold;
Z_bin = double(Z_bin);
G_bin = digraph(Z_bin,'OmitSelfLoops');
G_bin.Nodes.Var1(1:123,1) = name;
indices = find(G_all.Edges.Weight<threshold);
G_rm = G_all;
G_rm = rmedge(G_rm,indices);


label = name;


figure('Name','WD Network');
p = plot(G_rm,'NodeLabel',label,'MarkerSize',5,'Layout','circle');
title('WD Network low freq');
%highlight(p,interest_channels,'Node  Color','green');

% n102 = neighbors(G_rm,102);
% n101 = neighbors(G_rm,101);
% highlight(p,n102,'NodeColor','red');
% highlight(p,n101,'NodeColor','yellow');
% highlight(p,102,n102,'EdgeColor','red','LineWidth',2);
% highlight(p,101,n101,'EdgeColor','red','LineWidth',2);

s = sparse(Z_bin);
 
figure('Name','Betweenness Centrality Scores');
betweenness = centrality(G_bin,'betweenness','Cost',G_bin.Edges.Weight);
p = plot(G_bin,'NodeLabel',label,'MarkerSize',5,'Layout','circle');
p.NodeCData = betweenness;
colormap(flip(autumn,1));
title('Betweenness Centrality Scores')
[betweenness_sorted, betweenness_sorted_ROI] = sort(betweenness,'descend');
betweenness_sorted_ROI_name = cell(123,1);
for i = 1:123
    betweenness_sorted_ROI_name{i,1} = name{betweenness_sorted_ROI(i),1};
end

% highlight(p,betweenness_sorted_ROI_low(1:top),'NodeColor','red');
% 
% for i = 1:top
%     nei = neighbors(G_low_bin, betweenness_sorted_ROI_low(i));
%     highlight(p, betweenness_sorted_ROI(i),nei,'EdgeColor','red','LineWidth',2);
% end


figure('Name','Clustering coefficient Scores');
C = full(clustering_coef_bd(s));
p = plot(G_rm,'NodeLabel',label,'MarkerSize',5,'Layout','circle');
p.NodeCData = C;
colormap(flip(autumn,1));
title('Clustering coefficient Scores')
[clustering_coeff_sorted, clustering_coeff_sorted_ROI] = sort(C,'descend');
clustering_coeff_sorted_ROI_name = cell(123,1);
for i = 1:123
   clustering_coeff_sorted_ROI_name{i,1} = name{clustering_coeff_sorted_ROI(i),1};
end
%highlight(p,clustering_coeff_sorted_ROI_low(1:top));
% nei = [];
% cluster_sorted_ROI_name = [];
% for i = 1:top
%     nei = neighbors(G_rm, clustering_coeff_sorted_ROI(i));
%     highlight(p, clustering_coeff_sorted_ROI(i),nei,'EdgeColor','red','LineWidth',2);
%     
% 
% end


figure('Name','Degree Scores');
deg = full(degrees_dir(s));
p = plot(G_bin,'NodeLabel',label,'MarkerSize',5,'Layout','circle');
p.NodeCData = deg;
colormap(flip(autumn,1));
title('Degree Scores')
[deg_sorted, deg_sorted_ROI] = sort(deg,'descend');
deg_sorted_ROI_name = cell(123,1);
for i = 1:123
   deg_sorted_ROI_name{i,1} = name{deg_sorted_ROI(i),1};
end
% highlight(p,deg_sorted_ROI(1:top));
% nei = [];
% deg_sorted_ROI_name = [];
% for i = 1:top
%     nei = neighbors(G_rm, deg_sorted_ROI(i));
%     highlight(p, deg_sorted_ROI(i),nei,'EdgeColor','red','LineWidth',2);
%    
% 
% end


% figure('Name','Distance matrix');
% distance = full(distance_bin(s));
% imagesc(distance);
% title('Distance matrix');

figure('Name','InCloseness Centrality Scores - Weighted');
incloseness = centrality(G_rm,'incloseness','Cost',G_rm.Edges.Weight);
p = plot(G_rm,'NodeLabel',label,'MarkerSize',5,'Layout','circle');
p.NodeCData = incloseness;
colormap(flip(autumn,1));
title('InCloseness Centrality Scores - Weighted');
[incloseness_sorted, incloseness_sorted_ROI] = sort(incloseness,'descend');
incloseness_sorted_ROI_name = cell(123,1);
for i = 1:123
   incloseness_sorted_ROI_name{i,1} = name{incloseness_sorted_ROI(i),1};
end
%highlight(p,incloseness_sorted_ROI_high(1:top));


% nei = [];
% close_sorted_ROI_name = [];
% for i = 1:top
%     nei = neighbors(G_rm, closeness_sorted_ROI(i));
%     highlight(p, closeness_sorted_ROI(i),nei,'EdgeColor','red','LineWidth',2);
%    
% end
% betweenness_sorted_ROI_name = cell(1,10);
% cluster_sorted_ROI_name=cell(1,10);
% deg_sorted_ROI_name= cell(1,10);
% close_sorted_ROI_name= cell(1,10);
% for i = 1:top
% betweenness_sorted_ROI_name{i} = names{betweenness_sorted_ROI(i)};
% cluster_sorted_ROI_name{i} = names{clustering_coeff_sorted_ROI(i)};
% deg_sorted_ROI_name{i} = names{deg_sorted_ROI(i)};
% close_sorted_ROI_name{i} = names{closeness_sorted_ROI(i)};
% end
% out = zeros(129,1);
% in = zeros(129,1);
% for i = 1:129
%  
%         out(i,1) = sum(inte_pcgc(i,:));
%         in(i,1) = sum(inte_pcgc(:,i));
% end

[out_in_sorted, out_in_sorted_ROI] = sort(out+in,'descend');
out_in_sorted_ROI_name = cell(123,1);
for i = 1:123
   out_in_sorted_ROI_name{i,1} = name{out_in_sorted_ROI(i),1};
end