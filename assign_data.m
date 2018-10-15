% extract global mean raw time-series data
% load('REX.mat')
% for s = 1:132
% load('REX.mat');
% params.s = s;params.gui = 0;
% rex(params); clear;
% end
%after manualy move the data into the folder below
 Nt = 600;
 Ntr = 1; %1 trials;
 Nc = 129; %130 ROIs
 X = zeros(Nt,Ntr,Nc);
 Files=dir('D:\data\EPK185_rsfMRI_DW\20150514\WD-19_10mm_all_noconn\data\*.*');
 xyz = load('D:\data\EPK185_rsfMRI_DW\20170315\WD-23_10mm_all_noconn\xyz\xyz.mat');
 names = cell(Nc,1);
 xyz = xyz.xyz;
for k=1:129
    
  
    K = int2str(k);
    n = k;
    FileNames = Files(n);
    %names{k} = FileNames;
    content = load(char(strcat('D:\data\EPK185_rsfMRI_DW\20150514\WD-19_10mm_all_noconn\data\',FileNames)));
    index = xyz{k}(:,4);
%     [~, voxel] = size(content.R);
%     unit = floor(voxel/Ntr);
    
    m = content.R;    

    for j = 1:1
        X(:,j,k) = m(:,index(j));
    end
 
  
end
%X= cell2mat(X);

% 
% X_or = cell(Nc,1);
% names_or = cell(Nc,1);
% for j = 1:Nc
%     B = all(X(:,:,j));
%    X_or{j} = zeros(Nt,Ntr);
%     if all(B) == 1
%        X_or{j}(:,:) = X(:,:,j);
%        %X_or{j}(X_or{j}==0) = [];
%        names_or{j} = names{j};
%     end
% end
% %manually delete empty names_or cell: WD23: 6 ROI is zero
% 
% X_or = X_or(~cellfun('isempty',X_or));  
% X_or = cell2mat(X_or);
% X_or = reshape(X_or,[600,4,131]);
% p = 2;
% fs = 1/0.72;
% % freq = linspace(0,fs/2,Nt);
% freq = fs/Nt; 
% %M = cell(length(names_or));

%            
% %plot(M_p.freq,M_p.gc);hold on,  plot(M.freq,M.gc);title(sprintf('Causality [%s] -> [%s]:\n',num2str(i1),num2str(i2))); 

%de-time_mean and detrend
fs = 500;
%X_reshape =  reshape(X,size(X,1)/8,8*size(X,2),size(X,3));
time = 0:1/fs:(size(X_reshape,1)-1)/fs;
X = new_x(:,128:140);
tmeanX = mean(X,1);
tmeanX = repmat(tmeanX,size(X,1),1);
X = X-tmeanX;
Xr_detrend = zeros(size(X));

for i = 1:Nc
    Xr_detrend(:,:,i) = detrend(Xr(:,:,i));
end

% % %de-emean
% % emeanXr = mean(Xr,3);
% % Xr = permute(Xr,[3 1 2]);
% % emeanXr = repmat(mean(Xr,1),size(Xr,1),1);
% % Xrr = Xr-emeanXr;
% % Xrr = permute(Xrr,[2,3,1]);
% 
% tot = zeros(Nc,1);
% dif = zeros(Nc,1);
 Xr_p = permute(Xr_detrend,[3 1 2]);
 X_cell = cell(Nc,1);
 for i = 1:Nc
     X_cell{i} = Xr_detrend(:,:,i);
 end
 
 M_p = cell(Nc);
 p = 2;
 freq=0.001:0.001:fs/2;
 freq_low = 0.001:0.0005:0.1;
for i1 = 1:Nc
    for i2 = 1:Nc
        if i1 ~= i2
%         
%         M{i1,i2} = compute_npCGCblock(X_or,fs,i1,i2);
%        
            
                    M_p{i1,i2} = compute_pCGCblock(X_cell,fs,freq_low,p,i1,i2);
                    S = sprintf('%d and %d has done! \n',i1,i2);
                    fprintf(S);
        end
       
    end
end

% 
morder = zeros(10,1);
ii = randi(98,1,90);
for i = 1:90
   [p, powdiff, aic, bic, mdiff] = find_morder(Xr_detrend(:,:,ii(i):ii(i)+1),fs);
   morder(i) = mdiff;
end


 
% 
% for i = 1:Nc
%     x = Xr_p(i,:,:);
%     ii = setdiff(1:Nc,i);
%     y = Xr_p(ii,:,:);
%     
%    [blc] = block_connectivity2(x,y,p,fs,freq);
%    dif_ = blc.Fx2y-blc.Fy2x; tot_ = blc.Fx2y+blc.Fy2x;
%    [~, idif] = find(isnan(dif_));
%    [~, itot] = find(isnan(tot_));
%    tdif = trapz(blc.f(1:(idif(1)-1)),dif_(1:(idif(1)-1)))*2/fs;
%    ttot = trapz(blc.f(1:(itot(1)-1)),tot_(1:(itot(1)-1)))*2/fs;
%    tot(i) = ttot; 
%    dif(i) = tdif;
% end
% p = 12;
% freq_low = 0:0.001:0.2;
% freq_high = 50:1:300;
% fs = 1024;
% F_high = zeros(Nc);
% for i = 1:Nc
%     for j = 1:Nc
%         if i ~= j
%             x = Xr_p(i,:,:);
%     
%             y = Xr_p(j,:,:);
%     
%             [blc] = block_connectivity2(x,y,p,fs,freq_high);
%             xy = trapz(blc.f,blc.Fx2y);
%             yx = trapz(blc.f,blc.Fy2x);
%             F_high(i,j) = xy;
%             F_high(j,i) = yx;
%         end
%     end
% end
colormap('jet');imagesc(GC);set(gca, 'XTick', [1:112]);
set(gca, 'XTickLabel',new_name(1:112));xtickangle(90);set(gca, 'FontSize', 6);set(gca, 'YTick', [1:112]);
set(gca, 'YTickLabel', new_name(1:112));

colormap('jet');bar(out+in);set(gca, 'XTick', [1:112]);
set(gca, 'XTickLabel', new_name);xtickangle(90);set(gca, 'FontSize', 9);