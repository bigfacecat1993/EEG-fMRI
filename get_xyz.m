%Files=dir('D:\data\EPK185_rsfMRI_DW\20150514\WD-19_10mm_all_noconn\xyz\*.*');
Files = natsort(cellstr({Files.name}));
Nc = 129;
xyz = cell(Nc,1);
k = 1;
for k=1:Nc
    if k == 50 || k == 51
        xyz{k} = zeros(1,2);
        index = k;
        FileNames = Files(index);
    
        content = load(char(strcat('D:\data\EPK185_rsfMRI_DW\20150514\WD-19_10mm_all_noconn\xyz\',FileNames)));
    
        avg_x = round(median(content(:,1)));
        avg_y = round(median(content(:,2)));
        avg_z = round(median(content(:,3)));
        avg = [avg_x avg_y avg_z];
   
    
        [len,~] = size(content);
        content_xyz = zeros(len,1);
        
        for i = 1:len

            content_xyz(i) = norm(content(i,:) - avg);
        
        end
        [~, Index] = sort(content_xyz,'ascend');
   
        xyz{k}(1,1:3) = content(Index(1),1:3);
        xyz{k}(1,4) = Index(i);
    end
    
        
end