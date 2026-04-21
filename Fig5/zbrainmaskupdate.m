
load('D:\OneDrive - The University of Melbourne\Research\SCN1\region_shortlist_variables.mat')


combinelist = {[5:10],[11,40],[14:15],[32:33],[16:21],[49,52,53],[58:59],[68:77],[62,85:89],[94:100],[105:110],[111:114],[117:123],[42:48]};
templist = Zbrain_shortlist;
orderList = {'Telencephalon','Diencephalon','Mesencephalon','Rhombencephalon'};

% remove not wanted regions
for i = 1:length(combinelist)
    templist(combinelist{i},:) = {''};
end

% combine regions together
for i = 1:length(combinelist)
    newcol = Zbrain_shortlist(combinelist{i}(1),1:2);
    newcol{3} = vertcat(Zbrain_shortlist{combinelist{i},3});
    templist(length(templist)+1,:) = newcol;
end

finallist = [];
for i = 1:length(orderList)

    % move whole region index to the last
    allregind = find(ismember(templist(:,1),orderList{i}));
    rind = find(ismember(templist(:,1),orderList{i}) & strcmp(templist(:,2),' '));
    temprow = templist(allregind(end),:);
    templist(allregind(end),:) = templist(rind,:);
    templist(rind,:) = temprow;

    % reorder the cell file list
    finallist = [finallist;templist(ismember(templist(:,1),orderList{i}),:)];
end



Zbrain_shortlist = finallist;