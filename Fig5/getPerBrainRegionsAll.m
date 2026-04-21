function [ PerBrainRegions ] = getPerBrainRegionsAll(Zbrain_Masks, ROI_correct, regionname)

% Need to make sure ROIs are appropriately rotated to match the zbrain
% masks

if nargin>2
    Index = contains(Zbrain_Masks(:,1),regionname);
    Zbrain_Masks = Zbrain_Masks(Index,:);
end

% Round just in case the ROIs weren't already rounded
ROI_correct = round(ROI_correct);

%% Assign all ROIs to a region
PerBrainRegions=struct();
for i=1:length(Zbrain_Masks)
    subRegionName = ['r' Zbrain_Masks{i,2}];

    regionName = Zbrain_Masks{i,1};
    Mask=Zbrain_Masks{i,3};
    Mask=unique(Mask,'rows');
    IsInBrainRegion=ismember(ROI_correct,Mask,'rows');
    subRegionName = strrep(subRegionName,' ','_');
    subRegionName = strrep(subRegionName,'-','_');
    subRegionName = strrep(subRegionName,'.','_');
    subRegionName = strrep(subRegionName,'(','');
    subRegionName = strrep(subRegionName,')','');
    subRegionName = strrep(subRegionName,',','');
    subRegionName = strrep(subRegionName,"'",'_');
    PerBrainRegions.(regionName).(subRegionName).idx=find(IsInBrainRegion==1);


end