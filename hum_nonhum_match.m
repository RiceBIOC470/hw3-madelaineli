function [hum,nonhum] = hum_nonhum_match(accession)
% Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 
nums = top_blast_hits(accession, 20);
pos = [0,strfind(nums,',')];
hum = '';
nonhum = '';
for ii = 1:(length(pos)-1)
    num = nums((pos(ii)+1):(pos(ii+1)-1));
    data = getgenbank(num);
    if strcmp(data.Source,'Homo sapiens (human)')
        hum = strcat(hum,',',num);
    else
        nonhum = strcat(nonhum,',',num);
    end
end
if ~isempty(hum)
    hum_pos = strfind(hum,',');
    hum = hum(2:(hum_pos(2)-1));
else
    hum = 'no human match found';
end
if ~isempty(nonhum)
    nonhum_pos = strfind(nonhum,',');
    nonhum = nonhum(2:(nonhum_pos(2)-1));
else
    nonhum = 'no nonhuman match found';
end
end