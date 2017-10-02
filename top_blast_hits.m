function acc_num = top_blast_hits(accession,N)
%Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 
gb_data = getgenbank(accession);
[requestID, requestTime] = blastncbi(gb_data.Sequence,'blastn');
blast_data = getblast(requestID,'WaitTime',requestTime*1.5);
num = '';
for ii = 1:N
    name = blast_data.Hits(ii).Name;
    pos=strfind(name,'|');
    num = [num,',',name((pos(3)+1):(pos(4))-1)];
end
acc_num = num(2:end);
end