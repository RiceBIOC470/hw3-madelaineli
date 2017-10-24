% GB comments
1.	100
2a. 70 Fraction is incorrectly calculated. Need to calculate the total number of aligned nucleotides and divide this number by the total length of the coding sequence form either gb_1.Sequence or gb_2.Sequence. 
2b. 70 Same issue as 2a
2c. 70 Same issue as 2a. 
3a 100 
3b. 100
3c. 100  	
Overall: 87


%HW3

%% Problem 1 - Smith-Waterman alignment
% Consider two sequences 'GTAATCC' and 'GTATCCG'
% Construct the scoring matrix for this with the parameters:
% match value = 2, mismatch value = -1, and gap penalty = -1. Use your
% solution to get the optimal alignment. If you prefer, it is acceptable to do this with
% pencil and paper, you can then take a snapshot of your solution and
% include it in your repository. 

img = imread('SW.JPG');
imshow(img)

%% Problem 2 - using the NCBI databases and sequence alignments

% Erk proteins are critical signal transducers of MAP kinase signaling.
% Accessions numbers for ERK1 (also called MAPK3) and ERK2 (also called MAPK1) human mRNA are NM_002746 and
% NM_002745, respectively. 

% Part 1. Perform an alignment of the coding DNA sequences of ERK1 and
% ERK2. What fraction of base pairs in ERK1 can align to ERK2? 
gb_1 = getgenbank('NM_002746');
gb_2 = getgenbank('NM_002745');

[score, align, start] = swalign(gb_1.Sequence,gb_2.Sequence,'Alphabet','nt','Showscore',false);
mismatch = sum(double(isspace(align(2,:))));
total = size(align,2);
frac = 1- (mismatch / total);
disp('fraction of bp in ERK1 aligned to ERK2')
disp(frac)
% Part2. Perform an alignment of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?

aa_1 = getgenpept(gb_1.CDS.protein_id);
aa_2 = getgenpept(gb_2.CDS.protein_id);

[score_aa, align_aa, start_aa] = swalign(aa_1.Sequence,aa_2.Sequence,'Showscore',false);
match_aa = count(align_aa(2,:),'|');
frac_aa = match_aa / size(align_aa,2);
disp('fraction of aa sequence of ERK1 aligned to ERK2')
disp(frac_aa)
% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they?
mm_1 = getgenbank('BC029712'); %not sure
mm_2 = getgenbank('BC058258'); %correct
[score_m1, align_m1, start_m1] = swalign(mm_1.Sequence,gb_1.Sequence,'Alphabet','nt','Showscore',false);
mismatch_m1 = sum(double(isspace(align_m1(2,:))));
frac_m1 = 1- (mismatch_m1 / size(align_m1,2));
disp('fraction of mouse ERK1 aligned to human ERK1')
disp(frac_m1)

[score_m2, align_m2, start_m2] = swalign(mm_2.Sequence,gb_2.Sequence,'Alphabet','nt','Showscore',false);
mismatch_m2 = sum(double(isspace(align_m2(2,:))));
frac_m2 = 1-(mismatch_m2 / size(align_m2,2));
disp('fraction of mouse ERK2 aligned to human ERK2')
disp(frac_m2)

mp_1 = getgenpept(mm_1.CDS.protein_id);
mp_2 = getgenpept(mm_2.CDS.protein_id);
[score_mp1, align_mp1, start_mp1] = swalign(mp_1.Sequence,aa_1.Sequence,'Alphabet','nt','Showscore',false);
match_mp1 = count(align_mp1(2,:),'|');
frac_mp1 = match_mp1/size(align.mp1,2);
disp('fraction of mouse ERK1 aa sequce aligned to human ERK1 aa sequence')
disp(frac_mp1)

[score_mp2, align_mp2, start_mp2] = swalign(mp_2.Sequence,aa_2.Sequence,'Alphabet','nt','Showscore',false);
match_mp2 = count(align_mp2(2,:),'|');
frac_mp2 = match_mp2/size(align.mp2,2);
disp('fraction of mouse ERK1 aa sequce aligned to human ERK1 aa sequence')
disp(frac_mp2)
% 
%% Problem 3: using blast tools programatically

% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 

% written in top_blast_hits.m

% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 

% written in hum_nonhum_match.m

% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 

% human gene
[hum,nonhum] =hum_nonhum_match('NM_000372')
% both human and non human match are found, but the human match is a
% variant of the gene given. the non human match, on the other hand, comes
% from Pan paniscus. Thus, we conclude that from the information given, Pan
% paniscus cdk7 gene has closer resemblance to human cdk7 gene than other 
% species in the ncbi database.
% nonhuman gene
[hum_2,nonhum_2] = hum_nonhum_match('NM_009874')
% no human match found within the top 20 hits, which means that human is
% probably less resembling to this mouse gene than other organisms
% available in the gene database. 
