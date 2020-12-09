# === Lab 4: SNP Characterization, Analysis, and Discovery ===
# author: Nilesh Domah



# === Part 1: Analysis of SNP Variations in a Population ===
# 1.1:
# Open and read in chromosomal sequences in slim_chr2_seq.fasta for each of the 90 individuals.
infile = open('slim_chr2_seq.fasta', 'r')

# Read Reference header (1st line) and Sequence (2nd line)
tmp = infile.readline() # Reading first line
Ref_Header = tmp.strip() # Stripping any starting or/and ending whitespaces
tmp = infile.readline() # Reading second line
Ref_Seq = tmp.strip()

# This will read all lines (starting from line 3 until end)
# and save into two lists
Header = []
Seq = []
for tmp in infile:
    if tmp[0] == '>': # if tmp.find('>') != -1 # This will work too 
        Header.append(tmp.strip())
    else:
        Seq.append(tmp.strip())

# We don't need infile anymore (the data are stored), so close it
infile.close()


# 1.2:
# Open and read in SNP information stored in slim_chr2_SNPS.vcf (which can be opened as a text file).
infile2 = open("slim_chr2_SNPS.vcf", "r")
header = infile2.readline()
info = infile2.readlines()

# Store the information in this file in seven different lists, one for each column of data 
# (open the file in a text editor before writing your code to understand its format).
chrom_num = []
rsid = []
slim_pos = []
genome_pos = []
ref = []
alt = []
gene_info = []

for line in info:
    chrom_num_temp, rsid_temp, slim_pos_temp, genome_pos_temp, ref_temp, alt_temp, gene_info_temp = line.split(sep = "\t") 
    chrom_num.append(chrom_num_temp)
    rsid.append(rsid_temp)
    slim_pos.append(slim_pos_temp)
    genome_pos.append(genome_pos_temp)
    ref.append(ref_temp)
    alt.append(alt_temp)
    gene_info.append(gene_info_temp)
    
# Remember to close any files you open!
infile2.close()


# 1.3:
# For each provided variant, calculate the frequency of that SNP in the population of 90 individuals.
# In other words, how many of the 90 individuals have the corresponding base at the specified location for that variant?
# Store these variant frequencies in a separate list.
# 2.1 added below
NHfile = open("lymphoma_variants.txt", "w")
print("Individual ID", "variantID", "GeneDiseaseInfo", sep = "\t", file = NHfile)

var_freq = []

for i, pos in enumerate(slim_pos):
    real_pos = int(pos) - 1
    var = alt[i]
    count = 0
    for base in Seq:
        if (base[real_pos] == var):
            count += 1
            if "Non-Hodgkin Lymphoma" in gene_info:
                IndID = rsid
                varID = alt
                GDinfo = gene_info
                print(IndID, varID, GDinfo, sep = "\t", file = NHfile)
    freq = count / len(Seq)
    var_freq.append(freq)
    
NHfile.close()


# 1.4: 
# Find and print the ids of the SNPs in the population that occur with the highest and lowest frequency.
print("RSIDs\tfrequency")
for i, ID in enumerate(rsid):
    if var_freq[i] == max(var_freq):
        print("Highest Frequency: ", max(var_freq), "ID: ", ID)
    if var_freq[i] == min(var_freq):
        print("Lowest Frequency: ", min(var_freq), " ID: ", ID)
    

# 1.5:
# Open a new output file called slim_chr2_SNPs_withfrequencies.vcf and print the following information for each SNP:
#       Chromosome Number
#       SNP Name
#       SNP “Slim” Position
#       SNP Reference SNP
#       SNP Alternative SNP
#       SNP Frequency
#       Gene Disease Information
# Your file should be tab delimited, with the information for one SNP on each line. 
# Please include column titles, as you did in Lab 3! 
# You are provided with a sample output file called ‘slim_chr2_SNPs_withfrequencies_small.vcf’ to show the expected format.
writefile = open("slim_chr2_SNPs_withfrequencies_small.vcf", "w")
print("Chromosome Number", "SNP Name", 'SNP "Slim" Position', "SNP Reference SNP", "SNP Alternate SNP", "SNP Frequency", "Gene Disease Info", sep = "\t", end = "\n", file = writefile)

for i in range(len(chrom_num)):
    print(chrom_num[i], rsid[i], slim_pos[i], genome_pos[i], ref[i], alt[i], gene_info[i], sep = "\n", file = writefile)
    
writefile.close()

# Calculate mean of a list (input: a list, output: that list’s mean)
# Calculate the mean across all SNP frequencies in the population.
def mean_list(list_num):
    # Input: list_num: a list of numbers
    total_sum = 0
    for data in list_num:
        total_sum += data

    return total_sum / len(list_num)


# === Part 2: Characterizing SNPs ===
syn_freqs = []
for i, info in enumerate(gene_info):
    if 'synonymous' in info.split('|'):
        syn_freqs.append(freq[i])
        
nonsyn_freqs = []
for i, info in enumerate(gene_info):
    if 'nonsynonymous' in info.split('|'):
        nonsyn_freqs.append(freq[i])
        
print('Average Frequency (synonymous):', mean_list(syn_freqs))
print('Average Frequency (non-synonymous):', mean_list(nonsyn_freqs))



# === Part 3: Discovering New Variants ===
# 3.1:
# For each chromosome in the population, scan the nucleotides, comparing each base against the reference sequence.
dict_SNP = {}

for i, seq_data in enumerate(Seq):
    for j in range(len(seq_data)):
        pos = j + 1

        # Finding SNPs
        if seq_data[j] != Ref_Seq[j]:
            key = str(pos) + ':' + seq_data[j]

            ## 2. If you find a SNP (index where the individual’s sequence differs 
            #  from the reference), store its position in the sequence as well as 
            #  the variant base.
            
            ## Way 1: Storing list of Individual Headers as value
            #if not key in dict_SNP:
            #    dict_SNP[key] = []
            #dict_SNP[key].append(Header[i])
            
            ## Way 2: Storing Count of Individuals as value (need integer)
            if not key in dict_SNP.keys():
                dict_SNP[key] = 1
            else:
                dict_SNP[key] += 1
                

# 3.2:
# If you find a SNP (index where the individual’s sequence differs from the reference), 
# store its position in the sequence as well as the variant base.
# 3.3:
# Iterate over your newly discovered SNPs and print out the following information to a file named novel_variants.txt 
# (consult to novel_variants_small.txt for the formatting):
#       The SNP’s position in the sequence
#       The reference SNP
#       The alternative SNP
#       The frequency of the alternative SNP in this population
novel_out = open('novel_variants.txt', 'wt')
print ('Position', 'Ref base', 'Alt base', 'SNP frequency', sep='\t', file = novel_out)

positions = []

for SNP, count in dict_SNP.items():
    novel_pos, novel_alt = SNP.split(':')
    novel_pos = int(novel_pos)
    novel_ref = Ref_Seq[int(novel_pos)-1] # Getting Ref base using position
    novel_freq =  count / len(Seq)
    print(novel_pos, novel_ref, novel_alt, novel_freq, sep = '\t', file = novel_out)
    positions.append(novel_pos)
    
novel_out.close()
