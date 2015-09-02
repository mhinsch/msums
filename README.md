# msums

Program for efficient computation of summary statistics from ms outputs including any arbitrary number of populations.
___

## Description of statistics
### Manual conventions
- The suffixes 'avg' and 'std' stand for average and standrad deviation, respectively.
- Populations are referred by indexes *i* and *j* and can range from 0 to n, where n is the number of sampled populations.
  - As an example 'FST_*i*x*j*_mean' refers to the mean Fst across loci between populations *i* and *j*.
- The coding of ancestral/derived alleles follow the ms convention -> 0: ancestral, 1: derived.
- sampled alleles at a given SNP are represented by strings of 0s and 1s where populations are separated by a slash. E.g., the string 0100/110101 represents the alleles of a SNP for which we have sampled 4 and 6 chromosomes from population *i* and *j*, respectively.


### Statistics of genetic diversity
#### Within population statistics
- __sumpairdif__: sum of pairwise difference? [**to be checked - Martin can you confirm?**]
- __theta__: Watterson *Theta*.
- __D__: Tajima's *D* [(pi-theta)/variance] ? [**to be checked - Martin can you confirm?**]
- __pi__: nucleotide diversity (Tajima's *Theta*).


#### Between population statistics
- __d_*i*x*j*__: Raw nucleotidic divergence between species *i* and *j* (Nei's *Dxy*, eq.12.66, Nei and Kumar 2000).
- __dn_*i*x*j*__: Net nucleotidic divergence between species *i* and *j*  (Nei's *DA*, eq.12.67, Nei and Kumar 2000).
- __FST_*i*x*j*__: (1-(pi_*i*+pi_*j*)/2)/pi_total; also called FST sensu Hudson or Nst.
- __bialsites_*i*x*j*__: mean number of bi-allelic sites per population (infinite site model). I.e. (number of bi-allelic sites in pop*i* + number of bi-allelic sites in pop*j*)/2
- __sfA_*i*x*j*__: number of sites fixed for the derived allele in species A  and fixed for the ancestral allele in species B (1111/0000), where A=*i* and B=*j*.
- __sfB_*i*x*j*__: number of sites fixed for the ancestral allele in species A and fixed for the derived allele in species B (0000/1111), where A=*i* and B=*j*.
- __sxA_*i*x*j*__: number of sites that are polymorphic in species A and fixed for the ancestral allele in species B (0101/0000), where A=*i* and B=*j* .
- __sxB_*i*x*j*__: number of sites that are fixed for the ancestral allele in species A and polymorphic in species B (0000/0101), where A=*i* and B=*j* .
- __sfout_*i*x*j*__: ? [**to be checked**]
- __sxAfB_*i*x*j*__: number of sites polymorphic in species A and fixed for the derived allele in species B (0101/1111), where A=*i* and B=*j*.
- __sxBfA_*i*x*j*__: number of sites fixed for the derived allele in species A and polymorphic in species B (1111/0010), where A=*i* and B=*j*.
- __ss_*i*x*j*__	: number of sites with shared derived alleles between pop*i* and pop*j* (1010/1110). Mean over populations.  [**to be checked - this stats should be symetrical between pops, so that the mean equal the value in each pop actually unless we divide by the total number of pops (that might be different from 2 but that would be surprising.**]
- __Wald_*i*x*j*__	: statistics based on intra locus recombination, see Navascues *et al. BMC Evol. Biol* 2014.


#### Patterson's statistics (requuire more than two populations)
d statistics and f statistics..
NOT IN USED FOR THE MOMENT
TO BE DEFINED.