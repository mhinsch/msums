# msums

Program for efficient computation of summary statistics from ms outputs.
___

## Description of statistics
### Conventions
- The suffixes 'avg' and 'std' stand for average and standrad deviation respectively.

- Populations are referred by indexes *i* and *j* ranging from 0 to n.

- The coding of ancestral/derived alleles follow ms convention -> 0: ancestral, 1: derived.

- SNPs genotyped at 4 chromomes sampled from each panmictic population are represented as 0100/1101, where the slash delimits populations.

As an example 'FST_*i*x*j*_mean' refers to the mean Fst across loci between populations *i* and *j*.


### Statistics of genetic diversity
#### Within population statistics
- __sumpairdif__: sum of pairwise difference? [to be checked]

- __theta__: Watterson *Theta*.

- __D__: Tajima's *D* [(pi-theta)/variance] ? [to be checked]

- __pi__: nucleotide diversity (Tajima's *Theta*).


#### Between population statistics
- __d_*i*x*j*__: Raw nucleotidic divergence between species *i* and *j* (Nei's *Dxy*, eq. Khumar and Nei, 200x) ?
 [to be checked]

- __dn_*i*x*j*__: Net nucleotidic divergence between species *i* and *j*  (Nei's *DA*, eq. Khumar and Nei, 200x) ?
 [to be checked]

- __FST_*i*x*j*__: (1-(pi_*i*+pi_*j*)/2)/pi_total; also called FST sensu Hudson or Nst.

- __bialsites_*i*x*j*__: number of bi-allelic sites in a sample made of all individuals from pop*i* + pop*j*. [to be checked]

- __sfA_*i*x*j*__: number of sites fixed for the derived allele in species A  and fixed for the ancestral allele in species B (1111/0000), where A=*i* and B=*j*.

- __sfB_*i*x*j*__: number of sites fixed for the ancestral allele in species A and fixed for the derived allele in species B (0000/1111), where A=*i* and B=*j*.

- __sxA_*i*x*j*__: number of polymorphic sites in species A? [to be checked]

- __sxB_*i*x*j*__: number of polymorphic sites in species B? [to be checked]

- __sfout_*i*x*j*__: ? [to be checked]

- __sxAfB_*i*x*j*__: number of sites polymorphic in species A and fixed for the derived allele in species B (0101/1111), where A=*i* and B=*j*.

- __sxBfA_*i*x*j*__: number of sites fixed for the derived allele in species A and polymorphic in species B (1111/0010), where A=*i* and B=*j*.

- __ss_*i*x*j*__	: number of sites with shared derived alleles between pop*i* and pop*j* (1010/1110). Mean over populations.

- __Wald_*i*x*j*__	: statistics based on intra locus recombination, see Navascues *et al.* 201x.


#### Patterson's statistics
TO BE DEFINED.