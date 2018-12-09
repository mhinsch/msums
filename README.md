# msums

A program for the efficient computation of a number of population genetics summary statistics. msums can read ms-format data 
on (nearly) arbitrary numbers of populations.

___

## Download
===
You need to have Git installed, then:
```
git clone https://github.com/mhinsch/msums
```

## Installation
===
The installation requires the boost and boost-dev libraries. For 
Manjaro/Archlinux, type:

```bash
sudo pacman -S boost boost-libs 
```


```bash
touch Makefile.dep
make
```

Do not pay attention to warnings like:

```
`stats_multi.h:365:8: warning: unused parameter ‘stop2’ [-Wunused-parameter]
 double patterson_f4
        ^
stats_multi.h:365:8: warning: unused parameter ‘stop3’ [-Wunused-parameter]
stats_multi.h:365:8: warning: unused parameter ‘stop4’ [-Wunused-parameter]
...
stats_multi.h:343:8: warning: unused parameter ‘stop2’ [-Wunused-parameter]
 double patterson_f3
        ^
stats_multi.h:343:8: warning: unused parameter ‘stop3’ [-Wunused-parameter]
...
stats_multi.h:325:8: warning: unused parameter ‘stop2’ [-Wunused-parameter]
 double patterson_D(
        ^
```


## Description of statistics
===
### Conventions
In the following text:
- The suffixes 'mean' and 'std' stand for average and standard deviation, respectively.
- Populations are referred by indexes *i* and *j* and can range from 0 to n, where n is the number of sampled populations.
  - As an example 'FST_*i*x*j*_mean' refers to the mean Fst across loci between populations *i* and *j*.
- The coding of ancestral/derived alleles follow the ms convention -> 0: ancestral, 1: derived.
- sampled alleles at a given SNP are represented by strings of 0s and 1s where populations are separated by a slash. E.g., the string 0100/110101 represents the alleles of a SNP for which we have sampled 4 and 6 chromosomes from population *i* and *j*, respectively.


### Statistics of genetic diversity
#### Within population statistics
- __pairdif__: sum of pairwise allele differences
- __segr__: number of segregating sites(i.e. SNPs) per locus
- __singlet__: overall number of singleton alleles (across all sites)
- __thpi__: Tajima's *Theta*, i.e. nucleotide diversity
- __thW__: Watterson *Theta*.
- __flDstar__: Fu & Li’s _D*_
- __flFstar__: Fu & Li’s _F*_
- __tD__: Tajima's *D* [(pi-theta)/variance] 
- __R2__: Ramos-Onsins *R2* test (Ramos-Onsins & Rozas, *Mol.Biol.Evol.* 2002)

#### Between population statistics
- __d_*i*x*j*__: Raw nucleotidic divergence between species *i* and *j* (Nei's *Dxy*, eq.12.66, Nei and Kumar 2000).
- __dn_*i*x*j*__: Net nucleotidic divergence between species *i* and *j*  (Nei's *DA*, eq.12.67, Nei and Kumar 2000).
- __FST_*i*x*j*__: (1-(pi_*i*+pi_*j*)/2)/pi_total; also called FST sensu Hudson or Nst.
- __bialsites_*i*x*j*__: mean number of bi-allelic sites per population (infinite site model). I.e. (number of bi-allelic sites in pop*i* + number of bi-allelic sites in pop*j*)/2
- __multisites_*i*x*j*__: mean number of multi-allelic (more than 2 alleles) sites per population (infinite site model). I.e. (number of bi-allelic sites in pop*i* + number of bi-allelic sites in pop*j*)/2 [**to be checked**]
- __sfA_*i*x*j*__: number of sites fixed for the derived allele in species A  and fixed for the ancestral allele in species B (1111/0000), where A=*i* and B=*j*.
- __sfB_*i*x*j*__: number of sites fixed for the ancestral allele in species A and fixed for the derived allele in species B (0000/1111), where A=*i* and B=*j*.
- __sfout_*i*x*j*__: ? [**to be checked**]
- __sxA_*i*x*j*__: number of sites that are polymorphic in species A and fixed for the ancestral allele in species B (0101/0000), where A=*i* and B=*j* .
- __sxB_*i*x*j*__: number of sites that are fixed for the ancestral allele in species A and polymorphic in species B (0000/0101), where A=*i* and B=*j* .
- __sxAfB_*i*x*j*__: number of sites polymorphic in species A and fixed for the derived allele in species B (0101/1111), where A=*i* and B=*j*.
- __sxBfA_*i*x*j*__: number of sites fixed for the derived allele in species A and polymorphic in species B (1111/0010), where A=*i* and B=*j*.
- __ss_*i*x*j*__	: number of sites with shared derived alleles between pop*i* and pop*j* (1010/1110). Mean over populations.  [**to be checked - this stat should be symetrical between pops, so that the mean equal the value in each pop actually unless we divide by the total number of pops (that might be different from 2 but that would be surprising.**]
- __Rf_*i*x*j*__: see Navascues *et al. BMC Evol. Biol* 2014. [**to be detailed**]
- __Rs_*i*x*j*__: see Navascues *et al. BMC Evol. Biol* 2014. [**to be detailed**]
- __Wx2s1_*i*x*j*__	: see Navascues *et al. BMC Evol. Biol* 2014. [**to be detailed**]
- __Wx1s2_*i*x*j*__	: see Navascues *et al. BMC Evol. Biol* 2014. [**to be detailed**]
- __pattD_*i*x*j*__: Patterson's *D* statistic used in the "ABBA-BABA" test (Patterson *et al. Genetics* 2012). [How it is implemented, does it not need 4 pops? Is it F2 maybe? See Patterson *et al. Genetics* 2012]

#### Multi-population statistics (more than two populations)
Those are Patterson's test described in Patterson *et al. Genetics* 2012.
- __f3__: is this test really implemented? [**Martin can you confirm?**]
- __f4__: is this test really implemented? [**Martin can you confirm?**]
- __pattD_*i*x*j*__: Patterson's D, see *insert paper here*. [How it is implemented, does it not need 4 pops? Is it F2 maybe? Patterson *et al. Genetics* 2012] Is this test really implemented? [**Martin can you confirm?**]

