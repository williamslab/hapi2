HAPI version 2
==============

HAPI (HAPlotype Inference) 2 performs phasing (i.e., haplotype inference) on nuclear families (parents and children) and applies to families with data for both, one, or neither parent (i.e., a set of siblings only). Note that in the latter case, HAPI2 itself does not detect which parent is the same across chromosomes and, for certain ambiguous haplotype transmission patterns, does not determine linkage of distinct segments on either side of the ambiguity. Providing HAPI2's output to [HAPI-RECAP](https://github.com/23andMe/hapi-recap) along with IBD segments to the siblings' relatives determines a set of parent DNA likely to belong to the two parents (rather than being mixed between them, as in HAPI2's output for sibling-only families). HAPI-RECAP can also infer the sexes of the parents using sex-specific genetic maps and the crossovers HAPI2 detects. Technical details on HAPI2 and HAPI-RECAP are available in [the HAPI-RECAP preprint](https://doi.org/10.1101/2024.05.10.593578).

Compiling HAPI2
---------------

HAPI2 requires zlib, including developmental headers (e.g., the `zlib1g-dev` package on Ubuntu).

Clone the repository by running

    git clone --recurse-submodules https://github.com/williamslab/hapi2.git

(Alternatively, `git clone [repo]` followed by `git submodule update --init` in the cloned directory will do the same as the above.)

Next, compile genetio and then HAPI2 by running

    cd hapi2/genetio
    make
    cd ..
    make

To pull HAPI2 updates, use

    git pull
    git submodule update --remote

Running HAPI2
-------------

HAPI2 takes [PLINK binary format data](https://www.cog-genomics.org/plink/1.9/formats#bed) (`.bed`, `.bim`, and `.fam` format files) as input. (Running [PLINK](https://www.cog-genomics.org/plink2/) on the input data with the `--make-bed` option converts to this format.)

HAPI2 reads the family relationships from the `.fam` file. To specify a nuclear family, this file must include parent ids for each child even if the dataset does not contain genotypes for one or both parents. To phase a set of siblings only, use the same father and mother id for each child in the `.fam` file.

The method phases each family independently, so a person that is a child in one family and a parent in another may have different phasing results across the two families.

To run HAPI2 on PLINK format files `data.bed`, `data.bim`, and `data.fam` and print the phase in VCF format, execute

    ./hapi2  -g data.bed -s data.bim -i data.bim  -o data-phase  --vcf

or, equivilently,

    ./hapi2 -p data -o data-phase --vcf

This prints one VCF file for each family into the directory `data-phase`. HAPI2 will create that directory if necessary.

HAPI2 output
------------

HAPI2 writes its phasing results to a directory given by `-o [output_prefix]` option. It also writes a log file named `[output_prefix].log`. It requires one or more output type to run, as specified with the following command-line options:

* `--vcf`, print phased haplotypes in [VCF format](https://samtools.github.io/hts-specs/VCFv4.2.pdf). Prints one file per family.
* `--json`, print phased haplotypes and accompanying metadata (SNP positions, etc.) in JSON format, which HAPI-RECAP reads. Prints data for _all_ families in one file.
* `--json_par`, as above, but only includes data for the parents. Prints data for _all_ families in one file.
* `--ped`, print phased haplotypes in [PLINK ped format](https://www.cog-genomics.org/plink/1.9/formats#ped). This is an individual-major format (each line corresponds to one person) and it lists the two alleles for each site successively. To encode phase, the first allele for each site is from the first haplotype and the second allele is from the second haplotype. Prints one file per family.
* `--txt`, print phased haplotypes in a text format that is intended to be human-readable. Informative markers include the inheritance vectors HAPI2 inferred (as described below). Prints one file for each family and each chromosome.
* `--iv`, print inheritance vectors (see below). One file for each family and each chromosome.
* `--detect_co`, print the locations of crossovers. One file for each family and each chromosome.

HAPI2 infers inheritance vectors that encode which haplotype the two parents transmitted to each child. While the [paper](https://doi.org/10.1101/2024.05.10.593578) describing HAPI2 uses numbers 0 and 1 (and HAPI2 uses numbers interally), HAPI2 prints these haplotypes using the letters `A` and `B`. This is so that ambiguous values can be indicated. Specifically, if the inheritance vector element could be take on either value, the output uses a lowercase `a` or `b`. The printed letter corresponds to the state path HAPI2 arbitrarily chose, but any lowercase element indicates that there is one or more other state paths the have a different inheritance vector value.

Other HAPI2 command-line options
--------------------------------

The following command-line options control which input data HAPI2 analyzes:

* `-c` or `--chr`, only analyze the given chromosome.
* `--start`, analyze physical positions greater or equal to the given integer.
* `--end`, analyze physical positions less than or equal to the given integer.
* `--min_par`, only analyze families with at least the given number of parents. Can be 2, 1, or 0. Default: 0.
* `--min_child`, only analyze families with at least the given number of children. Default: 2.
* `--set_miss_par_bits`, treat one or both parents as missing even if the input contains genotypes for them. Can be 1 (set father as missing), 2 (set mother as missing), 3 (set both as missing), or 0 (use all input data). Default: 0.

The following options control HAPI2's error model:

* `--err_min`, minimum number of recombinations that must occur within a span of `--err_dist` markers to flag those markers as errors. For analyzing sibling-only data, our tests used a value of 1. Default: 2. 
* `--err_dist`, how many sequential informative markers should be considers as part of a single error. The penalty for the error increases with the number of markers. Default: 5.

Crossovers at the beginning of the chromosome sometimes occur after a small number of informative markers. The `--detect_co` option calls crossovers only after a given number of informative markers occur on a single haplotype followed by that many more on a different haplotype. To reduce the number required to support the first called haplotype on a chromosome, use `--edge_co`.

Options controlling HAPI2's output:

* `--force`, overwrite a file if it exists. By default HAPI2 prints an error and does not overwrite existing files.
* `--no_family_id`, PLINK's `.fam` format includes a family id, and by default HAPI2's output encodes ids as a string `[family_id]:[individual_id]`. Using this option removes the family id and `:` from each person's id.
* `--verbose`, print information about each site to the log during processing. Can be useful for debugging.

Citing HAPI2
------------

If you use, HAPI2, please cite our preprint:

Qiao Y, Jewett EM, McManus KF, Freeman WA, Curran JE, Williams-Blangero S, Blangero J, The 23andMe Research Team, Williams AL. Reconstructing parent genomes using siblings and other relatives. bioRxiv. [https://doi.org/10.1101/2024.05.10.593578](https://doi.org/10.1101/2024.05.10.593578).

