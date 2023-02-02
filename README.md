# MetaMLST #
A computational pipeline for MLST typing from metagenomic data

MetaMLST performs an *in-silico* Multi Locus Sequence Typing (MLST) Analysis on metagenomic samples, directly from raw reads. MetaMLST achieves cultivation- and assembly- free strain level tracking and profiles all the species to which the standard MLST protocol is applicable.

**Current database version**: [ May 2022 ]

### What can I use MetaMLST for? ###

* Strain-level identification
* Microbial population analysis

### How can I use it? ###

This repo contains sub-repositories, to clone it use:
```
git clone --recurse-submodules https://github.com/SegataLab/metamlst.git
```

Check out the [**Quick Start**](https://github.com/SegataLab/metamlst/wiki#-quick-start) or refer to the [**Wiki**](https://github.com/SegataLab/metamlst/wiki/) for the *documentation*. You can also try the [**examples**](https://zenodo.org/record/4399251/files/metamlst_examples.zip?download=1), or take a look to the [**Examples**](https://github.com/SegataLab/metamlst/wiki/Examples) section

### How it works? ###

MetaMLST reconsctructs the MLST loci-sequences using the closest reference from the publicly available datsets (PubMLST) and traces the most abundant strain of each species:
![MetaMLST pipeline schema](http://segatalab.github.io/images/metamlst_working_concept.jpg)

### Databases ###

MetaMLST automatically downloads the latest version of its database. You can also manually download databases from [here](https://zenodo.org/record/4399251#.X-uTwVn0muU). 
The latest database is: **metamlstDB_2022**.

### Publication ###

* Moreno Zolfo, Adrian Tett, Olivier Jousson, Claudio Donati and Nicola Segata - **[MetaMLST: multi-locus strain-level bacterial typing from metagenomic samples](http://nar.oxfordjournals.org/content/early/2016/09/19/nar.gkw837.full)** - *Nucleic Acids Research, 2016* DOI: 10.1093/nar/gkw837

### Where can I get support? ###

* Bugs and Support: [**MetaMLST Users Support Group**](https://groups.google.com/forum/#!forum/metamlst)
* Project page at [SegataLab](http://segatalab.cibio.unitn.it/tools/metamlst/index.html)
* Project [Wiki](https://github.com/SegataLab/metamlst/wiki/)
