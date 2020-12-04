metabarcoding with dada2: environment installation
================

  - [update VM configuration](#update-vm-configuration)
  - [package install](#package-install)
      - [Pour l’arbre phylogénétique -
        phyloseq](#pour-larbre-phylogénétique---phyloseq)
  - [rmarkdown](#rmarkdown)
  - [knitr](#knitr)

# update VM configuration

``` bash
sudo apt-get update -y 
sudo apt-get install -y libbz2-dev
sudo apt-get install -y liblzma-dev
sudo apt-get install -y libglpk-dev
```

    ## sudo: unable to resolve host bbffce91a3f0: Name or service not known
    ## Hit:1 http://archive.ubuntu.com/ubuntu focal InRelease
    ## Get:2 http://archive.ubuntu.com/ubuntu focal-updates InRelease [114 kB]
    ## Get:3 http://archive.ubuntu.com/ubuntu focal-backports InRelease [101 kB]
    ## Get:4 http://security.ubuntu.com/ubuntu focal-security InRelease [109 kB]
    ## Get:5 http://archive.ubuntu.com/ubuntu focal-updates/restricted amd64 Packages [136 kB]
    ## Get:6 http://archive.ubuntu.com/ubuntu focal-updates/universe amd64 Packages [874 kB]
    ## Get:7 http://archive.ubuntu.com/ubuntu focal-updates/main amd64 Packages [876 kB]
    ## Get:8 http://archive.ubuntu.com/ubuntu focal-backports/universe amd64 Packages [4,248 B]
    ## Get:9 http://security.ubuntu.com/ubuntu focal-security/universe amd64 Packages [643 kB]
    ## Get:10 http://security.ubuntu.com/ubuntu focal-security/restricted amd64 Packages [103 kB]
    ## Get:11 http://security.ubuntu.com/ubuntu focal-security/main amd64 Packages [490 kB]
    ## Fetched 3,449 kB in 1s (2,423 kB/s)
    ## Reading package lists...
    ## sudo: unable to resolve host bbffce91a3f0: Name or service not known
    ## Reading package lists...
    ## Building dependency tree...
    ## Reading state information...
    ## libbz2-dev is already the newest version (1.0.8-2).
    ## 0 upgraded, 0 newly installed, 0 to remove and 39 not upgraded.
    ## sudo: unable to resolve host bbffce91a3f0: Name or service not known
    ## Reading package lists...
    ## Building dependency tree...
    ## Reading state information...
    ## liblzma-dev is already the newest version (5.2.4-1ubuntu1).
    ## 0 upgraded, 0 newly installed, 0 to remove and 39 not upgraded.
    ## sudo: unable to resolve host bbffce91a3f0: Name or service not known
    ## Reading package lists...
    ## Building dependency tree...
    ## Reading state information...
    ## libglpk-dev is already the newest version (4.65-2).
    ## 0 upgraded, 0 newly installed, 0 to remove and 39 not upgraded.

# package install

Following instruction on
<https://benjjneb.github.io/dada2/dada-installation.html>

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = '3.11')
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

``` r
BiocManager::install("dada2", version = "3.11")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'dada2'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

``` r
BiocManager::install("phangorn")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'phangorn'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

``` r
BiocManager::install("DECIPHER")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'DECIPHER'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

``` r
BiocManager::install("DESseq2")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'DESseq2'

    ## Warning: package 'DESseq2' is not available for this version of R
    ## 
    ## A version of this package for your version of R might be available elsewhere,
    ## see the ideas at
    ## https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

``` r
BiocManager::install("ggplot2", version = "3.11")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'ggplot2'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

``` r
BiocManager::install("phyloseq", version = "3.11")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'phyloseq'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

## Pour l’arbre phylogénétique - phyloseq

``` r
install.packages("gridExtra")
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

``` r
.cran_packages <- c( "shiny","miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest", "dplyr", "ggrepel", "nlme", "devtools",
                  "reshape2", "PMA", "structSSI", "ade4",
                  "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("genefilter", "impute")
```

``` r
install.packages(.cran_packages)
```

    ## Installing packages into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

    ## Warning: package 'structSSI' is not available for this version of R
    ## 
    ## A version of this package for your version of R might be available elsewhere,
    ## see the ideas at
    ## https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages

``` r
devtools::install_github(.github_packages)
```

    ## Skipping install of 'phyloseqGraphTest' from a github remote, the SHA1 (3fb6c274) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
BiocManager::install(.bioc_packages)
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'genefilter', 'impute'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

``` bash
wget https://cran.r-project.org/src/contrib/Archive/structSSI/structSSI_1.1.1.tar.gz
```

    ## --2020-12-04 18:30:13--  https://cran.r-project.org/src/contrib/Archive/structSSI/structSSI_1.1.1.tar.gz
    ## Resolving cran.r-project.org (cran.r-project.org)... 137.208.57.37
    ## Connecting to cran.r-project.org (cran.r-project.org)|137.208.57.37|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 25591 (25K) [application/x-gzip]
    ## Saving to: ‘structSSI_1.1.1.tar.gz.1’
    ## 
    ##      0K .......... .......... ....                            100% 1.11M=0.02s
    ## 
    ## 2020-12-04 18:30:13 (1.11 MB/s) - ‘structSSI_1.1.1.tar.gz.1’ saved [25591/25591]

``` r
library(devtools)
```

    ## Loading required package: usethis

``` r
install_local("./structSSI_1.1.1.tar.gz")
```

    ## Skipping 1 packages not available: multtest

    ##      checking for file ‘/tmp/RtmpYoUpDF/remotes638571a9cfb/structSSI/DESCRIPTION’ ...  ✓  checking for file ‘/tmp/RtmpYoUpDF/remotes638571a9cfb/structSSI/DESCRIPTION’
    ##   ─  preparing ‘structSSI’:
    ##    checking DESCRIPTION meta-information ...  ✓  checking DESCRIPTION meta-information
    ##   ─  checking for LF line-endings in source and make files and shell scripts
    ##   ─  checking for empty or unneeded directories
    ## ─  looking to see if a ‘data/datalist’ file should be added
    ## ─  building ‘structSSI_1.1.1.tar.gz’
    ##      
    ## 

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

# rmarkdown

``` r
install.packages("rmarkdown")
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

# knitr

``` r
install.packages("knitr")
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)
