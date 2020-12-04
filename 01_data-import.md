01\_data-import
================

  - [Importation des données et des
    IP](#importation-des-données-et-des-ip)
      - [Downloading the Standard Operating Procedure example
        data](#downloading-the-standard-operating-procedure-example-data)
      - [decompress the data](#decompress-the-data)
  - [On importe la base de donnée
    silva](#on-importe-la-base-de-donnée-silva)
  - [On importe une autre base de donnée
    silva](#on-importe-une-autre-base-de-donnée-silva)

# Importation des données et des IP

## Downloading the Standard Operating Procedure example data

``` bash
wget https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip -nv
```

    ## 2020-12-04 21:25:14 URL:https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip [36910055/36910055] -> "miseqsopdata.zip.1" [1]

## decompress the data

``` bash
unzip miseqsopdata.zip
```

    ## Archive:  miseqsopdata.zip
    ##    creating: MiSeq_SOP/
    ##   inflating: MiSeq_SOP/F3D0_S188_L001_R1_001.fastq  
    ##    creating: __MACOSX/
    ##    creating: __MACOSX/MiSeq_SOP/
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D0_S188_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D0_S188_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D0_S188_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D141_S207_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D141_S207_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D141_S207_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D141_S207_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D142_S208_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D142_S208_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D142_S208_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D142_S208_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D143_S209_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D143_S209_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D143_S209_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D143_S209_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D144_S210_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D144_S210_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D144_S210_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D144_S210_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D145_S211_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D145_S211_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D145_S211_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D145_S211_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D146_S212_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D146_S212_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D146_S212_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D146_S212_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D147_S213_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D147_S213_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D147_S213_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D147_S213_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D148_S214_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D148_S214_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D148_S214_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D148_S214_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D149_S215_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D149_S215_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D149_S215_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D149_S215_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D150_S216_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D150_S216_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D150_S216_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D150_S216_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D1_S189_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D1_S189_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D1_S189_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D1_S189_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D2_S190_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D2_S190_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D2_S190_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D2_S190_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D3_S191_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D3_S191_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D3_S191_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D3_S191_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D5_S193_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D5_S193_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D5_S193_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D5_S193_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D6_S194_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D6_S194_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D6_S194_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D6_S194_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D7_S195_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D7_S195_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D7_S195_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D7_S195_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D8_S196_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D8_S196_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D8_S196_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D8_S196_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/F3D9_S197_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D9_S197_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/F3D9_S197_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._F3D9_S197_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/HMP_MOCK.v35.fasta  
    ##   inflating: __MACOSX/MiSeq_SOP/._HMP_MOCK.v35.fasta  
    ##   inflating: MiSeq_SOP/Mock_S280_L001_R1_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._Mock_S280_L001_R1_001.fastq  
    ##   inflating: MiSeq_SOP/Mock_S280_L001_R2_001.fastq  
    ##   inflating: __MACOSX/MiSeq_SOP/._Mock_S280_L001_R2_001.fastq  
    ##   inflating: MiSeq_SOP/mouse.dpw.metadata  
    ##   inflating: __MACOSX/MiSeq_SOP/._mouse.dpw.metadata  
    ##   inflating: MiSeq_SOP/mouse.time.design  
    ##   inflating: __MACOSX/MiSeq_SOP/._mouse.time.design  
    ##   inflating: MiSeq_SOP/stability.batch  
    ##   inflating: __MACOSX/MiSeq_SOP/._stability.batch  
    ##   inflating: MiSeq_SOP/stability.files  
    ##   inflating: __MACOSX/MiSeq_SOP/._stability.files

# On importe la base de donnée silva

``` bash
wget https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz -nv
```

    ## 2020-12-04 21:25:17 URL:https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz [81840166/81840166] -> "silva_species_assignment_v138.fa.gz" [1]

# On importe une autre base de donnée silva

``` bash
wget https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz -nv
```

    ## 2020-12-04 21:25:20 URL:https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz [137973851/137973851] -> "silva_nr99_v138_train_set.fa.gz" [1]
