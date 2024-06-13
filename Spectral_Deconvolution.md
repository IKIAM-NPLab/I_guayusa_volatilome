Spectral_Deconvolution
================
Edison Gonzales , Jefferson Pastuna
2024-06-09

- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#before-to-start" id="toc-before-to-start">Before to start</a>
- <a href="#erah-package-workflow" id="toc-erah-package-workflow">eRah
  package workflow</a>

### Introduction
En los últimos años se han venido desarrollado tecnologías basadas en espectrometría de masas que han permitido dar un auge a la metabolómica, un ejemplo claro es la cromatografía gaseosa acoplada a espectrometría de masas GC-MS, técnica que permite el análisis de metabos de sustancias volátiles y semi volátiles. La base dicha técnica es la ionización a partir de 70v. 
Este proceso de ionización en GC es altamente reproducible sin embargo la coelución de muestras complejas junto a la extensa fragmentación de iones da como resultado un conjunto de datos grandes y complejos, unido también a la falta de herramientas computacionales integradas a la metabolómica no dirigida en GC-MS dificulta en gran medida en la identificación y cuantificación de metabolitos. 
El análisis de datos de GC-MS se divide en dos categorías principales, la primera es “peak-picking”, que implica la desconvolución para fragmentar e identificar los picos más relevantes en los espectros entre diferentes muestras, permitiéndonos descubrir variaciones estadísticas en los picos entre grupos experimentales. Las herramientas que permiten este tipo de procesamiento de datos son MZmine,MetAligny XCMS. Sin embargo, estos métodos no se basan en los espectros de los compuestos, sino en el valor m/z el tiempo de retención y el área de los picos de iones fragmentados. Lo que ocasiona dificultad en la identificación de los compuestos. 
La segunda categoría se basa en la cuantificación e identificación mediante procesos de descovolucion multivariable que extrae y construye compuestos puros a partir de datos brutos. Las herramientas aplicadas para este proceso son TNO-DECO o ADAP-GC. De igual forma existen sofwares libres como AMDIS or BinBase cada uno de estos sofwares representan limitaciones especificas que impide un análisis amplio de la data obtenida por medio del GC-MS.
Por ende, a base de las limitaciones y parámetros a seguir, se han desarrollado softwares que permiten un análisis amplio y preciso de los datos. Entre ellos destaca eRah, un software gratuito y de código abierto diseñado para procesar datos en metabolómica no dirigida basada en GC-MS. Este paquete de R se basa en la deconvolución central, centrándose principalmente en la separación ciega de fuentes (BSS), la cuantificación y la identificación automatizada de espectros de muestra mediante la comparación con bibliotecas espectrales.

### Before to start
eRah, es un paquete de R gratuito el cual incorpora un método de deconvulcion central, de modo que utiliza técnicas multivariadas las cuales se basan en la separación ciega de fuentes  (BSS) el cual es un proceso que permite la alineación, cuantificación y la identificación de metabolitos a través de la comparación de bibliotecas espectrales, lo que a su vez, permite la obtención de una tabla con los nombres de los compuestos, las puntuaciones de coincidencia y el área integrada del compuesto para cada muestra.

### eRah package workflow

How to star.

``` r
# eRah package installation
#install.packages('erah')
# eRah library call
library(erah)
```

Folder files

``` r
# Delete all file that are not in folders
unlink('Data/Data_to_eRah/*')
# Data folder path
createdt('Data/Data_to_eRah/')
```

New experiment.

``` r
instrumental <- read.csv('Data/Metadata_to_eRah/Metadata_inst.csv')
phenotype <- read.csv('Data/Metadata_to_eRah/Metadata_pheno.csv')

raw_data <- newExp(instrumental = instrumental,
                   phenotype = phenotype,
                   info = 'I. guayusa volatilome')
```

Compound deconvolution

Parameter

``` r
dec_par <- setDecPar(min.peak.width = 3,
                     min.peak.height = 500,
                     noise.threshold = 50,
                     avoid.processing.mz = c(50:69,73:75,147:149))
```

parallel processing

``` r
plan(future::multisession,
     workers = 14)
```

Deconvolution

``` r
dec_data <- deconvolveComp(raw_data,
                           dec_par)
```

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Process_Blank/4_Blank.mzXML ... Processing 1 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/109_QC.mzXML ... Processing 2 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/110_QC.mzXML ... Processing 3 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/23_QC.mzXML ... Processing 4 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/29_QC.mzXML ... Processing 5 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/35_QC.mzXML ... Processing 6 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/50_QC.mzXML ... Processing 7 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/56_QC.mzXML ... Processing 8 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/62_QC.mzXML ... Processing 9 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/7_QC.mzXML ... Processing 10 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/77_QC.mzXML ... Processing 11 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/8_QC.mzXML ... Processing 12 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/83_QC.mzXML ... Processing 13 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/90_QC.mzXML ... Processing 14 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/96_QC.mzXML ... Processing 15 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/10_GC_Neg.mzXML ... Processing 16 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/100_MB_CPos2.mzXML ... Processing 17 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/101_GC_Neg.mzXML ... Processing 18 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/102_GC_0.mzXML ... Processing 19 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/103_GC_Pos.mzXML ... Processing 20 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/104_GC_2.mzXML ... Processing 21 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/105_GC_A.mzXML ... Processing 22 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/106_GC_1.mzXML ... Processing 23 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/107_GC_C.mzXML ... Processing 24 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/108_GC_B.mzXML ... Processing 25 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/11_GC_0.mzXML ... Processing 26 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/12_GC_Pos.mzXML ... Processing 27 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/13_GC_2.mzXML ... Processing 28 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/14_GC_A.mzXML ... Processing 29 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/15_GC_1.mzXML ... Processing 30 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/16_GC_C.mzXML ... Processing 31 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/17_GC_B.mzXML ... Processing 32 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/18_MB_CNeg0.mzXML ... Processing 33 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/19_MB_CNeg2.mzXML ... Processing 34 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/20_MB_APos2.mzXML ... Processing 35 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/21_MB_BPos1.mzXML ... Processing 36 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/22_MB_APos1.mzXML ... Processing 37 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/24_MB_ANeg0.mzXML ... Processing 38 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/25_MB_CPos2.mzXML ... Processing 39 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/26_MB_BPos0.mzXML ... Processing 40 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/27_MB_CNeg0.mzXML ... Processing 41 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/28_MB_APos2.mzXML ... Processing 42 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/30_MB_ANeg1.mzXML ... Processing 43 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/31_MB_APos2.mzXML ... Processing 44 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/32_MB_APos0.mzXML ... Processing 45 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/33_MB_CNeg0.mzXML ... Processing 46 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/34_MB_BPos0.mzXML ... Processing 47 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/37_GC_Neg.mzXML ... Processing 48 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/38_GC_0.mzXML ... Processing 49 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/39_GC_Pos.mzXML ... Processing 50 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/40_GC_2.mzXML ... Processing 51 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/41_GC_A.mzXML ... Processing 52 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/42_GC_1.mzXML ... Processing 53 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/43_GC_C.mzXML ... Processing 54 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/44_GC_B.mzXML ... Processing 55 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/45_MB_APos0.mzXML ... Processing 56 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/46_MB_BPos2.mzXML ... Processing 57 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/47_MB_ANeg1.mzXML ... Processing 58 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/48_MB_CNeg2.mzXML ... Processing 59 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/49_MB_CPos0.mzXML ... Processing 60 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/51_MB_CNeg2.mzXML ... Processing 61 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/52_MB_CNeg1.mzXML ... Processing 62 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/53_MB_BPos1.mzXML ... Processing 63 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/54_MB_APos1.mzXML ... Processing 64 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/55_MB_BNeg0.mzXML ... Processing 65 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/57_MB_APos1.mzXML ... Processing 66 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/58_MB_ANeg2.mzXML ... Processing 67 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/59_MB_BPos2.mzXML ... Processing 68 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/60_MB_CPos2.mzXML ... Processing 69 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/61_MB_BNeg1.mzXML ... Processing 70 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/64_GC_Neg.mzXML ... Processing 71 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/65_GC_0.mzXML ... Processing 72 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/66_GC_Pos.mzXML ... Processing 73 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/67_GC_2.mzXML ... Processing 74 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/68_GC_A.mzXML ... Processing 75 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/69_GC_1.mzXML ... Processing 76 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/70_GC_C.mzXML ... Processing 77 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/71_GC_B.mzXML ... Processing 78 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/72_MB_ANeg2.mzXML ... Processing 79 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/73_MB_BPos2.mzXML ... Processing 80 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/74_MB_CPos0.mzXML ... Processing 81 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/75_MB_BNeg1.mzXML ... Processing 82 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/76_MB_ANeg1.mzXML ... Processing 83 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/78_MB_APos0.mzXML ... Processing 84 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/79_MB_BNeg1.mzXML ... Processing 85 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/80_MB_CPos1.mzXML ... Processing 86 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/81_MB_BPos0.mzXML ... Processing 87 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/82_MB_ANeg0.mzXML ... Processing 88 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/85_MB_BPos1.mzXML ... Processing 89 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/86_MB_CNeg1.mzXML ... Processing 90 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/87_MB_BNeg2.mzXML ... Processing 91 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/88_MB_CPos1.mzXML ... Processing 92 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/89_MB_BNeg0.mzXML ... Processing 93 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/91_MB_CNeg1.mzXML ... Processing 94 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/92_MB_CPos0.mzXML ... Processing 95 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/93_MB_BNeg0.mzXML ... Processing 96 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/94_MB_ANeg0.mzXML ... Processing 97 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/95_MB_BNeg2.mzXML ... Processing 98 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/97_MB_ANeg2.mzXML ... Processing 99 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/98_MB_CPos1.mzXML ... Processing 100 / 101

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/99_MB_BNeg2.mzXML ... Processing 101 / 101

    ## Compounds deconvolved

Alignment

``` r
# Alignment parameters
alig_par <- setAlPar(min.spectra.cor = 0.9,
                     max.time.dist = 3,
                     mz.range = 50:500)
# Alignment
peak_alig <- alignComp(dec_data,
                       alParameters = alig_par)
```

Missing compound recovery

``` r
peak_find <- recMissComp(peak_alig,
                         min.samples = 3)
```

    ## 
    ##  Updating alignment table... 
    ## Model fitted!

Identification

``` r
# Loading NIST 20 (*.msp) library
nist.database <- importMSP(filename = "E:/NIST_20_Library/Result/NIST20EI_2eRah.MSP",
                           DB.name = "NIST",
                           DB.version = "NIST20",
                           DB.info = "NIST MS Search Export")
# Save library for a posterior faster loading
save(nist.database, file= "E:/NIST_20_Library/Result/NIST20EI_2eRah.rda")
# Load R library
load("E:/NIST_20_Library/Result/NIST20EI_2eRah.rda")
mslib <- nist.database
# Identification
peak_iden <- identifyComp(peak_find,
                          id.database = mslib,
                          mz.range = NULL,
                          n.putative = 20)
```

    ## Constructing matrix database... 
    ## Comparing spectra... 
    ## Done!

``` r
# Identified compound
id_list <- idList(peak_iden)
# Exporting identified compounds
write.csv(id_list,
          file = "Result/eRah_Result/NIST_Identification.csv")
```

Mirror plot of identified compounds

``` r
plotSpectra(peak_iden, 1, 1, draw.color = "red")
```

![](Spectral_Deconvolution_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Exporting spectra to NIST identification

``` r
export2MSP(peak_iden,
           store.path = "Result/eRah_Result",
           alg.version = 2)
```

    ## Spectra saved at: Result/eRah_Result/ExportMSP
