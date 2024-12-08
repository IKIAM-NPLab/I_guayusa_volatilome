*in-house* GC-MS library
================
Jefferson Pastuna
2024-12-02

- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#standard-preparation" id="toc-standard-preparation">Standard
  preparation</a>
- <a href="#gc-ms-instrument" id="toc-gc-ms-instrument">GC-MS
  instrument</a>
- <a href="#compound-deconvolution"
  id="toc-compound-deconvolution">Compound deconvolution</a>
- <a href="#library-assembly" id="toc-library-assembly">Library
  assembly</a>
  - <a href="#library-to-ms-dial" id="toc-library-to-ms-dial">Library to
    MS-DIAL</a>
  - <a href="#library-to-nist-ms-search"
    id="toc-library-to-nist-ms-search">Library to NIST MS Search</a>

# Introduction

The document aims to record the spectral Library assembly of authentic
compounds analyzed in the GC-EI (Q)MS. A brief explanation and graphics
are included for each step.

# Standard preparation

All standards (**Table 1**) were accurately weight and diluted in the
mass spectrometry-grade methanol. Working solutions of authentic
compounds were prepared by dilution of the stock solutions (**Table 1**)
with methanol to the final concentration of 1 µg/mL, 5 µg/mL, and 8
µg/mL.

**Table 1.** Authentic compounds description and weight.

| Name                  | Vendor        | Weight ± STD (mg)\* |
|-----------------------|---------------|---------------------|
| Theobromine           | Sigma-Aldrich | 1.3 ± 0.3           |
| Theophylline          | Sigma-Aldrich | 1.3 ± 0.1           |
| Caffeine              | Sigma-Aldrich | 3.9 ± 0.4           |
| Atropine              | Sigma-Aldrich | 2.0 ± 0.1           |
| Harmaline             | Sigma-Aldrich | 2.8 ± 0.1           |
| Linalool              | Sigma-Aldrich | 36.5 ± 0.1          |
| trans-Cinnamaldehyde  | Sigma-Aldrich | 32.2 ± 0.2          |
| Methyl cinnamate      | Sigma-Aldrich | 4.54 ± 0.2          |
| Methyl tetracosanoate | Sigma-Aldrich | 4.14 ± 0.3          |
| Nicotinic acid        | Sigma-Aldrich | 1.99 ± 0.1          |
| Resorcinol            | Sigma-Aldrich | 4.0 ± 0.1           |
| Coumarin              | Sigma-Aldrich | 3.2 ± 0.1           |
| Anthraquinone         | Sigma-Aldrich | 1.2 ± 0.1           |
| Eugenol               | Sigma-Aldrich | 27.1 ± 0.1          |
| n-Alkanes (C7-C33)    | Restek        | na                  |

> Were STD, standard desviation (n = 10), na, 100 uL of n-Alkanes were
> added to 900 μL of n-Hexane, \*, stock solutions was made by adding 1
> mL of MeOH.

# GC-MS instrument

**Column** Rtx-5 GC capillary column, 30 m, 0.25 mm ID, 0.25 µm (cat.#
10223)

**Injection**

Inj. Vol.: 1 µL splitless

Inj. Temp.: 270 °C

**Oven**

Oven Temp.: 70 °C to 200 °C (hold 10 min) at 5 °C/min to 300 °C (hold 15
min) at 6 °C/min

**Carrier Gas** He, constant flow

Flow Rate: 1 mL/min

**Detector** MS

Mode: Scan

Transfer Line Temp.: 280 °C

Analyzer Type: Quadrupole

Source Temp.: 230 °C

Ionization Mode: EI

Scan Range (amu): 50-500

**Instrument** Shimadzu GC-2030 & GCMS-QP2020 NX

# Compound deconvolution

The raw MS data [(.qgd)](Data/in-house_Library/qgd_raw_Files) were
converted to [(.abf)](Data/in-house_Library/abf_raw_Files) format using
Reifycs file converter. MS-DIAL 4.9.221218 was used to deconvolute the
raw spectral data using a pipeline that included peak detection,
alignment, filtering, and gap filling. Parameters used in MS-DIAL are
available in [(.txt)](Data/in-house_Library/MS-DIAL_Parameter) and
[(.med2)](Data/in-house_Library/MS-DIAL_Parameter) format.

# Library assembly

Each deconvoluted compounds was exported to MS-FINDER ver. 3.60 as see
in below imagen.

**Table 1.** Authentic compounds RT and RI analized in GC-MS.

| Name                  | Formula                                               | RT ± STD (min) | RI ± STD |
|-----------------------|-------------------------------------------------------|----------------|----------|
| Theobromine           | C<sub>7</sub>H<sub>8</sub>N<sub>4</sub>O<sub>2</sub>  | 25.500 ± 0.044 | 1877 ± 3 |
| Theophylline          | C<sub>7</sub>H<sub>8</sub>N<sub>4</sub>O<sub>2</sub>  | 28.193 ± 0.042 | 1999 ± 2 |
| Caffeine              | C<sub>8</sub>H<sub>10</sub>N<sub>4</sub>O<sub>2</sub> | 25.112 ± 0.006 | 1858 ± 0 |
| Atropine              | C<sub>17</sub>H<sub>23</sub>NO<sub>3</sub>            | 36.645 ± 0.013 | 2236 ± 0 |
| Harmaline             | C<sub>13</sub>H<sub>14</sub>N<sub>2</sub>O            | 36.303 ± 0.016 | 2228 ± 0 |
| Linalool              | C<sub>10</sub>H<sub>18</sub>O                         | 7.308 ± 0.006  | 1103 ± 0 |
| trans-Cinnamaldehyde  | C<sub>9</sub>H<sub>8</sub>O                           | 11.643 ± 0.003 | 1277 ± 0 |
| Methyl cinnamate      | C<sub>10</sub>H<sub>10</sub>O<sub>2</sub>             | 14.485 ± 0.005 | 1390 ± 0 |
| Methyl tetracosanoate | C<sub>25</sub>H<sub>50</sub>O<sub>2</sub>             | 47.878 ± 0.003 | 2735 ± 0 |
| Nicotinic acid        | C<sub>6</sub>H<sub>5</sub>NO<sub>2</sub>              | 10.405 ± 0.000 | 1229 ± 0 |
| Resorcinol            | C<sub>6</sub>H<sub>6</sub>O<sub>2</sub>               | 11.815 ± 0.010 | 1284 ± 0 |
| Coumarin              | C<sub>9</sub>H<sub>6</sub>O<sub>2</sub>               | 15.908 ± 0.003 | 1447 ± 0 |
| Anthraquinone         | C<sub>14</sub>H<sub>8</sub>O<sub>2</sub>              | 27.798 ± 0.006 | 1982 ± 1 |
| Eugenol               | C<sub>10</sub>H<sub>12</sub>O<sub>2</sub>             | 13.833 ± 0.003 | 1364 ± 0 |
| Nonane                | C<sub>9</sub>H<sub>20</sub>                           | 3.469 ± 0.004  | 900 ± 0  |
| Decane                | C<sub>10</sub>H<sub>22</sub>                          | 5.098 ± 0.004  | 1000 ± 0 |
| Undecane              | C<sub>11</sub>H<sub>24</sub>                          | 7.232 ± 0.005  | 1100 ± 0 |
| Dodecane              | C<sub>12</sub>H<sub>26</sub>                          | 9.677 ± 0.005  | 1200 ± 0 |
| Tridecane             | C<sub>13</sub>H<sub>28</sub>                          | 12.224 ± 0.008 | 1300 ± 0 |
| Tetradecane           | C<sub>14</sub>H<sub>30</sub>                          | 14.752 ± 0.005 | 1400 ± 0 |
| Pentadecane           | C<sub>15</sub>H<sub>32</sub>                          | 17.209 ± 0.006 | 1500 ± 0 |
| Hexadecane            | C<sub>16</sub>H<sub>34</sub>                          | 19.550 ± 0.006 | 1600 ± 0 |
| Heptadecane           | C<sub>17</sub>H<sub>36</sub>                          | 21.791 ± 0.006 | 1700 ± 0 |
| Octadecane            | C<sub>18</sub>H<sub>38</sub>                          | 23.927 ± 0.005 | 1800 ± 0 |
| Nonadecane            | C<sub>19</sub>H<sub>40</sub>                          | 25.966 ± 0.006 | 1900 ± 0 |
| Eicosane              | C<sub>20</sub>H<sub>42</sub>                          | 28.220 ± 0.008 | 2000 ± 0 |
| Heneicosane           | C<sub>21</sub>H<sub>44</sub>                          | 31.173 ± 0.008 | 2100 ± 0 |
| Docosane              | C<sub>22</sub>H<sub>46</sub>                          | 35.188 ± 0.015 | 2200 ± 1 |
| Tricosane             | C<sub>23</sub>H<sub>48</sub>                          | 39.210 ± 0.011 | 2300 ± 0 |
| Tetracosane           | C<sub>24</sub>H<sub>50</sub>                          | 41.914 ± 0.008 | 2400 ± 0 |
| Pentacosane           | C<sub>25</sub>H<sub>52</sub>                          | 44.032 ± 0.008 | 2500 ± 1 |
| Hexacosane            | C<sub>26</sub>H<sub>54</sub>                          | 45.816 ± 0.006 | 2600 ± 0 |
| Heptacosane           | C<sub>27</sub>H<sub>56</sub>                          | 47.384 ± 0.006 | 2700 ± 1 |
| Octacosane            | C<sub>28</sub>H<sub>58</sub>                          | 48.801 ± 0.006 | 2800 ± 1 |
| Nonacosane            | C<sub>29</sub>H<sub>60</sub>                          | 50.102 ± 0.007 | 2900 ± 1 |
| Triacontane           | C<sub>30</sub>H<sub>62</sub>                          | 51.328 ± 0.006 | 3000 ± 1 |
| Hentriacontane        | C<sub>31</sub>H<sub>64</sub>                          | 52.479 ± 0.005 | 3100 ± 1 |
| Dotriacontane         | C<sub>32</sub>H<sub>66</sub>                          | 53.638 ± 0.006 | 3200 ± 1 |
| Tritriacontane        | C<sub>33</sub>H<sub>68</sub>                          | 54.952 ± 0.009 | 3299 ± 1 |

> Were RT, retention time in minutes, RI, Semi-standard non-polar Van
> den Dool and Kratz linear retention index, STD, standard desviation (n
> ≥ 3).

## Library to MS-DIAL

This step is used to recover missing compounds in the samples, the
parameter only requires the number of minimum values for which a
compound wants to be ’re-searched’ in the samples.

Exporting alignment feature list for blank subtraction with
[“notame”](https://doi.org/10.3390/metabo10040135) R package. The
procedures used are available in the “Multivariate_Statistics” notebook
of GitHub repository
(<https://github.com/IKIAM-NPLab/I_guayusa_volatilome>). Only
high-quality features were used for metabolite identification.

## Library to NIST MS Search

Metabolite identification was by comparing all the spectra found against
a reference database. The eRah default database was used to metabolite
identification. Metabolite identification was improve by exporting all
the spectra found (.msp) in eRah to NIST MS Search 2.4 software.

Exporting spectra to NIST MS Search software identification with NIST-20
library.
