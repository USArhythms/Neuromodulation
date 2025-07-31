[Home](../README.md) | [Data Download](Data.md)

# Downloading and Loading the Data

The published dataset for 
[Rauscher et al. 2024. "Neurovascular Impulse Response Function (IRF) during spontaneous activity differentially reflects intrinsic neuromodulation across cortical regions". bioRXiv.](https://www.biorxiv.org/content/10.1101/2024.09.14.612514v1.full)
can be found at https://dandiarchive.org/dandiset/001543


## How to Download

Instructions on how to download data from dandiarchives can be found at
https://docs.dandiarchive.org/user-guide-using/accessing-data/downloading/

Note: not recommended to download entire dandiset! 6.2 TB total

To download a single session:
<pre>
dandi download https://dandiarchive.org/aip/dandisets/001543/versions/draft/assets/cdf5467c-22be-4257-b73d-05eccd54b644/
</pre>

To download other sessions, replace **cdf5467c-22be-4257-b73d-05eccd54b644** with session tag


## How to Load NWB Files in MATLAB

Use matNWB's function nwbRead to read nwb file:
<pre>
nwb = nwbRead(path_to_NWB)
</pre>

Use the included *f_extractNWB* function to extract data:
<pre>
f_extractNWB(nwb)
</pre>

Outputs:
rfp - raw &Delta;F/F Ca++
