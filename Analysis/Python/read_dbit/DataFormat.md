This is the directory structure for a single DBiT-seq sample for import into python as an AnnData object.

```
Sample_ID
|       count_file.tsv (samples x genes matrix of counts per pixel)
|
|_______spatial
|	|
|	| tissue_positions_list.csv
|	| scalefactors_json.json
|	| tissue_hires_image.png
|       | tissue_lowres_image.png
|       |_______masks
|               |
|               | tissue_hires_in_tissue_mask.png
|               | structure_annotation_01_mask.png
|               | structure_annotation_02_mask.png
|               | structure_annotation_03_mask.png
|
|_______experiment_metadata
|	|
|	| sample_info.txt
|	| missing_channels.txt
|	
|_______raw_image_data
	|
	| 4x_DBiT_TissueProfile_NoCHIP.tif
	| 4x_DBiT_TissueProfile_CHIP1.tif
	| 4x_DBiT_TissueProfile_CHIP2.tif
	| 10x_HandE_WholeSlideImage.svs


```

The main 'omics' data is in ``count_file.tsv``, a tab-delimited file of gene counts for each DBiT 'pixel' which has the following format:

```
	GENE1	GENE2	GENE3	GENE4	GENE5	...
23x28	7	0	3	38	1	...
28x28	16	1	2	82	7	...
37x28	5	0	0	15	1	...
17x28	1	0	1	12	0	...
47x28	3	0	0	1	0	...
5x28	7	0	0	21	4	...
20x28	6	0	1	43	5	...
20x27	1	0	0	2	0	...
28x25	3	0	0	0	1	...
.	.	.	.	.	.
.	.	.	.	.	.
.	.	.	.	.	.
```

The column labels are the pixel locations (``AxB``).

