# Template folder structure for a single sample


```
SAMPLE_ID
|
|_______raw
|         |_raw_molecular_data
|	        |                   |
|         |                   |_*.fastq.gz files, etc.
|         |
|	        |_______raw_image_data
|	                           |
|                            |_*.TIFF files
|                            |_(Whole Slide Images)
|                            |_(etc.)
|_______processed
|		    |_______molecular
|	      |       |
|       |       |_counts.tsv (output of st_pipeline)
|	      | 
|		    |_______spatial
|	              |
|               |_tissue_hires_image.png
|               |_tissue_lowres_image.png
|               |_tissue_positions_list.csv
|               |_scalefactors_json.json
|               |_intersections.csv
|               |_intersections_matx.txt
|               |_masks
|                     |_tissue_hires_intissuemask.png
|
|_______experiment_metadata

```
