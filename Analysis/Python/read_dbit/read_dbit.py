import os
from warnings import warn
import scanpy as sc
import anndata as ad
import seaborn as sns
import pandas as pd

import matplotlib.pyplot as plt

from pathlib import Path, PurePath
from typing import Union, Dict, Optional, Tuple, BinaryIO

import h5py
import json
import numpy as np

from matplotlib.image import imread
from skimage import io

import anndata
from anndata import (
    AnnData,
    read_csv,
    read_text,
    read_excel,
    read_mtx,
    read_loom,
    read_hdf,
)
from anndata import read as read_h5ad

class read_dbit:

    def __init__(self):
        self.docstring = """
    Read DBiT-seq dataset.
    In addition to reading DBiT-seq count matrices,
    this looks for the `spatial` folder and loads images,
    coordinates and scale factors.
    
    Parameters
    ----------
    path
        Path to directory for DBiT datafiles.
    genome
        Filter expression to genes within this genome.
    count_file
        Which file in the passed directory to use as the count file. Typically would be one of:
        'filtered_feature_bc_matrix.h5' or 'raw_feature_bc_matrix.h5'.
    library_id
        Identifier for the visium library. Can be modified when concatenating multiple adata objects.
    source_image_path
        Path to the high-resolution tissue image. Path will be included in
        `.uns["spatial"][library_id]["metadata"]["source_image_path"]`.
    Returns
    -------
    Annotated data matrix, where observations/cells are named by their
    barcode and variables/genes by gene name. Stores the following information:
    :attr:`~anndata.AnnData.X`
        The data matrix is stored
    :attr:`~anndata.AnnData.obs_names`
        Cell names
    :attr:`~anndata.AnnData.var_names`
        Gene names
    :attr:`~anndata.AnnData.var`\\ `['gene_ids']`
        Gene IDs
    :attr:`~anndata.AnnData.var`\\ `['feature_types']`
        Feature types
    :attr:`~anndata.AnnData.uns`\\ `['spatial']`
        Dict of spaceranger output files with 'library_id' as key
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['images']`
        Dict of images (`'hires'` and `'lowres'`)
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['scalefactors']`
        Scale factors for the spots
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['metadata']`
        Files metadata: 'chemistry_description', 'software_version', 'source_image_path'
    :attr:`~anndata.AnnData.obsm`\\ `['spatial']`
        Spatial spot coordinates, usable as `basis` by :func:`~scanpy.pl.embedding`.

"""

def read_dbit(
        path: Union[str, Path],
        genome: Optional[str] = None,
        *,
        count_file: str = "filtered_feature_bc_matrix.h5",
        library_id: str = None,
        load_images: Optional[bool] = True,
        source_image_path: Optional[Union[str, Path]] = None,
) -> AnnData:
    
    path = Path(path)
    
    p = pd.read_csv(os.path.join(path,count_file), sep='\t', header=0, index_col=0)
    adata = AnnData(X=p, dtype=np.float32)

    adata.uns["spatial"] = dict()

#     from h5py import File
#     with File(path / count_file, mode="r") as f:
#         attrs = dict(f.attrs)
#     if library_id is None:
#    library_id = str(attrs.pop("library_ids")[0], "utf-8")

    attrs = {"chemistry_description": "Spatial 3' v1", "software_version": "dbitranger-1.1.0"}

    adata.uns["spatial"][library_id] = dict()

    if load_images:
        files = dict(
            tissue_positions_file=path / 'spatial/tissue_positions_list.csv',
            scalefactors_json_file=path / 'spatial/scalefactors_json.json',
            hires_image=path / 'spatial/tissue_hires_image.png',
            lowres_image=path / 'spatial/tissue_lowres_image.png',
            tissue_mask=path / 'spatial/masks/tissue_hires_in_tissue_mask.png',
        )

        # check if files exists, continue if images are missing
        for f in files.values():
            if not f.exists():
                if any(x in str(f) for x in ["hires_image", "lowres_image"]):
                    warn(
                        f"You seem to be missing an image file.\n"
                        f"Could not find '{f}'."
                    )
                elif any(x in str(f) for x in ["tissue_positions"]):
                    warn(
                        f"You seem to be missing the tissue position list file:.\n"
                        f"Could not find '{f}'. Don't worry, I'll build one for you."
                    )
                else:
                    raise OSError(f"Could not find required file '{f}'")

        adata.uns["spatial"][library_id]['images'] = dict()
        for res in ['hires', 'lowres']:
            try:
                image_in = imread(str(files[f'{res}_image']))
            except Exception:
                raise OSError(f"Could not find '{res}_image'")

            if (len(np.shape(image_in)) == 2):
                image_in = np.repeat(image_in[:, :, np.newaxis], 3, axis=2)
            adata.uns["spatial"][library_id]['images'][res] = image_in

        # read json scalefactors
        try:
            adata.uns["spatial"][library_id]['scalefactors'] = json.loads(
                files['scalefactors_json_file'].read_bytes()
            )
        except Exception:
            # MMD temporary placeholder until I figure out what scale factors should be used for DBiT-seq
            adata.uns["spatial"][library_id]['scalefactors'] = {'fiducial_diameter_fullres': 288.25050000000005,
                                                                'spot_diameter_fullres': 178.44077999999996,
                                                                'tissue_hires_scalef': 3,
                                                                'tissue_lowres_scalef': 3}

        adata.uns["spatial"][library_id]["metadata"] = {
            k: (str(attrs[k], "utf-8") if isinstance(attrs[k], bytes) else attrs[k])
            for k in ("chemistry_description", "software_version")
            if k in attrs
        }

        ### Read coordinates

        # Check if tissue_positions_file exists, and build one if it doesn't.

        if not os.path.exists(files['tissue_positions_file']):
            tmi = None
            if os.path.exists(files['tissue_mask']):
                # If there is a tissue mask file, use this to populate the 'in_tissue' column of the position list file
                tmi = files['tissue_mask']
            pl = buildpositionlist(count_file=os.path.join(path,count_file), imfile=files['hires_image'], tissuemask_imfile=tmi)
            pl.to_csv(files['tissue_positions_file'], header=False)
            print('No position list file found (should be at \'spatial/tissue_positions_list.csv\'); Building one now.')
            
        positions = pd.read_csv(files['tissue_positions_file'], header=None)
        positions.columns = [
            'barcode',
            'in_tissue',
            'array_row',
            'array_col',
            'pxl_col_in_fullres',
            'pxl_row_in_fullres',
        ]
        positions.index = positions['barcode']

        adata.obs = adata.obs.join(positions, how="left")

        adata.obsm['spatial'] = adata.obs[
            ['pxl_row_in_fullres', 'pxl_col_in_fullres']
        ].to_numpy()
        adata.obs.drop(
            columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'],
            inplace=True,
        )

        # put image path in uns
        if source_image_path is not None:
            # get an absolute path
            source_image_path = str(Path(source_image_path).resolve())
            adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = str(
                source_image_path
            )

    return adata


def buildpositionlist(count_file=None, imfile=None, tissuemask_imfile=None):
    # load the image
    im = imread(imfile)

    # convert image to numpy array
    hires_image = np.asarray(im)

    # load the count file (just to get the rownames (e.g. '5x10', etc.)
    counts = pd.read_csv(count_file, sep='\t', header=0, index_col=0)
    rownames = counts.index

    # split the DBiT pixel names into list of indices, and reindex at 0
    test = [np.short(n.split('x'))-1 for n in rownames]

    # get step sizes for the high resolution image
    if (len(np.shape(hires_image)) == 3):
        w,l, _ = np.shape(hires_image)
    elif (len(np.shape(hires_image)) == 2):
        w,l = np.shape(hires_image)
    w1 = w/99
    w2 = w/198
    l1 = l/99
    l2 = l/198

    # Build and save the table of DBiT pixel positions on the high res image
    position_list = pd.DataFrame(test, index=rownames, columns=['DBiT_pixel_A','DBiT_pixel_B'])
    position_list['x'] = np.short(np.round(position_list['DBiT_pixel_A']*w1*2 + w2))
    position_list['y'] = np.short(np.round(position_list['DBiT_pixel_B']*l1*2 + l2))
    if (tissuemask_imfile is not None):
        mask = io.imread(tissuemask_imfile, as_gray=True)
        position_list['in_tissue'] = [int(mask[coords[0],coords[1]]>0.5) for coords in position_list[['x','y']].to_numpy()]
    else:
        position_list['in_tissue'] = np.short(np.ones(np.shape(rownames)))

    position_list = position_list[['in_tissue','DBiT_pixel_A','DBiT_pixel_B','x','y']]
        
    position_list_out = position_list.sort_values(by=['DBiT_pixel_A','DBiT_pixel_B'], ascending=[True, True])
    
    return position_list_out
