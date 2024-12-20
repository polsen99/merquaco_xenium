import os
from typing import Union
from pathlib import Path
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import merquaco.pixel_classification as pc
import merquaco.data_processing as data_processing
from merquaco.data_loss import FOVDropout, DropoutResult
import merquaco.figures as figures
import merquaco.z_plane_detection as zp
import merquaco.perfusion as perfusion
import merquaco.periodicity as periodicity
from merquaco.__init__ import __version__ as version


# Dictionary keys for output metadata
metrics_dict_keys = {'filtered_transcripts_count', 'transcript_density_um2', 'transcript_density_um2_per_gene',
                     'on_tissue_filtered_transcript_count', 'periodicity', 'counts_per_gene', 'n_dropped_fovs',
                     'n_dropped_genes', 'dropped_fovs_dict', 'dropped_genes', 'damage_area', 'transcripts_area',
                     'detachment_area', 'ventricle_area', 'total_area', 'damage_percent', 'transcripts_percent',
                     'detachment_percent', 'ventricle_percent', 'transcripts_mask_pixel_path',
                     'transcripts_mask_object_path', 'dapi_mask_pixel_path', 'dapi_mask_object_path',
                     'ventricle_mask_pixel_path', 'ventricle_mask_object_path', 'merquaco_version'}


class XeniumExperiment:

    def __init__(self,
                 transcripts_input: Union[pd.DataFrame, str, Path],
                 ilastik_program_path: Union[str, Path] = None,
                 dapi_high_res_image_path: Union[str, Path] = None,
                 output_dir: Union[str, Path] = None,
                 ventricle_genes_list: list = ["Crb2", "Glis3", "Inhbb", "Naaa", "Cd24a",
                                               "Dsg2",  "Hdc", "Shroom3", "Vit", "Rgs12", "Trp73"],
                 force_mask: bool = False):
        """
        Initialize a Xenium Experiment instance from transcripts table dataframe

        Parameters
        ----------
        transcripts_input : pd.DataFrame, str, or Path
            DataFrame of or path to transcripts table.
        ilastik_program_path : str or Path, optional
            Path to ilastik program. Default is None.
        dapi_high_res_image_path : str or Path, optional
            Path to high resolution DAPI tiff image. Default is None.
        output_dir : str or Path, optional
            Path to output directory. Default is None.
        ventricle_genes_list : list, optional
            List of genes that mark ventricle boundaries. Default is pre-set list.
        force_mask : bool, optional
            Whether to run ilastik workflow to generate mask. Default is False.
        """
        # Assign output dir if None is passed
        if output_dir is None:
            output_dir = Path(os.getcwd(), 'qc_output')

        # Create output directory if does not exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Assign ilastik output paths
        self.transcripts_image_path = Path(output_dir, 'transcripts.tiff')
        self.transcripts_mask_path = Path(output_dir, 'transcripts_mask.tiff')
        self.dapi_image_path = Path(output_dir, 'dapi.tiff')
        self.dapi_mask_path = Path(output_dir, 'dapi_mask.tiff')
        self.detachment_mask = Path(output_dir, 'detachment_mask.tiff')
        self.ventricle_image_path = Path(output_dir, 'ventricles.tiff')
        self.ventricle_mask_path = Path(output_dir, 'ventricles_mask.tiff')
        self.damage_mask_path = Path(output_dir, 'damage_mask.tiff')
        self.pixel_classification_path = Path(output_dir, 'pixel_classification.tiff')

        # Assign data paths as attributes
        self.ilastik_program_path = ilastik_program_path
        self.dapi_high_res_image_path = dapi_high_res_image_path
        self.output_dir = output_dir
        self.ventricle_genes_list = ventricle_genes_list

        # Paths to ilastik models assigned as attributes
        ilastik_models_dir = os.path.join(os.path.dirname(__file__), '..', 'ilastik_models')
        self.transcripts_mask_pixel_path = os.path.normpath(Path(ilastik_models_dir,
                                                             'TissueMaskPixelClassification_v1.0.ilp'))
        self.transcripts_mask_object_path = os.path.normpath(Path(ilastik_models_dir,
                                                                  'TissueMaskObjects_v1.1.ilp'))
        self.dapi_mask_pixel_path = os.path.normpath(Path(ilastik_models_dir,
                                                          'DapiMaskPixelClassification_Mouse.ilp'))
        self.dapi_mask_object_path = os.path.normpath(Path(ilastik_models_dir,
                                                           'DapiMaskObjectClassification_Mouse.ilp'))
        self.ventricle_mask_pixel_path = os.path.normpath(Path(ilastik_models_dir,
                                                               'VentriclesPixelClassification.ilp'))
        self.ventricle_mask_object_path = os.path.normpath(Path(ilastik_models_dir,
                                                                'VentriclesObjectClassification.ilp'))

        # Transcripts dataframe
        print('Processing transcripts datafame')
        transcripts = data_processing.process_input(transcripts_input)

        # Rename columns to match MERSCOPE convention
        transcripts.rename({'feature_name': 'gene', 'x_location': 'global_x', 'y_location': 'global_y',
                            'z_location': 'global_z', 'fov_name': 'fov'}, axis=1)
        transcripts['z_plane'] = transcripts['global_z'].astype(int)
        self.transcripts = transcripts

        # Adjust (x,y) locations
        self.transcript_counts = len(transcripts)

        # Counts per gene (including controls)
        self.counts_per_gene = self.transcripts.groupby('gene').size().to_dict()

        # Filter quality and controls
        self.filtered_transcripts = filter_low_quality_transcripts(self.transcripts)
        self.filtered_transcripts = filter_controls(self.filtered_transcripts, which='all')
        self.filtered_transcript_counts = len(self.filtered_transcripts)

        # Create FOVs dataframe
        print('Creating FOV dataframe')


@staticmethod
def scale_transcripts_xy()       


@staticmethod
def filter_low_quality_transcripts(transcripts: pd.DataFrame, val: int = 20):
    """
    Filters transcripts dataframe by Quality Value score

    Parameters
    ----------
    transcripts : pd.DataFrame
        Dataframe of detected transcripts
    val : int, optional
        QV score to filter by. Default is 20.

    Returns
    -------
    filtered_transcripts : pd.DataFrame
        Dataframe of detected transcripts with QV scores > 20

    Raises
    ------
    KeyError
        If `qv` is not a column in the transcripts dataframe
    """
    if 'qv' not in transcripts.columns:
        raise KeyError('transcripts dataframe must include "qv" column')

    return transcripts[transcripts['qv'] >= val]

@staticmethod
def filter_controls(transcripts: pd.DataFrame, which: str = 'all'):
    """
    Filters transcripts dataframe of control codewords

    Parameters
    ----------
    transcripts : pd.DataFrame
        Dataframe of detected transcripts
    which : ['all', 'NegControlCodeword', 'NegControlProbe', 'UnassignedCodeword', 'DeprecatedCodeword']
        Which control codewords to filter. Default is 'all'.

    Returns
    -------
    filtered_transcripts : pd.DataFrame
        Dataframe of detected trancsripts excluding specified control codewords

    Raises
    ------
    KeyError
        If `which` parameter is not a valid option
    """
    valid_which = ['all', 'NegControlCodeword', 'NegControlProbe', 'UnassignedCodeword', 'DeprecatedCodeword']
    if which not in valid_which:
        raise KeyError(f'"which" parameter is invalid. must be one of {valid_which}')

    if which == 'all':
        filtered_transcripts = transcripts[~transcripts['gene'].str.startswith(('NegControl',
                                                                                'Unassigned',
                                                                                'Deprecated'))]
    else:
        filtered_transcripts = transcripts[~transcripts['gene'].str.startswith(which)]

    return filtered_transcripts
