# umap_plots.py

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
from datetime import datetime
import numpy as np
import colorsys

class ColorPalettes:
    @staticmethod
    def create_extended_set2(n_colors=20):
        """Creates a color palette inspired by Set2 with better spacing."""
        # Get original Set2 colors to maintain the aesthetic
        original_set2 = sns.color_palette("Set2")
        
        # Calculate average HSV values from Set2 to maintain its character
        hsv_colors = [colorsys.rgb_to_hsv(*color) for color in original_set2]
        avg_saturation = sum(color[1] for color in hsv_colors) / len(hsv_colors)
        avg_value = sum(color[2] for color in hsv_colors) / len(hsv_colors)
        
        # Adjust these slightly to ensure colors are visible but not harsh
        base_saturation = avg_saturation * 0.95  # Slightly less saturated
        base_value = avg_value * 1.05  # Slightly brighter
        
        # Use golden ratio for optimal spacing
        golden_ratio = (1 + np.sqrt(5)) / 2
        golden_angle = 1 / golden_ratio
        
        colors = []
        for i in range(n_colors):
            # Calculate hue using golden ratio spacing
            hue = (i * golden_angle) % 1.0
            
            # Add slight variations to saturation and value
            # but keep them within the pleasing Set2-like range
            saturation = min(1.0, max(0.4, base_saturation + np.random.uniform(-0.05, 0.05)))
            value = min(1.0, max(0.7, base_value + np.random.uniform(-0.05, 0.05)))
            
            rgb = colorsys.hsv_to_rgb(hue, saturation, value)
            colors.append(rgb)
        
        # Optionally shuffle the colors to avoid similar colors being adjacent
        # while maintaining deterministic output
        rng = np.random.RandomState(42)
        rng.shuffle(colors)
        
        return colors

    @staticmethod
    def get_categorical_palette(data, minimum_colors=8):
        """
        Generate a color palette based on number of unique values.
        
        Parameters
        ----------
        data : pandas.Series
            Data to generate palette for
        minimum_colors : int
            Minimum number of colors to generate
        
        Returns
        -------
        list
            List of RGB colors
        """
        n_categories = len(data.unique())
        n_colors = max(n_categories, minimum_colors)
        return ColorPalettes.create_extended_set2(n_colors)

class UMAPPlotter:
    def __init__(self, 
                 dataset_name,
                 base_output_dir="figures",
                 default_continuous_cmap='viridis',
                 categorical_width_ratio=1.25,
                 categorical_height_ratio=1.0,
                 continuous_width_ratio=1.1,
                 continuous_height_ratio=1.0,
                 dpi = 300):
        """
        Initialize UMAP plotter with configurable parameters.
        
        Parameters
        ----------
        dataset_name : str
            Name of the dataset to be used in plot titles and file paths
        base_output_dir : str
            Base directory for saving figures
        default_continuous_cmap : str
            Default colormap for continuous variables
        categorical_width_ratio : float
            Width multiplier for categorical plots
        categorical_height_ratio : float
            Height multiplier for categorical plots
        continuous_width_ratio : float
            Width multiplier for continuous plots
        continuous_height_ratio : float
            Height multiplier for continuous plots
        """
        self.dataset_name = dataset_name
        self.base_output_dir = Path(base_output_dir)
        self.default_continuous_cmap = default_continuous_cmap
        self.categorical_width_ratio = categorical_width_ratio
        self.categorical_height_ratio = categorical_height_ratio
        self.continuous_width_ratio = continuous_width_ratio
        self.continuous_height_ratio = continuous_height_ratio
        self.dpi = dpi
        
        
    def _setup_output_dir(self, dataset_name):
        """Create dataset-specific output directory with timestamp."""
        timestamp = datetime.now().strftime("%Y%m%d")
        output_dir = self.base_output_dir / dataset_name / timestamp / "umaps"
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir

    def _get_value_range(self, adata, color_key, vmin, vmax):
        """Get value range for continuous variables."""
        if vmin == 'p01':
            vmin = np.percentile(adata.obs[color_key], 1)
        if vmax == 'p99':
            vmax = np.percentile(adata.obs[color_key], 99)
        return vmin, vmax

    def save_umap_plot(self,
                      adata,
                      color_key,
                      title,
                      filename_prefix,
                      base_figsize=3,
                      point_size=10,
                      cmap=None,
                      palette=None,
                      vmin=None,
                      vmax=None,
                      vcenter=None,
                      width_ratio=None,
                      height_ratio=None,
                      legend_loc=None,
                      save=True,
                      show=True):
        """
        Plot and optionally save/display UMAP with specified parameters.
        """
        sc.set_figure_params(dpi=self.dpi, frameon=False)
        
        if save:
            output_dir = self._setup_output_dir(self.dataset_name)
        
        # Check if the color_key exists in adata.obs
        if color_key not in adata.obs.columns:
            raise KeyError(f"Column '{color_key}' not found in adata.obs")
            
        is_categorical = isinstance(adata.obs[color_key].dtype, (object, 'category'))
        
        # Generate palette if categorical and no palette specified
        if is_categorical and palette is None:
            palette = ColorPalettes.get_categorical_palette(adata.obs[color_key])
        
        if width_ratio is None:
            width_ratio = self.categorical_width_ratio if is_categorical else self.continuous_width_ratio
        if height_ratio is None:
            height_ratio = self.categorical_height_ratio if is_categorical else self.continuous_height_ratio
        
        figsize = (base_figsize * width_ratio, base_figsize * height_ratio)
        
        if not is_categorical and cmap is None:
            cmap = self.default_continuous_cmap
            
        if legend_loc is None:
            legend_loc = 'right margin'
            
        if not is_categorical and (vmin == 'p01' or vmax == 'p99'):
            vmin, vmax = self._get_value_range(adata, color_key, vmin, vmax)
        
        fig = plt.figure(figsize=figsize)
        ax = plt.subplot(111)
        
        sc.pl.embedding(adata,
                       basis='X_umap',
                       color=color_key,
                       legend_loc=legend_loc,
                       legend_fontsize=8,
                       legend_fontoutline=2,
                       size=point_size,
                       alpha=0.5 if is_categorical else 0.7,
                       title=title,
                       ax=ax,
                       show=False,
                       cmap=cmap,
                       palette=palette,
                       vmin=vmin,
                       vmax=vmax,
                       vcenter=vcenter)
        
        plt.tight_layout()
        
        if save:
            for fmt in ['pdf', 'svg']:
                filepath = output_dir / f'{filename_prefix}.{fmt}'
                fig.savefig(filepath,
                           format=fmt,
                           dpi=self.dpi,
                           bbox_inches='tight',
                           transparent=True)
        
        if show:
            plt.show()
        else:
            plt.close()
            
        return str(output_dir) if save else None

    def plot_clusters(self,
                         adata,
                         cluster_key,
                         width_ratio=None,
                         **kwargs):
        """Plot cluster UMAP with specified cluster column."""
        palette = ColorPalettes.get_categorical_palette(adata.obs[cluster_key])
        
        # Use the instance variable if no width_ratio is provided
        if width_ratio is None:
            width_ratio = self.categorical_width_ratio
        
        return self.save_umap_plot(
            adata,
            color_key=cluster_key,
            title=f'{self.dataset_name} Clusters',
            filename_prefix='umap_clusters',
            palette=palette,
            width_ratio=width_ratio,  # Explicitly pass the width_ratio
            **kwargs
        )
    
    def plot_program(self,
                    adata,
                    program_key,
                    width_ratio=None,
                    height_ratio=None,
                    vmin='p01',
                    vmax='p99',
                    cmap=None,
                    **kwargs):
        """Plot program-specific UMAP."""
        if cmap is None:
            cmap = self.default_continuous_cmap
        
        # Use instance variables if no ratios are provided
        if width_ratio is None:
            width_ratio = self.continuous_width_ratio
        if height_ratio is None:
            height_ratio = self.continuous_height_ratio
        
        # Extract program number from the key if it exists
        try:
            program_number = program_key.split('_')[-1]
            title = f'Program {program_number}'
        except (AttributeError, IndexError):
            title = program_key
            
        return self.save_umap_plot(
            adata,
            color_key=program_key,
            title=title,
            filename_prefix=f'umap_{program_key}',
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            width_ratio=width_ratio,  # Explicitly pass the width_ratio
            height_ratio=height_ratio,  # Explicitly pass the height_ratio
            **kwargs
        )

    def plot_metadata(self,
                     adata,
                     metadata_key,
                     custom_title=None,
                     **kwargs):
        """Plot UMAP colored by specified metadata column."""
        title = custom_title or f'{self.dataset_name} - {metadata_key}'
        return self.save_umap_plot(
            adata,
            color_key=metadata_key,
            title=title,
            filename_prefix=f'umap_{metadata_key.lower().replace(" ", "_")}',
            **kwargs
        )