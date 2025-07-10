def TOF_parameters(instrument='TOF', avg_fwhm_rt=4):
    if instrument == 'TOF':
        TOF_parameters = {
            #
            # Instrument configuration.
            #
            'instrument_type': 'TOF',
            'resolution_ms1': 30000,
            'resolution_msn': 30000,
            'reference_mz': 400,
            'avg_fwhm_rt': avg_fwhm_rt,
            #
            # Meshing.
            #
            'num_samples_mz': 5,
            'num_samples_rt': 5,
            'smoothing_coefficient_mz': 0.4,
            'smoothing_coefficient_rt': 0.4,
            #
            # Warp2D.
            #
            'warp2d_slack': 30,
            'warp2d_window_size': 100,
            'warp2d_num_points': 2000,
            'warp2d_rt_expand_factor': 0.2,
            'warp2d_peaks_per_window': 100,
            #
            # MetaMatch.
            #
            'metamatch_fraction': 0.7,
            'metamatch_n_sig_mz': 1.5,
            'metamatch_n_sig_rt': 1.5,
            #
            # Feature detection.
            #
            'feature_detection_charge_states': [5, 4, 3, 2, 1],
            #
            # Other.
            #
            'max_peaks': 1000000,
            'polarity': 'both',
            'min_mz': 0,
            'max_mz': 100000,
            'min_rt': 0,
            'max_rt': 100000,
            #
            # Annotation linking.
            #
            'link_n_sig_mz': 3,
            'link_n_sig_rt': 3,
            #
            # Identification.
            #
            # Keep only the max rank PSM.
            'ident_max_rank_only': True,
            # Whether to filter the PSM that don't pass the FDR threshold.
            'ident_require_threshold': True,
            # Ignore PSMs that have been marked as decoy.
            'ident_ignore_decoy': True,
            #
            # Quality.
            #
            'similarity_num_peaks': 2000,
            # Options: Any 'seaborn' supported palette style, like:
            #          'husl', 'crest', 'Spectral', 'flare', 'mako', etc.
            'qc_plot_palette': 'husl',
            # Options: 'png', 'pdf', 'eps'...
            'qc_plot_extension': 'png',
            # Options: 'dynamic', [0.0-1.0]
            'qc_plot_fill_alpha': 'dynamic',
            'qc_plot_line_alpha': 0.5,
            'qc_plot_scatter_alpha': 0.3,
            'qc_plot_scatter_size': 2,
            'qc_plot_min_dynamic_alpha': 0.1,
            'qc_plot_per_file': False,
            # Options: 'fill', 'line'
            'qc_plot_line_style': 'fill',
            # Plot style config.
            'qc_plot_dpi': 300,
            'qc_plot_font_family': 'sans-serif',
            'qc_plot_font_size': 7,
            'qc_plot_fig_size_x': 7.08661,
            'qc_plot_fig_size_y': 7.08661/1.618034,
            'qc_plot_fig_legend': False,
            'qc_plot_mz_vs_sigma_mz_max_peaks': 200000,
            #
            # Quantitative table generation.
            #
            # Options: 'height', 'volume'
            'quant_isotopes': 'height',
            # Options: 'monoisotopic_height', 'monoisotopic_volume',
            #          'total_height', 'total_volume',
            #          'max_height', 'max_volume',
            'quant_features': 'max_height',
            # Whether to remove feature annotations where the charge state of
            # the detected feature doesn't match the one given by the
            # identification engine.
            'quant_features_charge_state_filter': True,
            # Options: 'theoretical_mz', 'msms_event'
            'quant_ident_linkage': 'msms_event',
            # Whether to obtain a consensus sequence and proteins on identifications.
            'quant_consensus': True,
            # Demand a minimum number of files with identification per cluster.
            'quant_consensus_min_ident': 2,
            # Whether to store all the annotations prior to cluster aggregation.
            'quant_save_all_annotations': True,
            # Minimum number of peptides necessary to consider a protein for
            # quantification.
            'quant_proteins_min_peptides': 1,
            # Wether to remove proteins whose peptides are entirely contined
            # within another with a longest number of evidence peptides.
            'quant_proteins_remove_subset_proteins': True,
            # In case a peptide can't be assigned to a unique protein as 'razor'
            # we can choose to use them regardless in all instances where they
            # would if they were to be considered razor or to ignore them.
            'quant_proteins_ignore_ambiguous_peptides': True,
            # Protein inference method:
            #     - 'unique': Only unique peptides will be used for
            #                 quantification.
            #     - 'razor': Unique and peptides assigned as most likely through
            #                the Occam's razor constrain.
            #     - 'all': All peptides will be used for quantification for all
            #              protein groups. Thus shared peptides will be used more
            #              than once.
            'quant_proteins_quant_type': 'razor',
        }
        return TOF_parameters