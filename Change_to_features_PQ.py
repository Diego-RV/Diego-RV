            _custom_log("Generating features quantitative table", logger)
            features_df = pd.DataFrame({
                'feature_id': [feature.id for feature in features],
                'average_mz': [feature.average_mz for feature in features],
                'average_mz_sigma': [feature.average_mz_sigma for feature in features],
                'average_rt': [feature.average_rt for feature in features],
                'average_rt_sigma': [feature.average_rt_sigma for feature in features],
                'average_rt_delta': [feature.average_rt_delta for feature in features],
                'total_height': [feature.total_height for feature in features],
                'total_volume': [feature.total_volume for feature in features],
                'monoisotopic_mz': [feature.monoisotopic_mz for feature in features],
                'monoisotopic_rt' : [feature.monoisotopic_rt for feature in features],
                'monoisotopic_height': [feature.monoisotopic_height for feature in features],
                'monoisotopic_volume': [feature.monoisotopic_volume for feature in features],
                'charge_state': [feature.charge_state for feature in features],
                'peak_id': [feature.peak_ids for feature in features],
            })