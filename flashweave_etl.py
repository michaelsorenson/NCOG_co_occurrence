import pandas as pd
import numpy as np
import re
import sys

EXCLUDE = [
    '201404_090.0_120.0_112',
    '201611_093.3_026.7_24_T',
    '201611_093.3_045.0_33_T',
    '201611_093.3_070.0_40_T',
    '201611_086.7_070.0_10_T',
    '201611_076.7_070.0_8_T',
    '201704_067.0_050.0_9.8',
    '201704_067.0_050.0_50.5',
    '201704_067.0_050.0_49.5',
    '201704_067.0_050.0_18.2',
    '201704_067.0_050.0_10.1',
    '201810_093.3_045.0_11',
    '201810_093.3_045.0_62',
    '202007_086.7_035.0_20'
]


def main(sixteens_path, eighteens_path, metadata_path, sixteens_outpath, eighteens_outpath, metadata_outpath, id_map_outpath):
    sixteenS = pd.read_csv(sixteens_path, sep='\t')
    eighteenS = pd.read_csv(eighteens_path, sep='\t')
    meta = pd.read_csv(metadata_path, sep='\t')
    sixteenS.columns = pd.Series(sixteenS.columns).apply(lambda x: x.replace('_16S_S2', '_S'))
    meta['sampleid'] = meta['sampleid'].apply(lambda x: "X" + x)
    # Drop rows and columns from 16S data
    sixteenSdroprows = []
    for i, row in sixteenS.iterrows():
        taxon = row['silva_Taxon'].lower()
        if 'eukaryota' in taxon or 'chloroplast' in taxon or 'mitochondria' in taxon:
            sixteenSdroprows.append(i)
    sixteenS_filtered = sixteenS.drop(sixteenSdroprows)
    sixteenSdropcols = []
    for col in sixteenS_filtered.columns:
        if re.fullmatch('X201[456].*\d\d?(_T)?', col) or 'blank' in col.lower() or 'ntc' in col.lower() \
        or 'mock' in col.lower() or col in EXCLUDE:
            sixteenSdropcols.append(col)
    sixteenSdropcols.extend(['silva_Taxon', 'silva_Confidence'])
    sixteenS_filtered = sixteenS_filtered.drop(sixteenSdropcols, axis=1)
    # Drop rows and columns from 18S data
    eighteenSdropcols = []
    for col in eighteenS.columns:
        if re.fullmatch('X201[456].*\d\d?(_T)?', col) or 'blank' in col.lower() or 'ntc' in col.lower() \
        or 'mock' in col.lower() or col in EXCLUDE:
            eighteenSdropcols.append(col)
    eighteenSdropcols.extend(['blank1_S', 'blank2_S', 'blank3_S', 'pr2_Taxon', 'pr2_Confidence', 'silva_Taxon',
                              'silva_Confidence'])
    if 'silva.full_Taxon' in eighteenS.columns:
        eighteenSdropcols.append('silva.full_Taxon')
    elif 'silva_full_Taxon' in eighteenS.columns:
        eighteenSdropcols.append('silva_full_Taxon')
    if 'silva.full_Confidence' in eighteenS.columns:
        eighteenSdropcols.append('silva.full_Confidence')
    elif 'silva_full_Confidence' in eighteenS.columns:
        eighteenSdropcols.append('silva_full_Confidence')
    eighteenS_filtered = eighteenS.drop(eighteenSdropcols, axis=1)
    # Drop rows from metadata
    metadroprows = []
    for i, sample in meta['sampleid'].iteritems():
        if re.fullmatch('X201[456].*\d\d?(_T)?', sample) or 'blank' in sample.lower() or 'ntc' in sample.lower() \
        or 'mock' in sample.lower() or sample in EXCLUDE:
            metadroprows.append(i)
    meta_filtered = meta.drop(metadroprows)
    # Sort columns
    sixteenS_filtered = sixteenS_filtered[np.sort(sixteenS_filtered.columns.values)]
    eighteenS_filtered = eighteenS_filtered[np.sort(eighteenS_filtered.columns.values)]
    # Drop exclusive columns
    nonexcl_cols = []
    for col in sixteenS_filtered.columns:
        if col not in eighteenS_filtered.columns:
            nonexcl_cols.append(col)
    sixteenS_filtered = sixteenS_filtered.drop(nonexcl_cols, axis=1)
    nonexcl_cols = []
    for col in eighteenS_filtered.columns:
        if col not in sixteenS_filtered.columns:
            nonexcl_cols.append(col)
    eighteenS_filtered = eighteenS_filtered.drop(nonexcl_cols, axis=1)
    # Drop rows with missing metadata
    dropids = list(set(sixteenS_filtered.columns[1:]).difference(set(meta_filtered['sampleid'])))
    sixteenS_filtered = sixteenS_filtered.drop(dropids, axis=1)
    eighteenS_filtered = eighteenS_filtered.drop(dropids, axis=1)
    meta_filtered = (
        meta_filtered.set_index('sampleid', drop=True).loc[sixteenS_filtered.columns[1:]]
    )
    # Drop OTUs that aren't above the Q3 in at least five samples
    def check_otus(df):# for each sample, compute Q3
        q3s = []
        for col in df.columns:
            sample = df[col]
            q3 = sample[sample != 0].describe()['75%']
            q3s.append(q3)
        # for each OTU, check if value > Q3 in at least 5 samples
        q3counts = []
        for i, row in df.iterrows():
            count = 0
            for j, val in enumerate(row):
                if val > q3s[j]:
                    count += 1
            q3counts.append(count)
        finalcounts = pd.Series(q3counts, index=df.index)
        return finalcounts[finalcounts >= 5]
    sixteenS_filtered = sixteenS_filtered.loc[check_otus(sixteenS_filtered.drop('Feature.ID', axis=1)).index.values]
    eighteenS_filtered = eighteenS_filtered.loc[check_otus(eighteenS_filtered.drop('Feature.ID', axis=1)).index.values]
    # Filter metadata to important columns
    meta_filtered = meta_filtered[['Cruise', 'Depthm', 'Distance', 'T_degC', 'Salnty', 'STheta', 'O2ml_L', 'PO4ug',
                                   'SiO3ug', 'NO3ug', 'NH3ug', 'ChlorA', 'Phaeop', 'MLD_Sigma', 'NCDepth', 'season',
                                   'year', 'sample_type', 'sample_type2', 'siex']]
    # Flip 16S and 18S and add potu prefix to 16S and eotu prefix to 18S (for "prokaryote" and "eukaryote" respectively)
    potus = pd.Series(sixteenS_filtered.index).apply(lambda x: 'potu_' + str(x)).values
    eotus = pd.Series(eighteenS_filtered.index).apply(lambda x: 'eotu_' + str(x)).values
    id_map = pd.DataFrame({'Feature.ID': sixteenS_filtered['Feature.ID'].values, 'otu_id': potus})
    id_map = pd.concat([id_map, pd.DataFrame({'Feature.ID': eighteenS_filtered['Feature.ID'].values,
                                              'otu_id': eotus})])
    id_map.to_csv(id_map_outpath, sep='\t', index=False)
    sixteenS_filtered = sixteenS_filtered.drop('Feature.ID', axis=1).T
    eighteenS_filtered = eighteenS_filtered.drop('Feature.ID', axis=1).T
    sixteenS_filtered.columns = potus
    eighteenS_filtered.columns = eotus
    # Drop rows with NA metadata
    droprows = np.array(list(set(meta_filtered.index).difference(set(meta_filtered.dropna().index))))
    sixteenS_filtered = sixteenS_filtered.drop(droprows)
    eighteenS_filtered = eighteenS_filtered.drop(droprows)
    meta_filtered = meta_filtered.drop(droprows)
    # re-encode sample_type columns to make it quantitative
    meta_filtered['sample_type'] = meta_filtered['sample_type'].replace({'515': 'T_515', '170': 'T_170'})
    meta_filtered['sample_type2'] = meta_filtered['sample_type2'].replace({'515': 'T_515', '170': 'T_170'})
    # save filtered datasets to csv's
    sixteenS_filtered.to_csv(sixteens_outpath, sep='\t')
    eighteenS_filtered.to_csv(eighteens_outpath, sep='\t')
    meta_filtered.to_csv(metadata_outpath, sep='\t')
    

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print('Usage: python etl.py sixteens_path eighteens_path metadata_path ' +
              'sixteens_outpath eighteens_outpath metadata_outpath')
        exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    