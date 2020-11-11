import pandas as pd
import requests
import io
import numpy as np

# Splits string entry to rows and expand the dataframe
def explode(df, lst_cols, fill_value='NA', preserve_index=False):
    # make sure `lst_cols` is list-alike
    if (lst_cols is not None
        and len(lst_cols) > 0
        and not isinstance(lst_cols, (list, tuple, np.ndarray, pd.Series))):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)
    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()
    # preserve original index values    
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    res = (pd.DataFrame({
                col:np.repeat(df[col].values, lens)
                for col in idx_cols},
                index=idx)
             .assign(**{col:np.concatenate(df.loc[lens>0, col].values)
                            for col in lst_cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        res = (res.append(df.loc[lens==0, idx_cols], sort=False)
                  .fillna(fill_value))
    # revert the original index order
    res = res.sort_index()
    # reset index if requested
    if not preserve_index:        
        res = res.reset_index(drop=True)
    return res

# Get the current vs previous gene symbol mapping from HGNC
url = "https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&"+ \
    "col=gd_prev_sym&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&" + \
    "order_by=gd_app_sym_sort&format=text&submit=submit"
current_symbols_df = []
previous_symbols = []

def build_df():
    global current_symbols_df
    global previous_symbols

    data = requests.get(url).content
    current_symbols_df = pd.read_csv(io.StringIO(data.decode('utf-8')), sep="\t")
    current_symbols_df = current_symbols_df.fillna("NA")

    # Split string entry(morethan one previous Gene symbols listed as a string) to separate rows
    current_symbols_df["Previous symbols"] = current_symbols_df["Previous symbols"].str.split(",")
    current_symbols_df = explode(current_symbols_df, ["Previous symbols"])

    # Convert Gene symbols to uppercase and remove space
    current_symbols_df["Previous symbols"] = current_symbols_df["Previous symbols"].str.upper().str.strip()
    current_symbols_df["Approved symbol"] = current_symbols_df["Approved symbol"].str.upper()

    previous_symbols = current_symbols_df["Previous symbols"].unique()

def get_current_symbol(gene):
    if len(current_symbols_df) == 0:
        build_df()
    gene = str(gene).upper()
    if gene in previous_symbols:
        result = current_symbols_df[current_symbols_df["Previous symbols"] == gene]["Approved symbol"].values
        if len(result) > 0:
            gene = result[0]
    return gene
    