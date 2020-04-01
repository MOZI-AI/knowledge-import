import pandas as pd
import sys
from datetime import date

def evaLink(elem1, elem2, predicate):
    return '(EvaluationLink \n' + '\t(PredicateNode "'+predicate+'")\n' + '\t(ListLink \n'+ '\t\t'+ elem1 + '\n\t\t'+ elem2 + '))\n'

def toAtomese(data):
    df =  pd.read_csv(data, sep='\t')
    df = df.dropna()
    print("Started importing ")
    with open("dataset/current_symbols_{}.scm".format(str(date.today())), "w") as f:
        for i in range(len(df)):
            prev = df.iloc[i]['Previous symbols']
            current = df.iloc[i]['Approved symbol'].strip().upper()
            name = df.iloc[i]['Approved name']
            if prev:
                for p in prev.split(","):
                    f.write(evaLink('(GeneNode "{}")'.format(p.strip().upper()), '(GeneNode "{}")'.format(current), "has_current_symbol"))
                f.write(evaLink('(GeneNode "{}")'.format(current), '(ConceptNode "{}")'.format(name), "has_name"))
    print("Done")
    
if __name__ == "__main__":
    # Source https://www.genenames.org/download/custom/ 
    # - Select columns Approved symbol, Approved name and Previous symbol 
    # - save into file, and give as an input
    try:
        data = sys.argv[1]
    except:
        data = "raw_data/custom_current.txt" 
    toAtomese(data)
