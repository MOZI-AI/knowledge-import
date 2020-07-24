# Finds GO node type from their namespace
import requests
from atomwrappers import *
import json

def find_type(go_term, go_ns_dict, go_ns=False):
    if not go_ns:
        if go_term in go_ns_dict.keys():
            go_ns = go_ns_dict[go_term]
        else:
            go_ns = request_api(go_term)
            go_ns_dict[go_term] = go_ns
    else:
        if not go_term in go_ns_dict.keys():
            go_ns_dict[go_term] = go_ns

    if go_ns in ["BP","biological_process"]:
        result = GoBPNode(go_term)
    elif go_ns in ["CC", "cellular_component"]:
        result = GoCCNode(go_term)
    elif go_ns in ["MF", "molecular_function"]:
        result = GoMFNode(go_term)
    else:
        result = False
    
    return go_ns_dict, result

def request_api(go_term):
    requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{}".format(go_term)
    result = requests.get(requestURL, headers={ "Accept" : "application/json"})
    if result.ok:
        result = json.loads(result.text)
        namespace = result['results'][0]['aspect']
        return namespace