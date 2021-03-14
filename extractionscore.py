import urllib.parse, urllib.request
import math
import os
import rdkit.Chem as Chem
from rdkit.Chem.EnumerateStereoisomers import GetStereoisomerCount, EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from rdkit.Chem.Draw import ReactionToImage
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
from PIL import Image
import psutil
from draw import ReactionStringToImage

def get_logP(smiles):
    """
    Input: SMILES representation of a molecule
    Output: float value of logP calculated using ALOGPS
    """
    url = "http://www.vcclab.org/web/alogps/calc?" + urllib.parse.urlencode({"SMILES":smiles})
    results = urllib.request.urlopen(url)
    html_results = results.read()
    html_results_tokens = html_results.split()
    logP_token = html_results_tokens[4].decode("utf-8") #this is where the logP value is in the results page

    try:
        logP = float(logP_token)
        return logP
    except Exception: #For ions or SMILES strings the website doesn't recognize
        return math.inf


def get_canonical_smiles(smiles):
    """
    Input: SMILES string for a molecule
    Output: SMILES string after converting to and from an rdkit molecule for canonicalization (empty string if input is invalid)
    """
    molecule = Chem.rdmolfiles.MolFromSmiles(smiles)

    if molecule is None:
        return ""

    return Chem.rdmolfiles.MolToSmiles(molecule,isomericSmiles=True)


def get_step_score(species_list,target_species):
    """
    Input: 
        List of SMILES strings of species present in a synthetic step, e.g. ["CCO","CNCC"]
        SMILES String for target species, i.e. desired product
    Output: ExtractionScore for synthetic step (float)
    """
    global extractionscore_detailed_log
    target_logP = get_logP(target_species)
    extractionscore_detailed_log.write("{} {:.2f} (target)\n".format(target_species,target_logP))
    logPs_below_target = []
    logPs_above_target = []
    separability_threshold = 3 #units of logP

    for species in species_list:
        if species != target_species:
            current_logP = get_logP(species)
            extractionscore_detailed_log.write("{} {:.2f}\n".format(species,current_logP))
            if current_logP < target_logP:
                logPs_below_target.append(current_logP)
            if current_logP > target_logP:
                logPs_above_target.append(current_logP)
    
    if len(logPs_below_target) >= len(logPs_above_target):
        for logP in logPs_below_target:
            if target_logP - logP >= separability_threshold:
                number_separable += 1
        fraction_separable = number_separable/len(logPs_below_target)
    else:
        for logP in logPs_above_target:
            if logP - target_logP >= separability_threshold:
                number_separable += 1
        fraction_separable = number_separable/len(logPs_above_target)

    return fraction_separable


def insert_spaces(smiles):
    """
    Input: SMILES string
    Output: SMILES string with inserted spaces to work with Molecular Transformer
    """
    list_of_two_letter_elements = ["Li","Na","Mg","Al","Si","Cl","Ar","Ca","Ti","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Sr","Zr","Mo","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","Cs","Ba","La","Ce","Nd","Sm","Eu","Gd","Tb","Er","Yb","Ta","Re","Os","Ir","Pt","Au","Hg","Pb","Bi","At"]
    spaced_smiles_list = []
    index = 0
    while index < len(smiles):
        if smiles[index] == "[":
            end_bracket_offset = smiles[index:].find("]")
            spaced_smiles_list.append(smiles[index:index+end_bracket_offset+1])
            index += end_bracket_offset + 1
        elif smiles[index:index+2] in list_of_two_letter_elements:
            spaced_smiles_list.append(smiles[index:index+2])
            index += 2
        else:
            spaced_smiles_list.append(smiles[index])
            index += 1
    
    return " ".join(spaced_smiles_list)


def get_predicted_products(reactants,product):
    """
    Input: SMILES representation of reactant species, delimited by '.', e.g. "CCO.CNCC" and product
    Output: List of SMILES strings of predicted products, e.g. ["CCO","CNCC"], using MT_STEREO
    """
    input_file = open("temp_input.txt","w")
    input_file.write(insert_spaces(reactants))
    input_file.close()
    os.system("python translate.py -model models/STEREO_mixed_augm_model_average_20.pt -src temp_input.txt -output temp_output.txt -batch_size 1 -replace_unk -max_length 200 -n_best 2")
    output_file = open("temp_output.txt","r")
    predictions = output_file.readlines()
    output_file.close()
    individual_reactants = reactants.split(".")
    canonicalized_reactants = [get_canonical_smiles(reactant) for reactant in individual_reactants]
    untokenized_predictions = ["".join(prediction.split()) for prediction in predictions]
    verified_predictions = []

    #Only accept predictions with valid SMILES that don't appear in the reactants
    for prediction in untokenized_predictions:
        if (Chem.rdmolfiles.MolFromSmiles(prediction) is not None) and (prediction not in canonicalized_reactants):
            verified_predictions.append(prediction)

    if get_canonical_smiles(product) not in verified_predictions and len(verified_predictions) > 1:
        #Since we use MT_STEREO and there is currently no way to get more than 1 prediction using IBM RXN,
        #sometimes MT_STEREO can't predict the major product even though IBM RXN can, which is what matters for
        #the reaction being available for retrosynthesis. In this case we just take the top MT_STEREO prediction.
        verified_predictions = [verified_predictions[0]]

    verified_predictions_with_isomers = get_verified_predictions_with_isomers(canonicalized_reactants,verified_predictions,product)

    return verified_predictions_with_isomers


def close_image():
    """
    Closes Image.show()
    """
    for proc in psutil.process_iter():
        if proc.name() == "display":
            proc.kill()


def are_mutual_stereoisomers(smiles1,smiles2):
    """
    Input: two SMILES strings
    Output: True if the two molecules are stereoisomers of each other, False otherwise
    """
    molecule1 = Chem.rdmolfiles.MolFromSmiles(smiles1)
    molecule2 = Chem.rdmolfiles.MolFromSmiles(smiles2)
    stereoambiguous_smiles1 = Chem.rdmolfiles.MolToSmiles(molecule1,isomericSmiles=False)
    stereoambiguous_smiles2 = Chem.rdmolfiles.MolToSmiles(molecule2,isomericSmiles=False)

    return stereoambiguous_smiles1 == stereoambiguous_smiles2


def get_stereoisomers(smiles):
    """
    Input: SMILES string for a molecule
    Output: set of SMILES strings of stereoisomers for the molecule
    """
    opts = StereoEnumerationOptions(tryEmbedding=True,unique=True)
    molecule = Chem.rdmolfiles.MolFromSmiles(smiles)
    
    return {Chem.rdmolfiles.MolToSmiles(i,isomericSmiles=True) for i in list(EnumerateStereoisomers(molecule,options=opts))}


def have_same_isomers(smiles1,smiles2):
    """
    Input: two SMILES strings
    Output: True if the sets of all possible stereoisomers of the two molecules are equal, False otherwise
    """
    molecule1_isomers_smiles = get_stereoisomers(smiles1)
    molecule2_isomers_smiles = get_stereoisomers(smiles2)

    return molecule1_isomers_smiles == molecule2_isomers_smiles


def get_distinct_prediction_isomers(reactants,prediction,product):
    """
    Input:
        -List of SMILES strings for reactants
        -SMILES string for prediction
        -SMILES string for intended product
    Output: List of stereochemically explicit SMILES strings for isomers 
        of the prediction that are not isomers of the product
    """
    prediction_isomers = get_stereoisomers(prediction)
    product_isomers = get_stereoisomers(product)
    distinct_isomer_list = []

    for isomer in prediction_isomers:
        if isomer not in product_isomers:
            proposed_reaction_smarts = ".".join(reactants) + ">>" + isomer
            illustration = ReactionStringToImage(proposed_reaction_smarts)
            illustration.show()
            accept = input("Accept this product? [y]/n\n")
            if accept == "" or "y":
                distinct_isomer_list.append(isomer)
            close_image()

    return distinct_isomer_list


def get_verified_predictions_with_isomers(reactants,predictions,product):
    """
    Input: 
        -List of SMILES strings for reactants
        -List of SMILES strings for predictions
        -Intended product as a single SMILES string (may have multiple species)
    Output: List of stereochemically explicit SMILES strings for predictions
    """
    #Convert the product to a list, breaking it up into several species if necessary
    product = product.replace("~",".")
    product_list = product.split(".")

    enumerated_product_list = []
    enumerated_product_list.extend(product_list)

    for prediction in predictions:
        prediction_matches_a_product = False
        for product_species in product_list:
            if are_mutual_stereoisomers(prediction,product_species):
                prediction_matches_a_product = True
                matched_product = product_species
                break
        if prediction_matches_a_product:
            if not have_same_isomers(prediction,matched_product):
                prediction_isomers = get_distinct_prediction_isomers(reactants,prediction,matched_product)
                enumerated_product_list.extend(prediction_isomers)
        else:
            enumerated_product_list.append(prediction)

    return enumerated_product_list


def get_major_product(product_list):
    """
    Input: list of product SMILES strings
    Output: SMILES of heaviest species

    This function is needed in cases where a counterion is included in the intended products.
    This often happens e.g. with an amine and HCl. Since the logP calculator can't parse
    SMILES strings with "~" or "." in them, we need to go from e.g. "OCCNCCO~Cl" to "OCCNCCO".
    After the former is broken into a list (["OCCNCCO","Cl"]), this function uses molecular weight
    to determine which species is the "primary" product.
    """
    max_molecular_weight = -math.inf
    for product in product_list:
        molecule = Chem.rdmolfiles.MolFromSmiles(product)
        molecular_weight = CalcExactMolWt(molecule)
        if molecular_weight > max_molecular_weight:
            max_molecular_weight = molecular_weight
            heaviest_species = product

    return heaviest_species


def get_route_score(steps):
    """
    Input: List of SMILES reaction representations for each step of a synthetic route, e.g. ["C.O>>CO","CO.C>>COC"]
    Output: ExtractionScore for the synthetic route (float)
    """
    global extractionscore_detailed_log
    number_of_steps = len(steps)
    step_scores = []
    
    for step in steps:
        extractionscore_detailed_log.write(step)
        reactants, product = step.strip().split(">>")
        predicted_products = get_predicted_products(reactants,product) #includes intended product
        reactant_list = reactants.split(".")
        full_species_list = list(set(reactant_list + predicted_products))
        #Handle counterions or other additional species in the product
        product = product.replace("~",".")
        product_list = product.split(".")
        if len(product_list) > 1:
            product = get_major_product(product_list)
        else:
            product = product_list[0]
        step_score = get_step_score(full_species_list,product)
        step_scores.append(step_score)
        extractionscore_detailed_log.write("Step score: {:.2f}\n".format(step_score))
    
    return 0.5 * ( sum(step_scores) / number_of_steps + 1 / math.sqrt(number_of_steps) )


input_file = open("route_pairs.txt","r")
input_lines = input_file.readlines()
input_file.close()
target_indices = [i for i,line in enumerate(input_lines) if "Target" in line]
extractionscore_detailed_log = open("extractionscorelog.txt","x")
for target_index in target_indices:
    remaining_lines = input_lines[target_index:]
    starts_of_remaining_routes = [i for i,line in enumerate(remaining_lines[1:]) if "Target" in line]

    if starts_of_remaining_routes:
        end_of_route = starts_of_remaining_routes[0] + 1
        lines_to_use = remaining_lines[:end_of_route]
    else:
        lines_to_use = remaining_lines

    target_id = " ".join(lines_to_use[0].split()[1:])
    print("Working on target " + target_id)
    extractionscore_detailed_log.write("Target " + target_id + "\n")
    end_of_industrial_route = [i for i,line in enumerate(lines_to_use) if "Suggested" in line][0]
    industrial_route = lines_to_use[2:end_of_industrial_route]
    suggested_route = lines_to_use[end_of_industrial_route+1:]
    extractionscore_detailed_log.write("Industrial:\n")
    industrial_extractionscore = get_route_score(industrial_route)
    extractionscore_detailed_log.write("IBM RXN:\n")
    suggested_extractionscore = get_route_score(suggested_route)

    extractionscore_detailed_log.write("For target {} industrial score = {:.2f}, IBM RXN score = {:.2f}\n".format(target_id,industrial_extractionscore,suggested_extractionscore))

extractionscore_detailed_log.close()