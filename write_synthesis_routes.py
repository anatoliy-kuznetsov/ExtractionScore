from rxn4chemistry import RXN4ChemistryWrapper
from api_info import API_KEY, PROJECT_ID
import rdkit.Chem as Chem
import time


def wait():
    if request_count % 5 == 0 and request_count != 0:
        time.sleep(60)
    else:
        time.sleep(2)


def collect_reactions(path):
    """
    Input: synthetic path represented as dictionary
    Output: list of SMILES representations of reactions in the path
    """
    reactions = []

    if 'children' in path and len(path['children']):
        reactions.append('{}>>{}'.format('.'.join([node['smiles'] for node in path['children']]),path['smiles']))
    
    for node in path['children']:
        reactions.extend(collect_reactions(node))

    return reactions


def get_maximum_path_score(paths):
    """
    Input: list of synthetic paths represented as dictionaries
    Output: maximum path confidence score
    """
    maximum_score = 0

    for path in paths:
        if path["confidence"] > maximum_score:
            maximum_score = path["confidence"]
    
    return maximum_score


def get_intermediates_to_exclude(industrial_route):
    """
    Input: industrially practiced synthesis route as a list of SMILES reaction strings
    Output: SMILES string of intermediates to exclude delimited by "."
    """
    intermediates_string = ""
    for step in industrial_route[:-1]:
        step_product = step.strip().split(">>")[1]
        intermediates_string += step_product + "."
    intermediates_string = intermediates_string[:-1] #remove trailing "."

    return intermediates_string


def get_retrosynthesis(industrial_route):
    """
    Input: industrially practiced synthesis route as a list of SMILES reaction strings
    Output: list of SMILES strings in the highest-scoring retrosynthetic route
        -> If multiple routes have the highest score, the one with the smallest index is returned
    """
    global request_count
    target = industrial_route[-1].split(">>")[1]
    intermediates_to_exclude = get_intermediates_to_exclude(industrial_route)
    results_successfully_retrieved = False
    response = rxn4chemistry_wrapper.predict_automatic_retrosynthesis(product=target,max_steps=10,exclude_smiles=intermediates_to_exclude)
    time.sleep(60)
    request_count = 0
    start_of_retrosynthesis_attempt = time.time()

    while not results_successfully_retrieved:
        try:
            results = rxn4chemistry_wrapper.get_predict_automatic_retrosynthesis_results(response['prediction_id'])
            if results["status"] == "SUCCESS":
                results_successfully_retrieved = True
            if time.time() - start_of_retrosynthesis_attempt > 3600:
                raise Exception("Couldn't return retrosynthetic route within one hour. Moving on to next target.")
            request_count += 1
            wait()
        except KeyError: #when the website doesn't respond
            print("Website didn't respond. Trying again in 5 minutes...")
            time.sleep(300)
            request_count = 0
        except Exception as e:
            print("Unexpected error on target " + target + ": " + "{0}".format(e))
            print("Suggested route could not be returned.")
            request_count += 1
            wait()
            return []

    maximum_path_score = get_maximum_path_score(results["retrosynthetic_paths"])

    for path in results["retrosynthetic_paths"]:
        if path["confidence"] == maximum_path_score:

            return collect_reactions(path)


def get_target_molecule_name(line):
    """
    Input: line with target molecule name from the file with industrial synthesis routes
    Output: string with just the target molecule name
    """
    tokens = line.split()
    name_with_colon = " ".join(tokens[1:-2])

    return name_with_colon[:-1]


def get_industrial_routes(lines):
    """
    Input: lines from the file with industrial synthesis routes
    Output: list of industrial routes. Each route is a list with the following format:
        Index 0: string name of target molecule
        Remaining indices: strings for each reaction in the synthesis route
    """
    industrial_routes = []
    industrial_route_indices = [i for i,line in enumerate(lines) if "step" in line]

    for index in industrial_route_indices:
        current_route = []
        target_molecule_name = get_target_molecule_name(lines[index])
        current_route.append(target_molecule_name)
        number_of_steps = int(lines[index].split()[-2])
        for line in lines[index+1:index+number_of_steps+1]:
            current_route.append(line.split()[-1])
        industrial_routes.append(current_route)
    
    return industrial_routes


def get_canonical_smiles(smiles):
    """
    Input: SMILES string for a molecule
    Output: SMILES string after converting to and from an rdkit molecule for canonicalization
    """
    molecule = Chem.rdmolfiles.MolFromSmiles(smiles)

    if molecule is None:
        print(smiles)
        exit()

    return Chem.rdmolfiles.MolToSmiles(molecule,isomericSmiles=True)


def are_all_products_correctly_predicted(industrial_route):
    """
    Input: industrial synthesis route as a list of SMILES reaction representations
    Output: True if IBM RXN correctly predicts the major product for each step, False otherwise
    """
    global request_count
    are_all_correct = True

    for reaction in industrial_route:
        reaction_tokens = reaction.split(">>")
        reactants = reaction_tokens[0]
        product = reaction_tokens[1]
        response = rxn4chemistry_wrapper.predict_reaction(reactants)
        request_count += 1
        wait()
        results = rxn4chemistry_wrapper.get_predict_reaction_results(response["prediction_id"])
        request_count += 1
        wait()
        predicted_product = results["response"]["payload"]["attempts"][0]["smiles"].split(">>")[1]
        if get_canonical_smiles(predicted_product) != get_canonical_smiles(product):
            are_all_correct = False
            print("Failed to correctly predict product of the following reaction:")
            print(reaction)

    return are_all_correct


def has_only_one_heavy_atom(smiles):
    """
    Input: SMILES string for a molecule
    Output: True if molecule has just one heavy atom, False otherwise
    """
    molecule = Chem.rdmolfiles.MolFromSmiles(smiles)

    return molecule.GetNumHeavyAtoms() == 1


def get_starting_material_paths(starting_materials_list):
    """
    Input: list of SMILES strings for starting materials in a synthetic route
    Output: list of SMILES reaction strings for materials which are not available in the emolecules database
    """
    global request_count
    additional_reactions = []

    for starting_material in starting_materials_list:
        if has_only_one_heavy_atom(starting_material): #Don't need to run retrosynthesis on these to know they are available
            break
        response = rxn4chemistry_wrapper.predict_automatic_retrosynthesis(product=starting_material,max_steps=10,exclude_target_molecule=False)
        time.sleep(60)
        request_count = 0
        results = rxn4chemistry_wrapper.get_predict_automatic_retrosynthesis_results(response['prediction_id'])
        request_count += 1
        wait()

        while results["status"] != "SUCCESS":
            results = rxn4chemistry_wrapper.get_predict_automatic_retrosynthesis_results(response['prediction_id'])
            request_count += 1
            wait()

        maximum_path_score = get_maximum_path_score(results["retrosynthetic_paths"])

        for path in results["retrosynthetic_paths"]:
            if path["confidence"] == maximum_path_score:
                path_reactions = collect_reactions(path)
                if not (len(path_reactions) == 1 and ">>" not in path_reactions[0]): #molecule is unavailable
                    additional_reactions.extend(path_reactions)
                break


    return additional_reactions


def expand_starting_materials(industrial_route):
    """
    Input: industrial synthesis route as a list of SMILES reaction representations
    Output: industrial synthesis route with appended reactions to form the starting materials,
        if they cannot be found in the eMolecules database
    """
    starting_materials_list = []
    products_list = []

    for step in industrial_route:
        reaction_tokens = step.split(">>")
        reactants = reaction_tokens[0]
        product = reaction_tokens[1]
        reactant_tokens = reactants.split(".")
        starting_materials_list.extend(reactant_tokens)
        products_list.append(get_canonical_smiles(product))

    starting_materials_list = list(set(starting_materials_list)) #remove duplicates

    for starting_material in starting_materials_list:
        if get_canonical_smiles(starting_material) in products_list:
            starting_materials_list.remove(starting_material)

    industrial_route.extend(get_starting_material_paths(starting_materials_list))

    return list(set(industrial_route))
    

if __name__ == "__main__":
    request_count = 0
    rxn4chemistry_wrapper = RXN4ChemistryWrapper(api_key=API_KEY)
    rxn4chemistry_wrapper.set_project(PROJECT_ID)
    request_count += 1
    wait()
    industrial_route_file = open("industrial_synthesis_routes.txt","r")
    industrial_route_lines = industrial_route_file.readlines()
    industrial_route_file.close()
    industrial_routes = get_industrial_routes(industrial_route_lines)

    output_file = open("route_pairs.txt","x")

    for index,industrial_route in enumerate(industrial_routes):
        if are_all_products_correctly_predicted(industrial_route[1:]):
            print("Working on route " + str(index+1) + ": " + industrial_route[0])
            suggested_route = get_retrosynthesis(industrial_route[1:])
            expanded_industrial_route = expand_starting_materials(industrial_route[1:])
            output_file.write("Route " + str(index+1) + ": " + industrial_route[0] + "\n")
            output_file.write("Industrial:\n")
            for step in industrial_route[1:]:
                output_file.write(step + "\n")
            output_file.write("IBM RXN:\n")
            for step in suggested_route:
                output_file.write(step + "\n")

    output_file.close()