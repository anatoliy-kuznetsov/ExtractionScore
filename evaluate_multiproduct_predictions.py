import rdkit.Chem as Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions, GetStereoisomerCount
enumeration_options = StereoEnumerationOptions(tryEmbedding=True,unique=True)


def get_individual_smiles_lists(untokenized_sequence,is_predictions):
    """
    Input: SMILES sequence of molecules separated by "."
    Output: 
        if is_predictions:
            List of lists of SMILES strings for all unique stereoisomers in the sequence
            e.g. if the top 5 predictions are "CO.CCO.FC(O)(C)Br.CCCO.CCCCO", the output should be
            [["CO"],["CCO"],["F[C@@](O)(C)Br","F[C@](O)(C)Br"],["CCCO"],["CCCCO"]]
            This way, both stereoisomers are included in the 3rd prediction
        else:
            List of SMILES strings for all unique stereoisomers in the products
            e.g. the above example returns ["CO","CCO","F[C@@](O)(C)Br","F[C@](O)(C)Br","CCCO","CCCCO"]

        This facilitates easier top-k coverage comparisons
    """
    individual_smiles_list = []
    original_individual_smiles_list = untokenized_sequence.split(".")
    for individual_smiles in original_individual_smiles_list:
        molecule = Chem.rdmolfiles.MolFromSmiles(individual_smiles)
        try:
            num_stereoisomers = GetStereoisomerCount(molecule)
            if num_stereoisomers > 1:
                isomers = tuple(EnumerateStereoisomers(molecule,options=enumeration_options))
                isomer_smiles_list = sorted(Chem.rdmolfiles.MolToSmiles(isomer,isomericSmiles=True) for isomer in isomers)
            else:
                isomer_smiles_list = [Chem.rdmolfiles.MolToSmiles(molecule,isomericSmiles=True)]
        except Exception:
            isomer_smiles_list = []

        if is_predictions:
            individual_smiles_list.append(isomer_smiles_list)
        else:
            individual_smiles_list.extend(isomer_smiles_list)

    return individual_smiles_list


def get_coverage_fraction(predictions,targets):
    """
    Input: predictions and targets as lists of SMILES strings
    Output: fraction of target species found in the predictions
    """
    number_predicted = 0
    for target in targets:
        if target in predictions:
            number_predicted += 1

    return number_predicted/len(targets)


def get_false_positive_fraction(predictions,targets):
    """
    Input: predictions and targets as lists of SMILES strings
    Output: fraction of predictions not found in the targets
    """
    number_of_false_positives = 0
    if len(predictions) == 0:
        return 0
    for prediction in predictions:
        if prediction not in targets:
            number_of_false_positives += 1

    return number_of_false_positives/len(predictions)


prediction_file = open("multiproduct_predictions_strict.txt","r")
predictions = prediction_file.readlines()
prediction_file.close()
targets_file = open("multiproduct_target_strict.txt","r")
targets = targets_file.readlines()
targets_file.close()

output_file = open("multiproduct_predictions_evaluation.csv","x")
output_file.write("num_prod,num_pred,cov1,fp1,cov2,fp2,cov3,fp3,cov4,fp4,cov5,fp5\n")

for index, target in enumerate(targets):
    output_string = ""
    coverage_fractions = [0,0,0,0,0]
    false_positive_fractions = [0,0,0,0,0]
    target_untokenized = "".join(target.split())
    individual_targets = get_individual_smiles_lists(target_untokenized,is_predictions=False)
    tokenized_prediction_list = predictions[5*index:5*(index+1)]
    predictions_untokenized = ".".join(["".join(prediction.split()) for prediction in tokenized_prediction_list])
    individual_predictions_as_lists = get_individual_smiles_lists(predictions_untokenized,is_predictions=True)

    individual_predictions = []
    for k in range(5):
        individual_predictions.extend(individual_predictions_as_lists[k])
        coverage_fractions[k] = get_coverage_fraction(individual_predictions,individual_targets)
        false_positive_fractions[k] = get_false_positive_fraction(individual_predictions,individual_targets)

    output_string = str(len(individual_targets)) + "," + str(len(individual_predictions))
    for k in range(5):
        output_string += "," + str(coverage_fractions[k]) + "," + str(false_positive_fractions[k])
    
    output_file.write(output_string + "\n")
    if index % 100 == 0:
        print("At index " + str(index))

output_file.close()

