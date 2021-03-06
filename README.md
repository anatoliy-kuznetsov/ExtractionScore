Described in our manuscript [ExtractionScore: A Quantitative Framework for Evaluating Synthetic Routes on Predicted Liquid–Liquid Extraction Performance](https://pubs.acs.org/doi/full/10.1021/acs.jcim.0c01426).

## Installation
This project uses [rxn4chemistry](https://github.com/rxn4chemistry/rxn4chemistry) as well as [MolecularTransformer](https://github.com/pschwllr/MolecularTransformer). To start, place ``write_synthesis_routes.py`` in the ``rxn4chemistry`` directory and ``extractionscore.py`` and the evaluation scripts in the ``MolecularTransformer`` directory. These projects have conflicting dependencies, so you should use separate environments for them.

## Using the test set
To test Molecular Transformer on the test set, download the version you want to test from [here](https://ibm.box.com/v/MolecularTransformerModels) and use ``multiproduct_source_strict.txt`` as the source file. The option ``-n_best`` corresponds to _k_, the number of predictions to include; integer values between 1 and 5 can be specified. To evaluate predictions, use ``evaluate_multiproduct_predictions.py``.

## Generating retrosyntheses
Create a Python file named ``api_info.py`` in the ``rxn4chemistry`` directory, and place your IBM RXN ``API_KEY`` and ``PROJECT_ID`` there. Write the industrially practiced synthesis routes following the format of ``industrial_synthesis_routes.txt``, then run ``write_synthesis_routes.py``. If IBM RXN is successfully able to predict the major products in each synthetic step for a target molecule, the retrosynthesis will be generated and written (along with the industrially practiced route) to ``route_pairs.txt``.

## Calculating ExtractionScore
Copy ``route_pairs.txt`` into your ``MolecularTransformer`` directory and activate the appropriate Python environment. Then, run ``extractionscore.py``. Detailed output will be written to ``extractionscorelog.txt``. Stereoisomers of the intended product that do not match the recorded product will result in a manual prompt so that you can verify that no stereocenters have been incorrectly switched. The code to draw reactions is based heavily on a [script](https://github.com/connorcoley/retrosim/blob/master/retrosim/utils/draw.py) by Connor Coley. If the VCC Lab website is down, one workaround is to write a file with SMILES for all your species and use the ALOGPS implementation hosted by [OCHEM](https://ochem.eu/home/show.do) to estimate partition coefficients instead.

## Visualizing synthesis routes
A convenient way to visualize synthesis routes is with the free [ChemDraw Direct](https://chemdrawdirect.perkinelmer.cloud/js/sample/index.html#) online tool. Select "Structure > Load SMILES" from the ribbon and paste a SMILES string into the field. This is the method we have used for our work. Alternatively, one could modify the reaction drawing code to visualize synthesis routes.
