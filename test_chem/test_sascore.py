from pathlib import Path, WindowsPath
from typing import Callable
from unittest import TestCase, main

from rdkit import Chem

from df.data_transfer import DataFunctionRequest, DataFunctionResponse
from df.sascorer import readFragmentScores, calculateScore
from ruse.rdkit.rdkit_utils import standardize_mol
from yaml_to_script import extract_script

PRINTED_GUFF = False


def run_script(in_file: str, execute: Callable[[DataFunctionRequest], DataFunctionResponse]) -> DataFunctionResponse:
    with open(in_file, 'r') as fh:
        request_json = fh.read()

    request = DataFunctionRequest.parse_raw(request_json)
    response = execute(request)
    return response


class ScriptTest(TestCase):

    def setUp(self) -> None:
        """
        It's important when doing production testing to use the script
        that's in the YAML, which is in the DataFxn repo, as that's the
        one that Spotfire will be running.  When developing, it's
        convenient to edit the script in this directory, transferring
        it to the YAML file once it's ready for production.
        If there's a script SyntheticAccessibilityScore_script_dev.py,
        use that, otherwise make one from the YAML file and use that.
        Arbitrary_script_dev.py should not be in the
        repo, only ever in a development branch so testing a production
        repo will use the script derived from the YAML.
        Doing it this way allows for people editing the YAML directly,
        which is handy for minor tweaks, and those tweaks always
        being tested.
        Assume the DataFxn repo has been cloned alongside this one.
        """
        global PRINTED_GUFF
        this_dir = Path(__file__).parent
        script_file = this_dir / 'SyntheticAccessibilityScore_script_dev.py'
        if Path(script_file).exists():
            if not PRINTED_GUFF:
                print(f'Using development script')
                PRINTED_GUFF = True
            from SyntheticAccessibilityScore_script_dev import execute
            self._script_file = None
        else:
            data_fxn_dir = this_dir.parent.parent / 'DataFxns'
            yaml_file = data_fxn_dir / 'python' / 'local' / 'SyntheticAccessibilityScore.yaml'
            if not PRINTED_GUFF:
                print(f'Using production script in {yaml_file}')
                PRINTED_GUFF = True
            self._script_file = this_dir / 'QEDScore_script.py'
            if yaml_file.exists():
                script_lines = extract_script(yaml_file)
                with open(self._script_file, 'w') as f:
                    f.write(''.join(script_lines))

            from SyntheticAccessibilityScore_script import execute
        self._execute = execute

    def tearDown(self) -> None:
        # if a script file was made, tidy it up
        if self._script_file is not None and self._script_file.exists():
            self._script_file.unlink()
            pass  # so we can comment the unlink() easily

    def test_script(self) -> None:
        file_in = Path(__file__).parent / 'resources' / 'test_sascore1.json'
        response = run_script(file_in, self._execute)
        scores = response.outputColumns[0].values
        self.assertEqual(15, len(scores))
        self.assertAlmostEqual(1.0, scores[0], 1)
        self.assertAlmostEqual(1.2, scores[-3], 1)
        self.assertTrue(response)

    def test_ertl_scores(self) -> None:
        ertl_file = Path(__file__).parent / 'resources' / 'ertl_test_structs.smi'
        suppl = Chem.SmilesMolSupplier(str(ertl_file), titleLine=False)
        exp_scores = {"Mol_1": 1.0000, "Mol_2": 2.8849, "Mol_3": 1.0000,
                      "Mol_4": 1.0874, "Mol_5": 1.0000, "Mol_6": 4.3597,
                      "Mol_7": 3.0356, "Mol_8": 1.0000, "Mol_9": 3.7411,
                      "Mol_10": 1.8176, "Mol_11": 4.4665, "Mol_12": 1.0000,
                      "Mol_13": 1.0000, "Mol_14": 1.0000, "Mol_15": 4.5446,
                      "Mol_16": 2.4960, "Mol_17": 3.8409, "Mol_18": 4.7418,
                      "Mol_19": 1.0000, "Mol_20": 1.2344, "Mol_21": 1.0000,
                      "Mol_22": 1.8705, "Mol_23": 1.0000, "Mol_24": 5.0739,
                      "Mol_25": 1.0000, "Mol_26": 1.0000, "Mol_27": 4.1296,
                      "Mol_28": 4.6406, "Mol_29": 2.2408, "Mol_30": 4.7943,
                      "Mol_31": 1.0000, "Mol_32": 3.1092, "Mol_33": 4.7043,
                      "Mol_34": 1.0000, "Mol_35": 1.1181, "Mol_36": 1.0000,
                      "Mol_37": 1.0000, "Mol_38": 2.1356, "Mol_39": 3.9957,
                      "Mol_40": 1.0000, }
        for mol in suppl:
            sm = standardize_mol(mol)
            score = calculateScore(sm)
            self.assertAlmostEqual(exp_scores[mol.GetProp("_Name")], score, 3)
            # print(f'"{mol.GetProp("_Name")}" : {score:.4f},')


if __name__ == '__main__':
    main()
