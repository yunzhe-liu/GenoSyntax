"""Regressions for slot counting and the implemented permutation null model."""

from contextlib import redirect_stdout, redirect_stderr
import importlib.util
from io import StringIO
from pathlib import Path
import random
import unittest

ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SlotAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dynamics = load_module('dynamics', 'genosyntax/slot_analysis/dynamics_analyzer.py')
        cls.stats = load_module('permutations', 'genosyntax/statistical_test/permutation_test.py')

    def test_partial_module_is_not_counted(self):
        modules = [{'module_id': 'M', 'module_type': 'multi-gene', 'num_genes': 2,
                    'members': ['x', 'y'], 'consensus_function': 'unknown'}]
        slots = {'g1': {'s': {'content': ['x'], 'slot_id_low_res': 'S',
                             'slot_id_high_res': 'S'}}}
        with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
            result = self.dynamics.analyze_module_dynamics(modules, slots, {}, ['g1'])
        self.assertEqual(result.iloc[0]['occurrence_count'], 0)
        self.assertEqual(result.iloc[0]['spi'], 0)

    def test_empty_module_report_keeps_headers(self):
        with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
            result = self.dynamics.analyze_module_dynamics([], {}, {}, [])
        self.assertTrue(result.empty)
        self.assertIn('spi', result.columns)

    def test_member_counts_and_seeded_permutations(self):
        modules = [{'members': ['x', 'y'], 'num_genes': 2}]
        slots = {'g1': {'a': {'content': ['x', 'y'], 'slot_id_low_res': 'A'},
                        'b': {'content': ['core', 'core'], 'slot_id_low_res': 'B'}}}
        with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
            random.seed(7)
            all_genes, multi = self.stats.run_permutation_test(slots, modules, 100)
            random.seed(7)
            repeat, _ = self.stats.run_permutation_test(slots, modules, 100)
        self.assertTrue(all_genes.equals(repeat))
        self.assertEqual(multi.set_index('slot_type').loc['A', 'observed_module_gene_count'], 2)
        self.assertEqual(all_genes['total_genes_in_slot'].sum(), 4)
        self.assertTrue((all_genes['p_value_enrichment'] > 0).all())


if __name__ == '__main__':
    unittest.main()
