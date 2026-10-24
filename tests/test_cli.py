"""CLI regressions that run without optional analysis dependencies."""

import importlib.util
import json
from contextlib import redirect_stdout
from io import StringIO
from types import ModuleType
from unittest.mock import Mock
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from argparse import Namespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('genosyntax_cli', ROOT / 'main.py')
cli = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cli)


class CommandLineTests(unittest.TestCase):
    def test_help_is_available_without_analysis_imports(self):
        result = subprocess.run([sys.executable, str(ROOT / 'main.py'), '--help'],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.assertIn('analyze_slots', result.stdout)

    def test_missing_inputs_exit_nonzero(self):
        with tempfile.TemporaryDirectory() as folder:
            result = subprocess.run([
                sys.executable, str(ROOT / 'main.py'), 'run_stats',
                '--slot_analysis_dir', str(Path(folder) / 'missing'),
                '--output_dir', str(Path(folder) / 'output'), '--dry_run'],
                capture_output=True, text=True)
            self.assertEqual(result.returncode, 1)
            self.assertIn('File not found', result.stderr)
            self.assertFalse((Path(folder) / 'output').exists())

    def test_dry_run_checks_artifacts_without_creating_outputs(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            inputs = root / 'intermediate'
            inputs.mkdir()
            for name in ('slot_database.json', 'accessory_module_definitions.json'):
                (inputs / name).write_text('{}')
            result = subprocess.run([
                sys.executable, str(ROOT / 'main.py'), 'run_stats',
                '--slot_analysis_dir', folder, '--output_dir', str(root / 'output'),
                '--dry_run'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse((root / 'output').exists())

    def test_pipeline_preflight_requires_mafft(self):
        with tempfile.TemporaryDirectory() as folder:
            matrix = Path(folder) / 'matrix.csv'
            matrix.touch()
            args = Namespace(command='run', panaroo_file=str(matrix),
                             pharokka_dir=folder, fasta_dir=folder, faa_dir=folder)
            with patch('shutil.which', return_value=None):
                with self.assertRaisesRegex(ValueError, 'MAFFT'):
                    cli.validate_inputs(args)

    def test_no_anchors_stops_before_slot_analysis(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)
            architecture = output / '1_define_architecture'
            architecture.mkdir()
            (architecture / 'architecture_definition.json').write_text(
                json.dumps({'anchor_definitions': []}))
            modules = {}
            for name in ('genosyntax.define_architecture', 'genosyntax.slot_analysis',
                         'genosyntax.statistical_test', 'genosyntax.mechanism_analysis'):
                module = ModuleType(name)
                module.run = Mock(return_value=True)
                modules[name] = module
            package = ModuleType('genosyntax')
            package.summarizer = ModuleType('genosyntax.summarizer')
            package.summarizer.run = Mock(return_value=True)
            modules['genosyntax'] = package
            file_ops = ModuleType('genosyntax.utils.file_ops')
            file_ops.load_json = lambda path: json.loads(Path(path).read_text())
            modules['genosyntax.utils.file_ops'] = file_ops
            with patch.dict(sys.modules, modules), redirect_stdout(StringIO()):
                result = cli.run_pipeline(Namespace(output_dir=folder, until_step=3))
            self.assertTrue(result)
            modules['genosyntax.slot_analysis'].run.assert_not_called()
            modules['genosyntax.statistical_test'].run.assert_not_called()

    def test_zero_permutations_rejected_before_analysis(self):
        result = subprocess.run([
            sys.executable, str(ROOT / 'main.py'), 'run_stats',
            '--slot_analysis_dir', '.', '--output_dir', '.', '--permutations', '0'],
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn('--permutations must be positive', result.stderr)


if __name__ == '__main__':
    unittest.main()
