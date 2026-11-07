"""Check SDK selection and compiler/wheel deployment-target agreement."""
import os
from pathlib import Path
import runpy
import unittest
from unittest.mock import patch

options = runpy.run_path(str(Path(__file__).resolve().parents[1] / 'build_support.py'))['macos_options']


class MacBuildTests(unittest.TestCase):
    @patch('sys.platform', 'darwin')
    @patch.dict(os.environ, {}, clear=True)
    @patch('pathlib.Path.is_dir', return_value=True)
    @patch('subprocess.check_output', side_effect=['13.7.2\n', '/SDK/MacOSX.sdk\n'])
    def test_actual_system_target_and_sdk(self, run, exists):
        flags, links, target = options(['-O3'], [])
        self.assertEqual(target, '11.0')
        self.assertIn('-mmacosx-version-min=11.0', flags)
        self.assertIn('-Wl,-syslibroot,/SDK/MacOSX.sdk', links)
        self.assertEqual(run.call_count, 2)

    @patch('sys.platform', 'darwin')
    @patch.dict(os.environ, {'MACOSX_DEPLOYMENT_TARGET': '13.0', 'SDKROOT': '/explicit-sdk'}, clear=True)
    @patch('pathlib.Path.is_dir', return_value=True)
    @patch('subprocess.check_output')
    def test_environment_overrides(self, run, exists):
        flags, links, target = options([], [])
        self.assertEqual((flags, links, target),
                         (['-mmacosx-version-min=13.0'], ['-Wl,-syslibroot,/explicit-sdk'], '13.0'))
        run.assert_not_called()

    @patch('sys.platform', 'darwin')
    @patch.dict(os.environ, {'MACOSX_DEPLOYMENT_TARGET': '15.0'}, clear=True)
    @patch('subprocess.check_output')
    def test_explicit_flags_take_precedence(self, run):
        flags = ['-O2', '-mmacosx-version-min=13.0']
        links = ['-Wl,-syslibroot,/user-sdk']
        self.assertEqual(options(flags, links), (flags, links, '13.0'))
        run.assert_not_called()

    @patch('sys.platform', 'darwin')
    @patch.dict(os.environ, {'MACOSX_DEPLOYMENT_TARGET': '13.0', 'SDKROOT': '/missing'}, clear=True)
    @patch('pathlib.Path.is_dir', return_value=False)
    def test_missing_sdk_fails_clearly(self, exists):
        with self.assertRaisesRegex(RuntimeError, 'SDK directory'):
            options([], [])

    @patch('sys.platform', 'linux')
    @patch('subprocess.check_output')
    def test_other_platforms_unchanged(self, run):
        self.assertEqual(options(['-O3'], ['-L/example']), (['-O3'], ['-L/example'], None))
        run.assert_not_called()


if __name__ == '__main__':
    unittest.main()
