'''unit tests for keeping pyxb and xsd up to date'''

import hashlib
import os
import unittest


class Pyxb(unittest.TestCase):
    def test_01(self):
        '''check that binding can be imported

        This will check that a recent enough version of pyxb is installed
        and that it matches the version that generated excalibur.runtime.binding
        '''
        try:
            import excalibur.runtime.binding

            self.assertTrue(True)
        except:
            self.assertTrue(False, 'problems with excalibur.runtime.binding')

    def test_02(self):
        '''check that the binding is newer than the schema'''
        basedir = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), '..', 'excalibur', 'runtime'
            )
        )
        with open(os.path.join(basedir, 'autogen.md5'), 'rt') as file:
            for line in file.readlines():
                cksum, fn = line.split()
                with open(os.path.join(basedir, fn), 'br') as f:
                    content = f.read()
                self.assertEqual(
                    cksum,
                    hashlib.md5(content).hexdigest(),
                    'please run tools/autogen.sh because checksums '
                    f'do not match for {fn}',
                )

    def test_03(self):
        '''check that excalibur._version_ is in good shape'''
        fn = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), '..', 'excalibur', '__init__.py'
            )
        )
        with open(fn, 'rt') as file:
            locs = file.readlines()
        found = "__version__ = '${ESP_GIT_REV}'\n" in locs
        self.assertTrue(
            found,
            "excalibur/__init__.py has wrong value for __version__. It needs to be `__version__ = '${ESP_GIT_REV}'`. Please correct and try again.",
        )
