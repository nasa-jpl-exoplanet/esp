'''unit tests for the presumptive despot'''

import collections
import os
import time
import unittest

from excalibur.presumptive import fsm
from excalibur.presumptive import mia
from requests.exceptions import HTTPError
from unittest.mock import patch
from urllib.parse import urlparse

NS = collections.namedtuple(
    'NS', ['email', 'url', 'threshold', 'nodes', 'replicas']
)
ARGS = NS('t@x.y', 'https://excalibur.jpl.nasa.gov:8080', 3, [1, 2, 3], 5)
PERFORMED = [False]
RESPONSE = [{}]


def fn():
    return urlparse(ARGS.url).netloc


def mock_request_get(url, **_kwds):
    mock_response = unittest.mock.Mock()
    path = urlparse(url).path
    unittest.TestCase().assertTrue('timeoout' in kwds)
    unittest.TestCase().assertTrue('verify' in kwds)
    if path in ['/app/pl/state', '/app/schedule/crew']:
        mock_response.status_code = 200
        mock_response.json.return_value = RESPONSE[0]
        mock_response.raise_for_status.return_value = None
    else:
        mock_reponse.status_code = 404
        mock_response.json.return_value = None
        mock_response.raise_for_status.side_effect = HTTPError(
            f"404 Client Error: Not Found for url: {url}"
        )
    return mock_response


def mock_email_send(args, msg):
    print(msg)  # generally do not care but when errors can see the output


def mock_subprocess_run(*cmd, check=False, shell=False, **kwds):
    unittest.TestCase().assertTrue(check)
    unittest.TestCase().assertTrue(shell)
    unittest.TestCase().assertFalse(kwds)
    print('executed command:', ' '.join(cmd))
    PERFORMED[0] += 1


class ValidateDespotism(unittest.TestCase):
    def setUp(self):
        PERFORMED[0] = 0
        RESPONSE[0] = {}
        for section in ['fsm', 'mia']:
            if os.path.exists(f'/tmp/{fn()}.{section}.pkl'):
                os.unlink(f'/tmp/{fn()}.{section}.pkl')

    @patch(
        'excalibur.presumptive.fsm.requests.get', side_effect=mock_request_get
    )
    @patch('excalibur.presumptive.fsm.email.send', side_effect=mock_email_send)
    @patch(
        'excalibur.presumptive.fsm.perform.subprocess.run',
        side_effect=mock_subprocess_run,
    )
    def test_check_running(self, mock_run, mock_email, mock_get):
        RESPONSE[0] = {'name': 'running', 'status': 'active'}
        retval = fsm.check(ARGS)
        self.assertEqual(0, retval)
        self.assertFalse(os.path.exists(f'/tmp/{fn()}.fsm.pkl'))
        self.assertFalse(PERFORMED[0])

    @patch(
        'excalibur.presumptive.fsm.requests.get', side_effect=mock_request_get
    )
    @patch('excalibur.presumptive.fsm.email.send', side_effect=mock_email_send)
    @patch(
        'excalibur.presumptive.fsm.perform.subprocess.run',
        side_effect=mock_subprocess_run,
    )
    def test_check_inactive(self, mock_run, mock_email, mock_get):
        RESPONSE[0] = {'name': 'running', 'status': 'inactive'}
        retval = fsm.check(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.fsm.pkl'))
        self.assertFalse(PERFORMED[0])
        time.sleep(1)
        retval = fsm.check(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.fsm.pkl'))
        self.assertFalse(PERFORMED[0])
        time.sleep(2)
        retval = fsm.check(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.fsm.pkl'))
        self.assertTrue(PERFORMED[0])

    @patch(
        'excalibur.presumptive.fsm.requests.get', side_effect=mock_request_get
    )
    @patch('excalibur.presumptive.fsm.email.send', side_effect=mock_email_send)
    @patch(
        'excalibur.presumptive.fsm.perform.subprocess.run',
        side_effect=mock_subprocess_run,
    )
    def test_check_loading(self, mock_run, mock_email, mock_get):
        RESPONSE[0] = {'name': 'loading', 'status': 'active'}
        retval = fsm.check(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.fsm.pkl'))
        self.assertFalse(PERFORMED[0])
        time.sleep(1)
        retval = fsm.check(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.fsm.pkl'))
        self.assertFalse(PERFORMED[0])
        time.sleep(2)
        retval = fsm.check(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.fsm.pkl'))
        self.assertFalse(PERFORMED[0])

    @patch(
        'excalibur.presumptive.fsm.requests.get', side_effect=mock_request_get
    )
    @patch('excalibur.presumptive.fsm.email.send', side_effect=mock_email_send)
    @patch(
        'excalibur.presumptive.fsm.perform.subprocess.run',
        side_effect=mock_subprocess_run,
    )
    def test_check_state_change(self, mock_run, mock_email, mock_get):
        RESPONSE[0] = {'name': 'apple', 'status': 'active'}
        retval = fsm.check(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.fsm.pkl'))
        self.assertFalse(PERFORMED[0])
        time.sleep(1)
        RESPONSE[0] = {'name': 'banana', 'status': 'active'}
        retval = fsm.check(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.fsm.pkl'))
        self.assertFalse(PERFORMED[0])
        time.sleep(2)
        retval = fsm.check(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.fsm.pkl'))
        self.assertFalse(PERFORMED[0])
        time.sleep(2)
        retval = fsm.check(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.fsm.pkl'))
        self.assertTrue(PERFORMED[0])

    @patch(
        'excalibur.presumptive.mia.requests.get', side_effect=mock_request_get
    )
    @patch('excalibur.presumptive.mia.email.send', side_effect=mock_email_send)
    @patch(
        'excalibur.presumptive.mia.perform.subprocess.run',
        side_effect=mock_subprocess_run,
    )
    def test_worker_all_idle(self, mock_run, mock_email, mock_get):
        RESPONSE[0] = {'busy': [], 'idle': '15'}
        retval = mia.worker(ARGS)
        self.assertEqual(0, retval)
        self.assertFalse(os.path.exists(f'/tmp/{fn()}.mia.pkl'))
        self.assertFalse(PERFORMED[0])

    @patch(
        'excalibur.presumptive.mia.requests.get', side_effect=mock_request_get
    )
    @patch('excalibur.presumptive.mia.email.send', side_effect=mock_email_send)
    @patch(
        'excalibur.presumptive.mia.perform.subprocess.run',
        side_effect=mock_subprocess_run,
    )
    def test_worker_all_present(self, mock_run, mock_email, mock_get):
        RESPONSE[0] = {
            'busy': ['a b c', 'b b c', 'c b c', 'd b c'],
            'idle': '11',
        }
        retval = mia.worker(ARGS)
        self.assertEqual(0, retval)
        self.assertFalse(os.path.exists(f'/tmp/{fn()}.mia.pkl'))
        self.assertFalse(PERFORMED[0])

    @patch(
        'excalibur.presumptive.mia.requests.get', side_effect=mock_request_get
    )
    @patch('excalibur.presumptive.mia.email.send', side_effect=mock_email_send)
    @patch(
        'excalibur.presumptive.mia.perform.subprocess.run',
        side_effect=mock_subprocess_run,
    )
    def test_worker_fallen(self, mock_run, mock_email, mock_get):
        RESPONSE[0] = {
            'busy': ['a b c', 'b b c', 'c b c'],
            'idle': '11',
        }
        retval = mia.worker(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.mia.pkl'))
        self.assertFalse(PERFORMED[0])
        retval = mia.worker(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.mia.pkl'))
        self.assertTrue(PERFORMED[0])

    @patch(
        'excalibur.presumptive.mia.requests.get', side_effect=mock_request_get
    )
    @patch('excalibur.presumptive.mia.email.send', side_effect=mock_email_send)
    @patch(
        'excalibur.presumptive.mia.perform.subprocess.run',
        side_effect=mock_subprocess_run,
    )
    def test_worker_fallen_with_change(self, mock_run, mock_email, mock_get):
        RESPONSE[0] = {
            'busy': ['a b c', 'b b c', 'c b c'],
            'idle': '11',
        }
        retval = mia.worker(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.mia.pkl'))
        self.assertFalse(PERFORMED[0])
        RESPONSE[0] = {
            'busy': ['a b c', 'b b c'],
            'idle': '11',
        }
        retval = mia.worker(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.mia.pkl'))
        self.assertFalse(PERFORMED[0])
        retval = mia.worker(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.mia.pkl'))
        self.assertTrue(PERFORMED[0])

    @patch(
        'excalibur.presumptive.mia.requests.get', side_effect=mock_request_get
    )
    @patch('excalibur.presumptive.mia.email.send', side_effect=mock_email_send)
    @patch(
        'excalibur.presumptive.mia.perform.subprocess.run',
        side_effect=mock_subprocess_run,
    )
    def test_worker_undead(self, mock_run, mock_email, mock_get):
        RESPONSE[0] = {
            'busy': ['a b c', 'b b c', 'c b c'],
            'idle': '15',
        }
        retval = mia.worker(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.mia.pkl'))
        self.assertFalse(PERFORMED[0])
        retval = mia.worker(ARGS)
        self.assertEqual(1, retval)
        self.assertTrue(os.path.exists(f'/tmp/{fn()}.mia.pkl'))
        self.assertTrue(PERFORMED[0])
