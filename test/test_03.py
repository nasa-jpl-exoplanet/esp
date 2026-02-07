'''unit tests for presumptive despot emails'''

import collections
import excalibur.presumptive.email
import os
import unittest

from unittest.mock import patch

ARGS = collections.namedtuple('ARGS', ['emails'])
MSG_A = 'a simple test message'
MSG_B = 'be simple test message'


class MockSMTP:
    sent = False

    def __enter__(self):
        MockSMTP.sent = False
        return self

    def __exit__(self, *_args, **_kwds):
        pass

    def __init__(self, *_args, **_kwds):
        pass

    def starttls(self):
        pass

    def send_message(self, *_args, **_kwds):
        MockSMTP.sent = True


class TestSpamReduction(unittest.TestCase):
    @patch('excalibur.presumptive.email.SMTP', side_effect=MockSMTP)
    def test_first(self, mock_smtp):
        if os.path.isfile(excalibur.presumptive.email.FN):
            os.unlink(excalibur.presumptive.email.FN)
        excalibur.presumptive.email.send(ARGS('apple'), MSG_A)
        self.assertTrue(MockSMTP.sent)
        MockSMTP.sent = False
        excalibur.presumptive.email.send(ARGS('apple'), MSG_B)
        self.assertTrue(MockSMTP.sent)
        MockSMTP.sent = False
        excalibur.presumptive.email.send(ARGS('apple'), MSG_A)
        self.assertTrue(MockSMTP.sent)
        MockSMTP.sent = False
        excalibur.presumptive.email.send(ARGS('apple'), MSG_A)
        self.assertFalse(MockSMTP.sent)
        MockSMTP.sent = False
        excalibur.presumptive.email.send(ARGS('apple'), MSG_A)
        self.assertTrue(MockSMTP.sent)
        MockSMTP.sent = False
        for _i in range(125):
            excalibur.presumptive.email.send(ARGS('apple'), MSG_A)
            self.assertFalse(MockSMTP.sent)
        excalibur.presumptive.email.send(ARGS('apple'), MSG_A)
        self.assertTrue(MockSMTP.sent)
        MockSMTP.sent = False
