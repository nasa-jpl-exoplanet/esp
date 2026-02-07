'''utility to send emails to a given list of people (tested)'''

import os
import pickle

from email.message import EmailMessage
from smtplib import SMTP

FN = '/tmp/presumptive.despot.last.email.pickle'


def send(args, msg: str):
    if os.path.isfile(FN):
        with open(FN, 'br') as file:
            count, mask, previous = pickle.load(file)
    else:
        count = -1
        mask = 0x1
        previous = ''

    if msg != previous:
        count = -1
        mask = 0x1
        previous = msg

    count += 1
    if not count & mask:
        if count:
            mask = 0x7F
            msg = f'''

After {count} messages, it is time to remind you:

{msg}
            '''
        print(msg)  # log this to the diaray.wannabe.despot
        em = EmailMessage()
        em.set_content(msg)
        em['Subject'] = 'ALERT: excalibur pipeline trouble'
        em['From'] = 'presumptive.despot@excalibur.jpl.nasa.gov'
        em['To'] = args.emails
        with SMTP('localhost', 25) as server:
            server.starttls()
            server.send_message(em)

    with open(FN, 'bw') as file:
        pickle.dump((count, mask, previous), file)
