'''utility to send emails to a given list of people (tested)
'''

from email.message import EmailMessage
from smtplib import SMTP

def send(args, msg:str):
    em = EmailMessage()
    em.set_content(msg)
    em['Subject'] = 'ALERT: excalibur pipeline trouble'
    em['From'] = 'presumptive.despot@excalibur.jpl.nasa.gov'
    em['To'] = args.emails
    with SMTP('localhost', 25) as server:
        server.starttls()
        server.send_message(em)
