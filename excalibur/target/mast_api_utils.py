'''
MAST API Utils from https://mast.stsci.edu/api/v0/MastApiTutorial.html
'''

# -- IMPORTS -- ------------------------------------------------------
import sys
import json
import requests

from urllib.parse import quote as urlencode

import logging

log = logging.getLogger(__name__)


# ------------- ------------------------------------------------------
def mast_query(request, maxwaittime=42):
    """
    Perform a MAST query. Parameters
    ----------
    request (dictionary): The MAST request json object
    maxwaittime: how long before timeout (in seconds)
    Returns head,content where head is the response HTTP headers, and content is the returned data
    """
    # Base API url
    request_url = 'https://mast.stsci.edu/api/v0/invoke'
    # Grab Python Version
    version = ".".join(map(str, sys.version_info[:3]))
    # Create Http Header Variables
    headers = {
        "Content-type": "application/x-www-form-urlencoded",
        "Accept": "text/plain",
        "User-agent": "python-requests/" + version,
    }
    # Encoding the request as a json string
    req_string = json.dumps(request)
    req_string = urlencode(req_string)
    # Perform the HTTP request
    try:
        resp = requests.post(
            request_url,
            data="request=" + req_string,
            headers=headers,
            timeout=maxwaittime,
        )
        # Pull out the headers and response content
        head = resp.headers
        content = resp.content.decode('utf-8')

    except requests.exceptions.ReadTimeout:
        log.warning(
            '--< TIMEDOUT on the mast query (%s seconds)  >--', maxwaittime
        )
        head = None
        content = None
    return head, content


def set_filters(parameters):
    '''set_filters ds'''
    return [{"paramName": p, "values": v} for p, v in parameters.items()]


def set_min_max(thismin, thismax):
    '''set_min_max ds'''
    return [{'min': thismin, 'max': thismax}]
