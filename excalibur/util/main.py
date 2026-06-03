'''util main ds'''

# -- IMPORTS --------------------------------------------------------
import os
import dawgie
import dawgie.db
import dawgie.pl.scan
import dawgie.security
import dawgie.util

# -------------------------------------------------------------------


def main_start():
    '''
    code that goes at the start of all __main__.py files
    '''

    fep = os.environ.get('FE_PORT', None)
    rid = int(os.environ.get('RUNID', None))
    tn = os.environ.get('TARGET_NAME', None)

    if fep:
        dawgie.util.set_ports(int(fep))

    dawgie.security.initialize(
        os.path.expandvars(
            os.path.expanduser(dawgie.context.guest_public_keys)
        ),
        myname=dawgie.context.ssl_pem_myname,
        myself=dawgie.context.ssl_pem_myself,
        system=dawgie.context.ssl_pem_file,
    )
    dawgie.pl.scan.for_factories(
        dawgie.context.ae_base_path, dawgie.context.ae_base_package
    )
    dawgie.db.reopen()
    return rid, tn
