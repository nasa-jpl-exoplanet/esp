#! /usr/bin/env python
'''intentionally cause the pipeline retrograde amnesia
'''

import argparse
import os

try:
    import dawgie.db.post
except ImportError:
    print('ERROR: it is best to be using a venv that has dawgie installed')
    print('       by installing esp/requirements.txt')

def cli():
    '''Command Line Interface'''
    cl = argparse.ArgumentParser(description=__doc__)
    cl.add_argument ('--dry-run', action='store_true',
                     help='do not do the actual deletions')
    cl.add_argument ('--nodes', default=[], nargs='*',
                     help='remove given nodes where the node is given in the form task.alg.sv.val. For instance, system.finalize.parameters.priors would remove all versions of priors. Shortening the name will remove all versions of everything under it as well. For instance, system.finalize would remove all versions of finalize and all of the state vectors and values below it.')
    cl.add_argument ('--targets', default=[], nargs='*',
                     help='remove given targets from memory')
    cl.add_argument ('--unique', action='store_true',
                     help='remove any duplicate rows from the prime table keeping those with the largest primary key')
    args = cl.parse_args()
    if not args.nodes and not args.targets and not args.unique:
        print ('WARN: you did not as me to do anything')
    else:
        print (args)
        confirm()
        connection = None
        cursor = None
        try:
            connection = dawgie.db.post._conn()
            cursor = dawgie.db.post._cur(connection)
            nodes(cursor, args.dry_run, args.nodes)
            target(cursor, args.dry_run, args.targets)
            if args.unique: unique(cursor, args.dry_run)
        except:
            if cursor is not None:
                cursor.rollback()
                cursor.close()
            if connection is not None:
                connection.close()
            raise

def confirm():
    '''check that some of the dawgie vars are defined in the environment'''
    passed = True
    for subname in ['DB_HOST', 'DB_IMPL', 'DB_NAME', 'DB_PATH', 'DB_PORT']:
        varname = 'DAWGIE_' + subname
        if varname not in os.environ:
            print (f'ERROR: {varname} is not in your environment')
            passed = False
    if not passed:
        raise AttributeError('Missing environment variables')
    if os.environ['DAWGIE_DB_IMPL'] != 'post':
        raise ValueError('DAWGIE_DB_IMPL must be post')
    pass

def nodes (cursor, dry, todo:[str]):
    '''process all of the nodes'''
    pass

def targets (cursor, dry, todo:[str]):
    '''process all of the targets'''
    for t in todo:
        print (f'INFO: removing target {t}')
    cursor.execute ('SELECT pk FROM Target WHERE name = ANY(%s);', todo)
    tnids = cursor.fetchall()

    if len(tnids) != len(todo):
        missing = todo.copy()
        cursor.execute ('SELECT name FROM Target WHERE pk = ANY(%s);', tnids)
        for name in cursor.fetchall():
            missing.remove (name)
        for m in missing:
            print (f'INFO: {m} is not in the database')

    if tnids not dry:
        cursor.execute('DELETE FROM Target WHERE pk = ANY(%s);')
        cursor.execute('DELETE FROM Prime WHERE tn_ID = ANY(%s);')
        cursor.commit()
    pass

def unique(cursor, dry):
    '''remove duplicate rows in the primary table keeping the largest PK'''
    cursor.execute ('SELECT pk,COUNT(*) FROM Prime GROUP BY '
                    'run_ID, tn_ID, task_ID, alg_ID, sv_ID, val_ID '
                    'HAVING COUNT(*) > 1;')
    duplicates = cursor.fetchall()
    total = sum(d[1]-1 for d in duplicates)
    print (f'INFO: found {len(duplicates)} duplicate groups and a total of {total} extra entries(rows)')
    if not dry:
        for group,_count in duplicates:
            group.sort()
            cursor.execute('DELETE FROM Prime WHERE pk = ANY(%s);' group[:-1])
        cursor.commit()
    pass

if __name__ == '__main__':
    cli()
