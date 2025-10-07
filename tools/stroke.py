#! /usr/bin/env python
'''--- intentionally cause the pipeline retrograde amnesia ---

The prime table is brain of the pipeline connecting all of the nodes, targets,
and run ids to specific data. Removing references from the prime table means
disassociating that run 12 for target GJ 1214 with '
system.finalize.parameters.priors' from the data that was computed at that time.
Bassically, the pipeline strokes causing retrograde amnesia. Just like strokes
in us, there can be hidden and complicated side effects that accumulate over
time causing systemic and/or catastrophic failure far into the future.

Other database tables mentioned help in minimizing the prime table size and add
to the potential side effects.


There are three forms of retrograde amnesia:
  1. remove targets
  2. remove nodes - task.algorithm.state_vactor.value
  3. remove duplicates - only useful on databases prior to dawgie 2.0.0

1. Removing targets is straight forward. If a pipeline contained the target
   names 'snafu' and 'foobar', then removing them would be done with
   `stroke.py --targets snafu foobar`. These names will be removed from the
   known target table and all references to them in the prime table will also
   be removed.

2. Removing nodes - computational elements like the bubbles in the dependency
   graph - is more complicated because there is a natural globbing (wildcarding)
   that takes place as well. Let us assume the pipeline contains four nodes:
   system.finalize.parameters.priors, system.finalize.parameters.autofill,
   system.finalize.madeup, and system.validate. Independent of their length,
   they are all nodes. The take should be that fewer '.' means more information
   is removed. So now lets look at what happens when we provide these nodes:
 a. system.finalize.parameters.priors
    All versions of system.finalize.parameters.priors are removed from the
    value table and any references to them in the prime table.
 b. system.finalize.parameters
    It globs like `ls` would glob 'system.finalize.parameters.*'. All version of
    system.finalize.parameters.priors and system.finalize.parameters.autofill
    are removed from the values table and those references from the prime table.
    All versions of system.finalize.parmeters are removed from state vector
    table with those references being removed from the prime table.
 c. system.finalize
    All versions of 'autofill' and 'priors' are removed from the value table
    with those references being remvoed from the prime table. All versions of
    'parameters' and 'madeup' are removed from the state vector table with those
    references being removed from the prime table. Finally, all versions of
    system.finalize are removed from the algorithm table with those references
    being removed from the prime table.
 d. system
    Everything is removed as stated in (c) and all versions of system,validate.
    'system' is removed from the task table with all of those references from
    the prime table.

3. Make the prime table unique. It is required for dawgie < 2 before upgrading.
   Otherwise this option is superfluous.
'''

import argparse
import os

try:
    import dawgie.db.post
except ImportError:
    print('ERROR: it is best to be using a venv that has dawgie installed')
    print('       by installing esp/requirements.txt')


def _resolve(cursor, name, parents, table):
    if name and parents[0]:
        cursor.execute(
            f'SELECT pk FROM {table} WHERE name = %s AND {parents[0]} = ANY(%s);',
            [name, parents[1]],
        )
    elif name:
        cursor.execute(f'SELECT pk FROM {table} WHERE name = %s;', [name])
    elif parents[1]:
        cursor.execute(
            f'SELECT pk FROM {table} WHERE {parents[0]} = ANY(%s);',
            [parents[1]],
        )
    else:
        ValueError('name and parents were None or empty')

    pks = cursor.fetchall()
    if pks:
        return [pk[0] for pk in pks]
    raise AttributeError(f'{name} could not be found in {table}')


def cli():
    '''Command Line Interface'''
    cl = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    cl.add_argument(
        '--dry-run',
        action='store_true',
        help='do not actually delete the information',
    )
    cl.add_argument(
        '--nodes', default=[], nargs='*', help='remove given nodes.'
    )
    cl.add_argument(
        '--targets',
        default=[],
        nargs='*',
        help='remove given targets from memory',
    )
    cl.add_argument(
        '--unique',
        action='store_true',
        help='remove any duplicate rows from the prime table keeping those with the largest primary key',
    )
    args = cl.parse_args()
    if not args.nodes and not args.targets and not args.unique:
        print('WARN: you did not as me to do anything')
    else:
        confirm()
        connection = None
        cursor = None
        try:
            connection = dawgie.db.post._conn()
            cursor = dawgie.db.post._cur(connection)
            nodes(cursor, args.dry_run, args.nodes)
            targets(cursor, args.dry_run, args.targets)
            if args.unique:
                unique(cursor, args.dry_run)
        finally:
            if cursor is not None:
                cursor.close()
            if connection is not None:
                connection.rollback()
                connection.close()


def confirm():
    '''check that some of the dawgie vars are defined in the environment'''
    passed = True
    for subname in ['DB_HOST', 'DB_IMPL', 'DB_NAME', 'DB_PATH', 'DB_PORT']:
        varname = 'DAWGIE_' + subname
        if varname not in os.environ:
            print(f'ERROR: {varname} is not in your environment')
            passed = False
    if not passed:
        raise AttributeError('Missing environment variables')
    if os.environ['DAWGIE_DB_IMPL'] != 'post':
        raise ValueError('DAWGIE_DB_IMPL must be post')
    pass


def nodes(cursor, dry, todo: [str]):
    '''process all of the nodes'''
    for node in todo:
        tn, an, svn, vn = (node.split('.') + [None, None, None])[:4]
        tids, aids, svids, vids = [], [], [], []
        parents = ['', []]
        for name, parent, pids, table in zip(
            (tn, an, svn, vn),
            ('', 'task_ID', 'alg_ID', 'sv_ID'),
            (tids, aids, svids, vids),
            ('Task', 'Algorithm', 'StateVector', 'Value'),
        ):
            parents[0] = parent
            pids.extend(_resolve(cursor, name, parents, table))
            parents[1] = pids
        cursor.execute(
            'SELECT pk FROM Prime WHERE '
            'task_ID = ANY(%s) AND alg_ID = ANY(%s) AND '
            'sv_ID = ANY(%s) AND val_ID = ANY(%s);',
            [tids, aids, svids, vids],
        )
        pks = [pk[0] for pk in cursor.fetchall()]
        print(f'INFO: For node {node}:')
        print(
            f'        From {tn} removing'
            if an
            else f'        Removing {len(tids)} task: {tn}'
        )
        print(
            f'        Form {an} removing'
            if svn
            else (
                f'        Removing {len(aids)} {an}(s)'
                if an
                else f'        Removing {len(aids)} algorithms'
            )
        )
        print(
            f'        Form {svn} removing'
            if vn
            else (
                f'        Removing {len(svids)} {svn}(s)'
                if svn
                else f'        Removing {len(svids)} state vectors'
            )
        )
        print(
            f'        Removing {len(vids)} {vn}(s)'
            if an
            else f'        Removing {len(vids)} values'
        )
        print(f'        Removing from PRIME {len(pks)} associated rows')

        if pks and not dry:
            sql = 'DELETE FROM {} WHERE pk = ANY(%s);'
            cursor.execute(sql.format('Prime'), [pks])
            cursor.execute(sql.format('Value'), [vids])
            if not vn:
                cursor.execute(sql.format('StateVector'), [svids])
            if not svn:
                cursor.execute(sql.format('Algorithm'), [aids])
            if not an:
                cursor.execute(sql.format('Task'), [tids])
    pass


def targets(cursor, dry, todo: [str]):
    '''process all of the targets'''
    for t in todo:
        print(f'INFO: removing target {t}')
    cursor.execute('SELECT pk FROM Target WHERE name = ANY(%s);', [todo])
    tnids = [pk[0] for pk in cursor.fetchall()]
    cursor.execute('SELECT pk from Prime WHERE tn_ID = ANY(%s);', [tnids])
    pks = cursor.fetchall()
    if pks:
        print(f'INFO: targets remove {len(pks)} rows from Prime')
    if len(tnids) != len(todo):
        missing = todo.copy()
        cursor.execute('SELECT name FROM Target WHERE pk = ANY(%s);', [tnids])
        for name, in cursor.fetchall():
            missing.remove(name)
        for m in missing:
            print(f'INFO: {m} is not in the database')

    if tnids and not dry:
        cursor.execute('DELETE FROM Prime WHERE tn_ID = ANY(%s);', [tnids])
        cursor.execute('DELETE FROM Target WHERE pk = ANY(%s);', [tnids])


def unique(cursor, dry):
    '''remove duplicate rows in the primary table keeping the largest PK'''
    cursor.execute(
        'SELECT COUNT(*) FROM Prime GROUP BY '
        'run_ID, tn_ID, task_ID, alg_ID, sv_ID, val_ID '
        'HAVING COUNT(*) > 1;'
    )
    basis = cursor.fetchall()
    total = sum(b[-1] - 1 for b in basis)
    print(
        f'INFO: found {len(basis)} duplicate groups and a total of {total} extra entries(rows)'
    )
    if total:
        cursor.execute(
            'SELECT pk, run_ID, task_ID, tn_ID, alg_ID, sv_ID, val_ID '
            'FROM ( '
            'SELECT pk, run_ID, task_ID, tn_ID, alg_ID, sv_ID, val_ID, '
            'COUNT(*) '
            'OVER (PARTITION BY run_ID, task_ID, tn_ID, alg_ID, sv_ID, val_ID) '
            'AS occurs '
            'FROM prime) '
            'AS counted_rows '
            'WHERE occurs > 1;'
        )
        rows = cursor.fetchall()
        duplicates = {}
        for row in rows:
            key = row[1:]
            if key not in duplicates:
                duplicates[key] = []
            duplicates[key].append(row[0])
        tot = sum(len(group) - 1 for group in duplicates.values())
        if len(duplicates) != len(basis):
            print(
                f'ERROR: expected {len(basis)} groupings but resolved {len(duplicates)}'
            )
            return
        if total != tot:
            print(f'ERROR: expected {total} duplicates but found {tot}')
            return
        if not dry:
            print('INFO:   Deleting rows from Prime')
            allpks = []
            for group in duplicates.values():
                group.sort()
                allpks.extend(group[:-1])
            cursor.execute('DELETE FROM Prime WHERE pk = ANY(%s);', [allpks])


if __name__ == '__main__':
    cli()
