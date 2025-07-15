'''target.monitor ds'''


# -- IMPORTS -- ------------------------------------------------------
from collections import defaultdict
import logging

# ------------- ------------------------------------------------------

log = logging.getLogger(__name__)

def create_identifier(frame_d):
    '''
    helper function that takes in a frame dict and
    returns the identifier, if it is HST or JWST
    '''
    if frame_d['observatory'] in ['HST', 'JWST']:
        identifier = (
            frame_d['observatory']
            + '-'
            + frame_d['instrument']
            + '-'
            + frame_d['detector']
            + '-'
        )
        # if filter is None: this frame was direct imaging
        if frame_d['filter'] is None:
            identifier += ''
        else:
            identifier += frame_d['filter'] + '-'
        if frame_d['mode'] is None:
            identifier += 'STARE'
        else:
            identifier += frame_d['mode']
        return identifier

def regress_for_frame_counts(
    data: {int: {str: int}},
    quality_dict: {int: {str: int, str: int}},
    tl: {str: {str: {str: object}}},
) -> ({str: float}, {str: []}, []):
    '''
    this functions regresses through the runids for target.scrape.databases.
    `data` currently contains:
    {observatory-instrument-detector-filter-mode: # of occurences}
    `quality` currently contains:
    {filter: {status: 1, -1}}
    '''



    for rid in tl:
        if rid in data:
            continue
        data[rid] = defaultdict(int)
        # only one product so loops only once
        for d in tl[rid].values():
            frames_d = d['name']
            for frame_d in frames_d.values():
                identifier = create_identifier(frame_d)
                if identifier:
                    data[rid][identifier] += 1

    # filter_name: [list of frame counts]
    history = defaultdict(list)

    # getting the most current frame counts and previous frames counts to compare
    sorted_keys = sorted(data.keys())

    quality_dict[sorted_keys[0]] = defaultdict(int)

    for filter_name in data[sorted_keys[0]]:
        history[filter_name].append(data[sorted_keys[0]][filter_name])
        quality_dict[sorted_keys[0]][filter_name] = 1

    for rid in sorted_keys[1:]:

        # if this rid already exists in quality_dict, skip
        if rid in quality_dict:
            continue

        quality_dict[rid] = defaultdict(int)

        for filter_name in data[rid]:
            # if this filter_name has already been seen, and frames exist
            # automatically status = 1
            if filter_name in history:
                history[filter_name].append(data[rid][filter_name])
                quality_dict[rid][filter_name] = 1
            # this is when a new filter_name is encountered
            else:
                # status automatically 1
                history[filter_name].append(data[rid][filter_name])
                quality_dict[rid][filter_name] = 1

        # set of filters that have been seen before but are not in this rid
        # these are dropped/not downloaded filters, so automatically -1
        dropped_filters = set(history).difference(set(data[rid]))
        for filter_name in dropped_filters:
            history[filter_name].append(0)
            quality_dict[rid][filter_name] = -1
