'''target.monitor ds'''

# -- IMPORTS -- ------------------------------------------------------
from collections import defaultdict
import logging
# ------------- ------------------------------------------------------

log = logging.getLogger(__name__)


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
                    data[rid][identifier] += 1

    # filter: [list of frame counts]
    history = defaultdict(list)

    # getting the most current frame counts and previous frames counts to compare
    sorted_keys = sorted(data.keys())

    quality_dict[sorted_keys[0]] = defaultdict(int)

    # unique_filter_names = set(data[sorted_keys[0]])
    for filter in data[sorted_keys[0]]:
        history[filter].append(data[sorted_keys[0]][filter])
        quality_dict[sorted_keys[0]][filter] = 1

    # this if is redundant, will remove later
    if len(sorted_keys) > 1:
        # i know im traversing the rids again, i feel like there's a good reason
        # to do this but if i find out there isn't ill try to conslidate the loops
        for rid in sorted_keys[1:]:

            # if this rid already exists in quality_dict, skip
            if rid in quality_dict:
                continue

            quality_dict[rid] = defaultdict(int)

            for filter in data[rid]:
                # if this filter has already been seen, and frames exist
                # automatically status = 1
                if filter in history:
                    history[filter].append(data[rid][filter])
                    quality_dict[rid][filter] = 1
                # this is when a new filter is encountered
                else:
                    # status automatically 1
                    history[filter].append(data[rid][filter])
                    quality_dict[rid][filter] = 1

            # set of filters that have been seen before but are not in this rid
            # these are dropped/not downloaded filters, so automatically -1
            dropped_filters = set(history).difference(set(data[rid]))
            for filter in dropped_filters:
                history[filter].append(0)
                quality_dict[rid][filter] = -1
