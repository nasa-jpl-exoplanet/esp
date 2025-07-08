'''target.monitor ds'''

# -- IMPORTS -- ------------------------------------------------------
from collections import defaultdict

# ------------- ------------------------------------------------------


def regress_for_frame_counts(
    data: {int: {str: int, str: int}},
    quality_dict: {int: int},
    tl: {str: {str: {str: object}}},
) -> ({str: float}, {str: []}, []):
    '''
    this functions regresses through the runids for target.scrape.databases.
    data currently contains:

    observatory-instrument-detector-filter-mode: # of occurences

    timeline iterates through descending order of runID. Therefore
    the newest runIDs go first, then code breaks once it hits a runID
    that has already been computed.
    #'''

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
                        + frame_d['filter']
                        + '-'
                    )
                    if frame_d['mode'] is None:
                        identifier += 'STARE'
                    else:
                        identifier += frame_d['mode']
                    data[rid][identifier] += 1

    # getting the most current frame counts and previous frames counts to compare
    sorted_keys = sorted(data.keys())
    # this if is redundant, will remove later
    if len(sorted_keys) > 1:
        # i know im traversing the rids again, i feel like there's a good reason
        # to do this but if i find out there isn't ill try to conslidate the loops
        for i, rid in enumerate(sorted_keys[1:], start=1):

            # if this rid already exists in quality_dict, skip
            if rid in quality_dict:
                print('key alr exists, continuing')
                continue
            # add quality flag entry for rid
            quality_dict[rid] = 1
            cur_frames = data[rid]
            prev_frames = data[sorted_keys[i - 1]]

            for identifier in prev_frames:
                if identifier not in cur_frames:
                    # checking if previously this key existed (ie, had some frames) but in current doesn't?
                    quality_dict[rid] = -1
                    break

    return data, quality_dict
