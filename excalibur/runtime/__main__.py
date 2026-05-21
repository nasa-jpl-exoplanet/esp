'''test the code without a pipeline'''

import dawgie.db
import dawgie.security
import excalibur.runtime.algorithms
import excalibur.runtime.bot
import os

from excalibur.util.main import main_start

if 'EXCALIBUR_PRIVATE_PIPELINE_INDEPENDENT' in os.environ:
    # need to fake some dawgie stuff so that can run the unit independent of any
    # pipeline since it is the data in the configuration file that matters not any
    # data or state of the pipeline itself. There going to

    class FakeDawgie:
        def __init__(self):
            self.__name = 'runtime'

        def _bot(self):
            return self

        def _runid(self):
            return self

        def name(self):
            return self.__name

        def ds(self):
            return self

        def update(self):
            return self

    def connect(_alg, _bit, _tn):
        return FakeDawgie()

    def targets():
        return ['a', 'b', 'c', 'GJ 3193', 'GJ 3193 (taurex sim @TS)']

    setattr(dawgie.db, 'connect', connect)
    setattr(dawgie.db, 'targets', targets)
    test = excalibur.runtime.algorithms.Create()
    test.run(FakeDawgie())
    pass
elif __name__ == "__main__":
    rid, tn = main_start()

    if tn in ['', '__all__']:
        NAME = 'create'
        subtasks = excalibur.runtime.bot.AnalysisTeam('runtime', 4, rid)
        pass
    else:
        NAME = 'autofill'
        subtasks = excalibur.runtime.bot.TaskTeam('runtime', 4, rid, tn)
        pass
    subtasks.do(NAME)
    dawgie.db.close()
    dawgie.security.finalize()
    pass
