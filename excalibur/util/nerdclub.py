'''
GMR: For Nerds
'''

# IMPORTS
import sys
import time
import datetime
import logging
import numpy as np

logger = logging.getLogger(__name__)


class Progressbar:
    """
    0.0.0:2024/09/20:GMR:Progress bar nerd club
    0.1.0:2024/09/21:GMR:Added timing in the log
    0.2.0:2024/09/21:GMR:Added nerdy colors + visual
    0.2.1:2024/09/25:GMR:End of loop fix for small progress bars
    0.2.2:2024/10/24:GMR:Fixed rounding error condition
    """

    def __init__(self, argsdict, title, iterobj):  # pylint: disable=R0902
        """
        argsdict['progbar']: True / False
        argsdict['progsizemax']: Max size of the progress bar in [command window prompts]
        title: Progress bar title, should be str-able
        iterobj: Iterative object, should be len-able
        """
        self._argsdict = argsdict
        self.title = title
        self.iterobj = iterobj
        self.progbsize = len(iterobj)
        self.maxindex = len(iterobj)
        self.scale = 1.0
        self.current = 0
        self.call = 0
        self.start = time.time()
        self.stop = None
        self.done = False
        self.percent = int(100 * self.call * self.scale / self.progbsize)
        self.colors = [
            "\033[1;31m",
            "\033[1;35m",
            "\033[1;33m",
            "\033[1;32m",
            "\033[1;36m",
            "\033[1;34m",
        ]
        if self.progbsize > self._argsdict["progsizemax"]:
            self.scale = 1e0 * self._argsdict["progsizemax"] / self.progbsize
            self.progbsize = self._argsdict["progsizemax"]
            pass
        if len(self.title) < self._argsdict["lbllen"]:
            self.title = self.title + " " * (
                self._argsdict["lbllen"] - len(self.title)
            )
            pass
        if len(self.title) > self._argsdict["lbllen"]:
            self.title = self.title[0 : self._argsdict["lbllen"] - 3] + "..."
            pass
        if self._argsdict["progbar"]:
            if self._argsdict["proginprompt"]:
                sys.stdout.write(f"{self.title} [{' ' * self.progbsize}]")
                sys.stdout.flush()
                sys.stdout.write("\b" * (self.progbsize + 1))
                pass
            else:
                print(f"\t{self.title}")
            pass
        pass

    def update(self):
        """update"""
        if self._argsdict["progbar"]:
            self.call += 1
            dotick = self.call * self.scale > self.current
            if dotick:
                if self._argsdict["proginprompt"]:
                    index = int(np.floor(self.percent / 100 * len(self.colors)))
                    sys.stdout.write(f"{self.colors[index]}")
                    sys.stdout.write("-")
                    sys.stdout.write("\033[0;0m")
                    sys.stdout.flush()
                    pass
                self.current += 1
                self.percent = int(
                    100 * self.call * self.scale / self.progbsize
                )
                pass
            if (not self.done) and (self.call == self.maxindex):
                if self._argsdict["proginprompt"]:
                    sys.stdout.write("]")
                    sys.stdout.flush()
                    self.done = True
                    pass
                pass
            pass
        return

    def close(self):
        """close"""
        self.stop = time.time()
        ftime = datetime.timedelta(seconds=self.stop - self.start)
        logger.info("\t\t\t%s: %s", self.title, ftime)
        if self._argsdict["progbar"]:
            if self._argsdict["proginprompt"]:
                fillme = " " * (self._argsdict["progsizemax"] - self.progbsize)
                sys.stdout.write(f" {fillme}{ftime}\n")
                pass
            pass
        return

    pass
