import re
import math
from future.utils import implements_iterator


@implements_iterator
class ColourCovgsReader(object):

    """
    Iterable over coverages in colour_covgs file
    """

    def __init__(self, file):
        self._reader = file

    def __iter__(self):
        return self

    def __next__(self):
        colour_covgs = []
        for i in range(4):
            line = next(self._reader)
            colour_covgs.append(line)
        read = ColourCovgsRead(colour_covgs)
        return read


class ColourCovgsRead(object):

    """
    Colour Coverage read object
    """

    def __init__(self, lines):
        self.lines = lines
        self.lines[3] = self.lines[3].replace(
            "  ",
            " ").strip()  # for mccortex format
        self.seq = self.lines[1].replace('\n', '')[3:-3]
        self.name = lines[0].replace('\n', '')[1:]
        self.covgs = [
            int(i) for i in lines[3].replace(
                '\n',
                '').split(' ')[
                3:-
                3] if i]

    def __str__(self):
        return "%s: %i coverage on %i%% of seq" % (
            self.name, self.median_non_zero_coverage, self.percent_non_zero_coverage)

    def median(self, lis):
        if not lis:
            return None
        else:
            sorts = sorted(lis)
            length = len(lis)
            if length == 1:
                return sorts[0]
            if not length % 2:
                return (
                    sorts[int(length / 2)] + sorts[int(length / 2 - 1)]) / 2.0
            return sorts[int(length / 2)]

    def mean(self, lis):
        if len(lis) == 0:
            return None
        else:
            return float(sum(lis)) / len(lis)

    @property
    def non_zero_covgs(self):
        return [i for i in self.covgs if i > 0]

    @property
    def percent_non_zero_coverage(self):
        return int(100 *
                   float(sum([i > 0 for i in self.covgs])) /
                   len(self.covgs))

    @property
    def mean_coverage(self):
        return self.mean(self.covgs)

    @property
    def mean_non_zero_coverage(self):
        if self.non_zero_covgs:
            return self.mean(self.non_zero_covgs)
        else:
            return 0

    @property
    def median_coverage(self):
        return self.median(self.covgs)

    @property
    def median_non_zero_coverage(self):
        if self.non_zero_covgs:
            return int(math.floor(self.median(self.non_zero_covgs)))
        else:
            return 0
