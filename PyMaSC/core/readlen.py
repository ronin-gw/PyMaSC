import logging
from collections import Counter

logger = logging.getLogger(__name__)


def _mean(c):
    return int(round(
        sum(length * freq for length, freq in c.items()) / float(sum(c.values()))
    ))


def _median(c):
    num = sum(c.values())
    target = num / 2
    _sum = 0

    if num % 2:
        for l in sorted(c):
            _sum += c[l]
            if target <= _sum:
                return l
    else:
        length = sorted(c)
        for i, l in enumerate(length):
            _sum += c[l]
            if target < _sum:
                return l
            elif target == _sum:
                return int(round((l + float(length[i+1])) / 2))


def _mode(c):
    return c.most_common(1)[0][0]


ESTFUNCTIONS = dict(
    MEAN=_mean, MEDIAN=_median, MODE=_mode, MIN=min, MAX=max
)


def estimate_readlen(parser, source, esttype, mapq_criteria):
    logger.info("Estimate read length: get {} from read length distribution...".format(esttype.lower()))
    estfunc = ESTFUNCTIONS[esttype]

    with parser(source) as stream:
        length = estfunc(Counter(
            readlen for _, is_duplicate, _, _, mapq, readlen in stream
            if not is_duplicate and mapq >= mapq_criteria
        ))

    logger.info("Estimated read length = {}".format(length))
    return length
