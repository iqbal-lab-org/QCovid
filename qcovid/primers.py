from collections import namedtuple, defaultdict

Primer = namedtuple(
    "Primer", ["amplicon", "name", "left", "forward", "rc", "pos", "length"]
)
Read = namedtuple("Read", ["name", "seq", "qual", "comment"])
Matched = namedtuple("Matched", ["r1", "p1", "r2", "p2"])


class ReadData:
    def __init__(self, name):
        self.amplicon = None
        self.good = False
        for data in name.split():
            d = data.split(":")
            if "qcovid-v0.1" in d:
                _, amplicon, good, primer_start = d
                if amplicon = "NO_MATCH":
                    self.amplicon = None
                else:
                    self.amplicon = amplicon
                self.good = bool(good)
                self.primer_start = int(primer_start)
                break

    def __str__(self):
        return f"qcovid-v0.1:{self.amplicon}:{self.good}:{self.primer_start}"


class Primers:
    def __init__(self, fn):
        self.name = fn
        self.min = 100
        self.max = 0
        self.seqs = {}
        _seqs = {}

        for l in open(fn):
            amplicon, name, seq, left, forward, pos = l.strip().split("\t")
            pos = int(pos)
            forward = forward.lower() in ["t", "+", "forward", "true"]
            left = left.lower() in ["left", "true", "t"]

            length = len(seq)
            if length > self.max:
                self.max = length
            if length < self.min:
                self.min = length

            _seqs[seq] = Primer(amplicon, name, left, forward, True, pos, length)
            _seqs[mp.revcomp(seq)] = Primer(
                amplicon, name, left, forward, False, pos, length
            )

        for k, v in _seqs.items():
            self.seqs[k[: self.min]] = v

    def match(self, seq):
        if seq[: self.min] in self.seqs:
            return self.seqs[seq[: self.min]]
        else:
            return None
