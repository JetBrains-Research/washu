#!/usr/bin/env python

"""
Simple bedtools / bash scripts wrapper with the following operations:
* commutative and associative UNION and INTERSECTION
* minus
* compare

NOTE: python3 required
> source activate py3.5

author oleg.shpynov@jetbrains.com
"""
import os
import subprocess
import tempfile
from pathlib import Path

UNION_SH = os.path.dirname(os.path.abspath(__file__)) + '/union.sh'
INTERSECT_SH = os.path.dirname(os.path.abspath(__file__)) + '/intersect.sh'
MINUS_SH = os.path.dirname(os.path.abspath(__file__)) + '/minus.sh'
COMPARE_SH = os.path.dirname(os.path.abspath(__file__)) + '/compare.sh'

TEMPFILES = []


class Bed:
    """Simple path of Bed file storage"""

    def __init__(self, path):
        self.path = path

    def compute(self):
        pass

    def __str__(self):
        return self.pp(0)

    def pp(self, indent):
        return '\t' * indent + os.path.basename(self.path)

    def count(self):
        self.compute()
        p1 = subprocess.Popen(["cat", self.path], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE)
        return int(p2.communicate()[0].strip())

    def save(self, path):
        self.compute()
        print(subprocess.Popen(["cp", self.path, path]).communicate())


class Operation(Bed):
    """Represents operations over Bed files and other Operations"""

    def __init__(self, operation=None, operands=None):
        super().__init__(None)
        self.operation = operation
        self.operands = operands

    def compute(self):
        raise Exception("Unknown operation: {}".format(self.operation))

    def __str__(self):
        return self.pp(0)

    def pp(self, indent):
        return '\t' * indent + self.operation + '\n' + '\n'.join([x.pp(indent + 1) for x in self.operands])


class Intersection(Operation):
    def __init__(self, operands):
        super().__init__("intersection", operands)

    def compute(self):
        # Do not compute twice
        if self.path is not None:
            return

        # Compute all the operands recursively
        for o in self.operands:
            o.compute()
        if len(self.operands) == 0:
            raise Exception("Illegal {}: {}".format(self.operation, str(self.operands)))
        self.path = self.intersect_files(*[x.path for x in self.operands])

    @staticmethod
    def intersect_files(*files):
        cmd = ["bash", INTERSECT_SH, *files]
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', prefix='metabeds', delete=False) as tmpfile:
            subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=tmpfile).communicate()
            TEMPFILES.append(tmpfile.name)
            return tmpfile.name


def intersect(*operands):
    return Intersection(operands)


class Minus(Operation):
    def __init__(self, operands):
        super().__init__("minus", operands)

    def compute(self):
        # Do not compute twice
        if self.path is not None:
            return

        # Compute all the operands recursively
        for o in self.operands:
            o.compute()
        if len(self.operands) != 2:
            raise Exception("Illegal minus: {}".format(str(self.operands)))
        self.path = self.minus_files(self.operands[0].path, self.operands[1].path)

    @staticmethod
    def minus_files(file1, file2):
        cmd = ["bash", MINUS_SH, file1, file2]
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', prefix='metabeds', delete=False) as tmpfile:
            subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=tmpfile).communicate()
            TEMPFILES.append(tmpfile.name)
            return tmpfile.name


def minus(*operands):
    return Minus(operands)


class Union(Operation):
    def __init__(self, operands):
        super().__init__("union", operands)

    def compute(self):
        # Do not compute twice
        if self.path is not None:
            return

        # Compute all the operands recursively
        for o in self.operands:
            o.compute()
        if len(self.operands) == 0:
            raise Exception("Illegal {}: {}".format(self.operation, str(self.operands)))
        self.path = self.union_files(*[x.path for x in self.operands])

    @staticmethod
    def union_files(*files):
        cmd = ["bash", UNION_SH, *files]
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', prefix='metabeds', delete=False) as tmpfile:
            subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=tmpfile).communicate()
            TEMPFILES.append(tmpfile.name)
            return tmpfile.name


def union(*operands):
    return Union(operands)


class Compare(Operation):
    def __init__(self, operands):
        super().__init__("compare", operands)
        self.cond1 = self.cond2 = self.common = None

    def compute(self):
        # Do not compute twice
        if self.path is not None:
            return

        if len(self.operands) != 2:
            raise Exception("Illegal compare: {}".format(str(self.operands)))
        self.path = self.compare(self.operands[0].path, self.operands[1].path)

    def compare(self, file1, file2):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', prefix='metabeds', delete=False) as tmpfile:
            prefix = tmpfile.name.replace('.txt', '')
            cmd = ["bash", COMPARE_SH, file1, file2, prefix]
            subprocess.Popen(cmd, shell=False, universal_newlines=True, stdout=tmpfile).communicate()

            self.cond1 = prefix + "_cond1.bed"
            self.cond2 = prefix + "_cond2.bed"
            self.common = prefix + "_common.bed"
            files = [self.cond1, self.cond2, self.common]
            Path(tmpfile.name).write_text('\n'.join(files))
            TEMPFILES.append(tmpfile.name)
            for f in files:
                TEMPFILES.append(f)
            return tmpfile.name


def compare(*operands):
    return Compare(operands)


def cleanup():
    for path in TEMPFILES:
        os.remove(path)
