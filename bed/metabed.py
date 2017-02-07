#!/usr/bin/env python

"""
Simple bedtools / bash scripts wrapper with the following operations:
* union                     multi argument union operation
* intersection              multi arguments intersection
* minus                     remove all the ranges from first file, which intersect with second one
* compare                   compares 2 files producing _cond1.bed, _cond2.bed and _common.bed files
* metapeaks                 compares multiple files and creates Venn diagram in case of 2 or 3 files
* process_pvalue            for each range show it's predecessors and compute min of p_values (-log 10p)

NOTE: python3 required
> source activate py3.5

author oleg.shpynov@jetbrains.com
"""
import os
import subprocess
import tempfile
from pathlib import Path
from matplotlib_venn import venn2
from matplotlib_venn import venn3

UNION_SH = os.path.dirname(os.path.abspath(__file__)) + '/union.sh'
INTERSECT_SH = os.path.dirname(os.path.abspath(__file__)) + '/intersect.sh'
MINUS_SH = os.path.dirname(os.path.abspath(__file__)) + '/minus.sh'
COMPARE_SH = os.path.dirname(os.path.abspath(__file__)) + '/compare.sh'
METAPEAKS_SH = os.path.dirname(os.path.abspath(__file__)) + '/metapeaks.sh'

TEMPFILES = []


class Bed:
    """Simple path of Bed file storage"""

    def __init__(self, path):
        self.path = path

    def compute(self):
        if not Path(self.path).is_file():
            raise Exception("File not found: {}".format(self.path))

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

    def collect_beds(self):
        return [self]

    def pvalue_position(self):
        if self.path.endswith('.broadPeak') or self.path.endswith('.narrowPeak'):
            return 8
        # Ignore position
        return 100


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

    def collect_beds(self):
        result = []
        for o in self.operands:
            result += o.collect_beds()
        return result

    def process_pvalue(self):
        """Method to process each row of resulting bed with the union of parents."""
        # Ensure that we've already computed result
        self.compute()
        p1 = subprocess.Popen(['head', '-1', self.path], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['awk', '{ print NF }'], stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        columns = int(p2.communicate()[0].decode('utf-8').strip())
        beds = self.collect_beds()
        with tempfile.NamedTemporaryFile(mode='w', suffix='_trace.bed', prefix='metabed', delete=False) as tmpfile:
            subprocess.Popen(['bedtools', 'intersect', '-wa', '-wb', '-a', self.path,
                              '-b', *[x.path for x in beds],
                              '-names', *[str(x) for x in beds], '-sorted'],
                             shell=False, universal_newlines=True, stdout=tmpfile).communicate()
            TEMPFILES.append(tmpfile.name)

            # Extract only pvalue_associated columns
            filtered_path = tmpfile.name.replace('_trace.bed', '_filtered.bed')
            TEMPFILES.append(filtered_path)
            with open(filtered_path, mode='a') as filtered_file:
                for bed in beds:
                    pvalue_position = bed.pvalue_position()
                    p1 = subprocess.Popen(['grep', str(bed)], stdin=open(tmpfile.name), stdout=subprocess.PIPE)
                    p2 = subprocess.Popen(['awk', "-v", "OFS=\\t",
                                           '{print $1,$2,$3,$' + str(columns + pvalue_position + 1) + '}'],
                                          stdin=p1.stdout, stdout=filtered_file)
                    p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
                    p2.communicate()

            # Finally merge pvalues with min function
            result_path = filtered_path.replace('_filtered.bed', '_pvalue.bed')
            TEMPFILES.append(result_path)
            with open(result_path, 'w') as result_file:
                p1 = subprocess.Popen(['sort', '-k1,1', '-k2,2n'],
                                      stdin=open(filtered_path), stdout=subprocess.PIPE)
                p2 = subprocess.Popen(['bedtools', 'merge', '-c', '4', '-o', 'min'],
                                      stdin=p1.stdout, stdout=subprocess.PIPE)
                p3 = subprocess.Popen(['sort', '-k4,4nr'],
                                      stdin=p2.stdout, stdout=result_file)
                p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
                p2.stdout.close()  # Allow p2 to receive a SIGPIPE if p3 exits.
                p3.communicate()
            print('Trace file', tmpfile.name)
            print('Pvalue filtered file', filtered_path)
            print('Result file', result_path)
            return result_path


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

    def collect_beds(self):
        return self.operands[0].collect_beds()


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


def metapeaks(filesmap):
    """Plot venn diagrams for 2 or 3 files"""
    VENN2_PATTERNS = ["0 1", "1 0", "1 1"]
    VENN3_PATTERNS = ["0 0 1", "0 1 0", "0 1 1", "1 0 0", "1 0 1", "1 1 0", "1 1 1"]

    def showvenn2(s1, s2, aB, Ab, AB):
        venn2(subsets=(Ab, aB, AB), set_labels=(s1, s2))

    def showvenn3(s1, s2, s3, abC, aBc, aBC, Abc, AbC, ABc, ABC):
        venn3(subsets=(Abc, aBc, ABc, abC, AbC, aBC, ABC), set_labels=(s1, s2, s3))

    if not isinstance(filesmap, dict):
        raise Exception("Map <name: bed> is expected")
    args = {}
    if len(filesmap) == 2:
        patterns = VENN2_PATTERNS
    elif len(filesmap) == 3:
        patterns = VENN3_PATTERNS
    else:
        raise Exception("Wrong number of files", len(filesmap))

    # Configure args for Venn diagram
    for p in patterns:
        args[p] = 0
    names = filesmap.keys()
    ps = subprocess.Popen(['bash', METAPEAKS_SH, *[filesmap[x].path for x in names]], stdout=subprocess.PIPE)
    out = ps.communicate()[0].decode("utf-8")
    for line in out.split('\n'):
        for p in patterns:
            if p in line:
                try:
                    args[p] = int(line[len(p):])
                except:
                    pass
    if len(filesmap) == 2:
        showvenn2(*names, *[args[x] for x in patterns])
    elif len(filesmap) == 3:
        showvenn3(*names, *[args[x] for x in patterns])


def cleanup():
    for path in TEMPFILES:
        if Path(path).is_file():
            os.remove(path)
