#!/usr/bin/env python

"""
Simple bedtools / bash scripts wrapper with the following operations:
* union           multi argument union operation
* intersection    multi arguments intersection
* minus           remove all the ranges from firs file, which intersect with
                  second one
* compare         compares 2 files producing _cond1.bed, _cond2.bed and
                  _common.bed files
* metapeaks       compares multiple files and creates Venn diagram in case of
                  2 or 3 files
* process_pvalue  for each range show it's predecessors and compute min of
                  p_values (-log 10p)

NOTE: it is not supposed to replace pybedtools, but add some missing
functionality.

Major differences are:
 * union and intersection operations are commutative and associative
 * TRACE available for all the operations, i.e. consider operation
        minus(minus(intersect(A1, A2), intersect(B1, B2)), C)
   you would like to find out what were the predecessors of result and sort
   result by p-value.


NOTE: python3 required
> source activate py3.5

author oleg.shpynov@jetbrains.com
"""
import os
import subprocess
import tempfile
import atexit

from matplotlib_venn import venn2
from matplotlib_venn import venn3
from pathlib import Path

from scripts.util import run

UNION_SH = os.path.dirname(os.path.abspath(__file__)) + '/union.sh'
INTERSECT_SH = os.path.dirname(os.path.abspath(__file__)) + '/intersect.sh'
MINUS_SH = os.path.dirname(os.path.abspath(__file__)) + '/minus.sh'
COMPARE_SH = os.path.dirname(os.path.abspath(__file__)) + '/compare.sh'
METAPEAKS_SH = os.path.dirname(os.path.abspath(__file__)) + '/metapeaks.sh'
JACCARD_SH = os.path.dirname(os.path.abspath(__file__)) + '/jaccard.sh'

TEMPFILES = []


def columns(path):
    stdout, _stderr = run([
        ['grep', 'chr', path], ['head', '-1'], ['awk', '{ print NF }']
    ])
    return int(stdout.decode('utf-8').strip())


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
        stdout, _stderr = run([["cat", self.path], ['wc', '-l']])
        return int(stdout.strip())

    def save(self, path):
        self.compute()
        print(subprocess.Popen(["cp", self.path, path]).communicate())

    def save3(self, path):
        """ Save as BED3 format """
        self.compute()
        with open(path, mode='w') as out:
            run([['awk', "-v", "OFS=\\t", '{print $1,$2,$3}']],
                stdin=open(self.path), stdout=out)

    def collect_beds(self):
        return [self]

    def cat(self):
        stdout, _stderr = run([['cat', self.path]])
        return stdout.decode('utf-8')

    def head(self, lines=5):
        print('HEAD')
        stdout, _stderr = run([['head', '-{}'.format(lines), self.path]])
        print(stdout.decode('utf-8'))

    def tail(self, lines=5):
        print('TAIL')
        stdout, _stderr = run([['tail', '-{}'.format(lines), self.path]])
        print(stdout.decode('utf-8'))

    def pvalue_position(self):
        # TODO[shpynov]: add peak caller option as in diffbind?
        path = self.path
        if path.endswith('.broadPeak') or path.endswith('.narrowPeak'):
            return 8
        if 'diffbind' in path:
            return 10
        if 'bdgdiff' in path:
            return 5
        raise Exception('PValue column is unknown for {}'.format(path))

    def process_pvalue(self):
        """
        Method to process each row of resulting bed with the union of parents.
        """
        # Ensure that we've already computed result
        self.compute()
        c = columns(self.path)
        beds = self.collect_beds()
        with tempfile.NamedTemporaryFile(mode='w', suffix='_trace.bed',
                                         prefix='bedtrace',
                                         delete=False) as tmpfile:
            run([['bedtools', 'intersect', '-wa', '-wb', '-a', self.path,
                  '-b', *[x.path for x in beds],
                  '-names', *[str(x) for x in beds], '-sorted']],
                stdout=tmpfile)
            TEMPFILES.append(tmpfile.name)

            # Extract only pvalue_associated columns
            filtered_path = tmpfile.name.replace('_trace.bed', '_filtered.bed')
            TEMPFILES.append(filtered_path)
            with open(filtered_path, mode='a') as filtered_file:
                for bed in beds:
                    pvalue_position = bed.pvalue_position()
                    if len(beds) > 1:
                        # Names are available for multiple -b intersection
                        #  files only
                        run([['grep', str(bed)],
                             ['awk', "-v", "OFS=\\t",
                              '{{print $1,$2,$3,${}}}'
                              ''.format(c + pvalue_position + 1)]],
                            stdin=open(tmpfile.name), stdout=filtered_file)
                    else:
                        run([['awk', "-v", "OFS=\\t",
                              '{{print $1,$2,$3,${}}}'
                              ''.format(c + pvalue_position)]],
                            stdin=open(tmpfile.name), stdout=filtered_file)

            # Finally merge pvalues with min function
            result_path = filtered_path.replace('_filtered.bed', '_pvalue.bed')
            TEMPFILES.append(result_path)
            with open(result_path, 'w') as result_file:
                run([['sort', '-k1,1', '-k2,2n'],
                     ['bedtools', 'merge', '-c', '4', '-o', 'min'],
                     ['sort', '-k4,4gr']],
                    stdin=open(filtered_path), stdout=result_file)
            print('TRACE', tmpfile.name)
            print('RESULT', result_path)
            result_bed = Bed(result_path)
            result_bed.head()
            result_bed.tail()
            os.remove(filtered_path)
            return result_path


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
        return '\t' * indent + self.operation + '\n'\
               + '\n'.join([x.pp(indent + 1) for x in self.operands])

    def collect_beds(self):
        result = []
        for o in self.operands:
            result += o.collect_beds()
        return result


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
            raise Exception("Illegal {}: {}".format(self.operation,
                                                    str(self.operands)))
        self.path = self.intersect_files(*[x.path for x in self.operands])

    @staticmethod
    def intersect_files(*files):
        with tempfile.NamedTemporaryFile(
                mode='w', suffix='.bed', prefix='bedtraces', delete=False
        ) as tmpfile:
            run([["bash", INTERSECT_SH, *files]], stdout=tmpfile)
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
        self.path = self.minus_files(self.operands[0].path,
                                     self.operands[1].path)

    @staticmethod
    def minus_files(file1, file2):
        with tempfile.NamedTemporaryFile(
                mode='w', suffix='.bed', prefix='bedtraces', delete=False
        ) as tmpfile:
            run([["bash", MINUS_SH, file1, file2]], stdout=tmpfile)
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
            raise Exception("Illegal {}: {}".format(self.operation,
                                                    str(self.operands)))
        self.path = self.union_files(*[x.path for x in self.operands])

    @staticmethod
    def union_files(*files):
        with tempfile.NamedTemporaryFile(
                mode='w', suffix='.bed', prefix='bedtraces', delete=False
        ) as tmpfile:
            run([["bash", UNION_SH, *files]], stdout=tmpfile)
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
        with tempfile.NamedTemporaryFile(
                mode='w', suffix='.txt', prefix='bedtraces', delete=False
        ) as tmpfile:
            prefix = tmpfile.name.replace('.txt', '')
            run([["bash", COMPARE_SH, file1, file2, prefix]], stdout=tmpfile)
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


def jaccard(file1, file2):
    stdout, _stderr = run([['bash', JACCARD_SH, file1, file2]])
    return float(stdout)


def metapeaks(filesmap):
    """Plot venn diagrams for 2 or 3 files"""
    VENN2_PATTERNS = ["0 1", "1 0", "1 1"]
    VENN3_PATTERNS = ["0 0 1", "0 1 0", "0 1 1", "1 0 0", "1 0 1", "1 1 0",
                      "1 1 1"]

    def showvenn2(s1, s2, aB, Ab, AB):
        venn2(subsets=(Ab, aB, AB), set_labels=(s1, s2))

    def showvenn3(s1, s2, s3, abC, aBc, aBC, Abc, AbC, ABc, ABC):
        venn3(subsets=(Abc, aBc, ABc, abC, AbC, aBC, ABC),
              set_labels=(s1, s2, s3))

    if not isinstance(filesmap, dict):
        raise Exception("Map <name: bed> is expected")
    args = {}
    venn_patterns = None
    if len(filesmap) == 2:
        venn_patterns = VENN2_PATTERNS
    elif len(filesmap) == 3:
        venn_patterns = VENN3_PATTERNS
    else:
        print("Cannot create Venn diagram, wrong number of files",
              len(filesmap))

    names = filesmap.keys()
    # Check everything is computed
    for x in filesmap.values():
        x.compute()
    stdout, _stderr = run([['bash', METAPEAKS_SH,
                            *[filesmap[x].path for x in names]]])
    out = stdout.decode("utf-8")
    if venn_patterns:
        # Configure args for Venn diagram
        for p in venn_patterns:
            args[p] = 0
        for line in out.split('\n'):
            for p in venn_patterns:
                if p in line:
                    try:
                        args[p] = int(line[len(p):])
                    except:
                        pass
        if len(filesmap) == 2:
            showvenn2(*names, *[args[x] for x in venn_patterns])
        elif len(filesmap) == 3:
            showvenn3(*names, *[args[x] for x in venn_patterns])
    else:
        print(out)


def _cleanup():
    for path in TEMPFILES:
        if Path(path).is_file():
            os.remove(path)


atexit.register(_cleanup())
