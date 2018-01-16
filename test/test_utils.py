import os
import subprocess

import pytest
from test.fixtures import test_data, tmp_dir

import scripts.util as su
from pipeline_utils import PROJECT_ROOT_PATH, run


def test_lcs():
    lcs1 = su.lcs('UW_CD14_RO01746_k4me3_1_1_ENCFF001FYS.fastq',
                  'Broad_CD14_2_input_ENCFF000CCW.fastq')
    lcs2 = su.lcs('UW_CD14_RO01746_k4me3_1_1_ENCFF001FYS.fastq',
                  'UW_CD14_input_ENCFF001HUV.fastq')
    assert lcs1 < lcs2


@pytest.mark.parametrize("input,donor", [
    ("", "40_donor6_input.bam"),
    ("40_donor6_input.bam", "37_donor6_k27ac.bam"),
    ("44_DONOR7_INPUT.bam", "41_donor7_k27ac.bam"),
    ("jcl320_wt1_gm_input_ctrl.1919_8.R1_mm10.bam",
     "jcl320_ko_gm_h3k27ac_chipd_dna.1919_8.R1_mm10.bam"),
    ("jcl320_wt1_gm_input_ctrl.1919_8.R1_mm10.bam",
     "jcl320_wt1_gm_bhlhe40_chipd_dna.1919_8.R1_mm10.bam"),
])
def test_find_input(test_data, input, donor):
    assert su.find_input(test_data("input/" + donor)) == input


@pytest.mark.parametrize("build,species", [
    ("hg18", "hs"),
    ("hg19", "hs"),
    ("hg37", "hs"),
    ("mm9", "mm"),
    ("mm10", "mm"),
])
def test_macs_input(build, species):
    assert su.macs_species(build) == species


@pytest.mark.parametrize("path,expected", [
    ("geo", "geo"),
    ("geo/.", "geo"),
    ("geo/tmp/..", "geo"),
    ("geo/../geo/tmp/doo/..", "geo/tmp"),
    ("geo/../geo/tmp/doo/file.txt", "geo/tmp/doo/file.txt"),
    ("geo/symlink", "geo/tmp"),
    ("geo/symlink/..", "geo"),
    ("geo/../geo/symlink/doo/file.txt", "geo/tmp/doo/file.txt"),
])
def test_expand_path(tmp_dir, capfd, path, expected):
    os.makedirs(os.path.join(tmp_dir, "geo/tmp/doo"))
    os.symlink(os.path.join(tmp_dir, "geo/tmp"),
               os.path.join(tmp_dir, "geo/symlink"))
    open(os.path.join(tmp_dir, "geo/tmp/doo/file.txt"), 'a').close()

    path = os.path.join(tmp_dir, path)
    util_sh = os.path.join(PROJECT_ROOT_PATH, "parallel/util/util.sh")
    subprocess.run(
        "bash -c 'source {}; echo $(expand_path \"{}\")'".format(util_sh,
                                                                 path),
        shell=True, check=True)
    out, _err = capfd.readouterr()
    assert out == os.path.join(tmp_dir, expected) + "\n"


def test_module_mock(tmp_dir, capfd):
    with open(os.path.join(tmp_dir, "modules.sh"), 'a') as f:
        f.write("module() { echo \"my test mock module $@\"; }\n")

    with open(os.path.join(tmp_dir, "foo.sh"), 'a') as f:
        f.write("source {}/modules.sh\n"
                "source {}/parallel/util/util.sh\n"
                "module load R\n".format(tmp_dir, PROJECT_ROOT_PATH))
    run("bash", "{}/foo.sh".format(tmp_dir))
    out, _err = capfd.readouterr()

    assert out.replace(tmp_dir, ".") == "bash ./foo.sh\n" \
                                        "my test mock module load R\n"


def test_run_parallel(tmp_dir, capfd):
    with open(os.path.join(tmp_dir, "foo.sh"), 'a') as f:
        f.write("""
source {}/parallel/util/util.sh
TASKS=""
for i in $(seq 1 100); do
    echo $i
    run_parallel << SCRIPT
#It is necessary to include LOG here because run_parallel use it as out/err
#PBS -o {}/file_$i.log
echo $i > {}/file_$i.txt
SCRIPT
    TASKS="$TASKS $QSUB_ID"
done
wait_complete $TASKS
""".format(PROJECT_ROOT_PATH, tmp_dir, tmp_dir))

    run("bash", "{}/foo.sh".format(tmp_dir))
    out, _err = capfd.readouterr()

    # Check expected stdout result
    assert out.replace(tmp_dir, ".") == "bash ./foo.sh\n" + "\n".join(
        [str(i) for i in range(1, 101)]
    ) + "\nLOCAL waiting for tasks...\nDone. LOCAL waiting for tasks\n"

    # Check that files and logs created successfully
    for i in range(1, 101):
        assert os.path.isfile('{}/file_{}.txt'.format(tmp_dir, i))
        assert os.path.isfile('{}/file_{}.log'.format(tmp_dir, i))


def test_run_stderr(capfd):
    stdout, stderr = su.run([["bash", "-c", "echo 'error!' 1>&2"]])

    assert stdout == b""
    assert stderr == b"error!\n"  # should be: b"error!\n"

    out, err = capfd.readouterr()
    assert out == ""
    assert err == ""


def test_run_stderr_piped(capfd):
    stdout, stderr = su.run([["echo", "ok1"],
                             ["bash", "-c", "echo \"error!\" 1>&2"],
                             ["echo", "ok2"]])

    out, err = capfd.readouterr()

    assert stdout == b"ok2\n"
    assert stderr == b""  # should be: b"error!\n"
    # assert stderr is None  # should be: b"error!\n"
    assert out == ""
    assert err == "error!\n"


def test_run_error_piped(capfd):
    stdout, stderr = su.run([["echo", "ok1"],
                             ["bash", "-c", "exit 1"],
                             ["echo", "ok2"]])

    out, err = capfd.readouterr()

    assert stdout == b"ok2\n"
    assert stderr == b""  # should be: b"error!\n"
    # assert stderr is None  # should be: b"error!\n"
    assert out == ""
    assert err == ""


def test_check_logs(tmp_dir, capfd):
    with open(os.path.join(tmp_dir, "macs2.log"), 'a') as f:
        f.write("""
ValueError: cannot resize this array: it does not own its data
Exception ValueError: 'cannot resize this array: it does not own its data' in ...
""")
    with open(os.path.join(tmp_dir, "foo.sh"), 'a') as f:
        f.write("cd {}\n"
                "source {}/parallel/util/util.sh\n"
                "check_logs\n".format(tmp_dir, PROJECT_ROOT_PATH))
    run("bash", "{}/foo.sh".format(tmp_dir))
    out, _err = capfd.readouterr()
    assert _err.replace(tmp_dir, ".") == "Local tasks WASHU_PARALLELISM=8\n"
    assert out.replace(tmp_dir, ".") == "bash ./foo.sh\n"
