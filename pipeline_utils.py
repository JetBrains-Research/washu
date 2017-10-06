import os
import os.path
from os.path import abspath, realpath, dirname
import argparse
import subprocess
import itertools
import shutil

from glob import glob

PROJECT_ROOT_PATH = abspath(os.path.join(dirname(realpath(__file__))))


class WritableDirectory(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if not os.path.isdir(values):
            raise argparse.ArgumentTypeError(
                "{0} is not a valid path".format(values)
            )
        if os.access(values, os.R_OK + os.W_OK):
            setattr(namespace, self.dest, values)
        else:
            raise argparse.ArgumentTypeError(
                "{0} is not a writable directory".format(values)
            )


def run_bash(script, *params):
    """
    Executes desired bash script and wait until it has finished
    :param script: script relative (to project root) path
    :param params: script args
    """
    cmdline = " ".join(["bash", os.path.join(PROJECT_ROOT_PATH, script),
                        *[str(p) for p in params]])
    print(cmdline)
    subprocess.run(cmdline, shell=True, check=True)


def run(*params):
    cmdline = " ".join([str(p) for p in params])
    print(cmdline)
    subprocess.run(cmdline, shell=True)


def move_forward(folder, new_folder, what_to_move, copy_files=False):
    # Necessary for correct copy behavior
    os.chdir(folder)
    if not os.path.exists(new_folder):
        os.mkdir(new_folder)
    for pattern in itertools.chain.from_iterable(map(glob, what_to_move)):
        if copy_files:
            shutil.copy(pattern, new_folder)
        else:
            shutil.move(pattern, new_folder)


def file_path_type(dir=False, exists=True, ext=None):
    """
    Argparse argument converter.
    Usage:
        parser.add_argument("tsv_file", type=file_path_type(ext="tsv"))

    :param dir: Whether directory or file expected
    :param exists: Ensure file or directory exists
    :param ext: Desired file extension. If None - do not check
    :return: Converter lambda: (arg: str) -> absolute path string
    """
    def inner(arg: str) -> str:
        msg = None
        if exists and not os.path.exists(arg):
            msg = "File not exists: " + arg
        elif not dir and not os.path.isfile(arg):
            msg = "File expected, by was: " + arg
        elif dir and (os.path.exists(arg) and not os.path.isdir(arg)):
            msg = "Dir expected, by was: " + arg
        elif ext and not arg.lower().endswith(".{}".format(ext)):
            msg = "File *.{} expected, by was: {}".format(ext, arg)

        if msg:
            raise argparse.ArgumentTypeError(msg)
        return os.path.abspath(arg)
    return inner
