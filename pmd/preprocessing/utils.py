import gzip
import shlex
import subprocess

def execute_command(command : str, *, stdin : any = None,
                    capture : bool = True) -> subprocess.CompletedProcess:
    """
    Execute a shell command with arguments and optional stdin value.

    Args:
        stdin    any: Values to pipe to the stdin.
        capture bool: Whether to capture the stdout and stderr.

    Returns: subprocess.CompletedProcess
    """

    cmdlist = shlex.split(command)
    if stdin is not None:
        stdin = str(stdin)
        stdin = bytes(stdin, encoding='utf-8')
    result = subprocess.run(cmdlist, shell=False,
                            capture_output=capture, input=stdin)
    return result


def extract_gzip(input_gz):
    """ Extract a gzip file. Returns the extracted file path. """
    if not input_gz.endswith(".gz"):
        raise ValueError("not a gzip file", input_gz)
    outfile = input_gz[:-3]

    # First extract the gzip file.
    fp = gzip.open(input_gz, 'rb')
    file_content = fp.read().decode("utf-8")
    fp.close()

    # Save
    with open(outfile, "w+") as fp:
        fp.write(file_content)
    return outfile
