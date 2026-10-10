"""Unit tests for aaftf/fcs_screen.py (the FCS-adaptor container is mocked)."""

import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from aaftf.fcs_screen import run
from aaftf.resources import FCSADAPTOR

pytestmark = pytest.mark.unit


def _which(available):
    return lambda name: f"/usr/bin/{name}" if name in available else None


def _run(tmp_path, engine="docker", available=("docker",), returncode=0, write_clean=True, **kw):
    infile = tmp_path / "asm.fasta"
    infile.write_text(">a\nACGT\n")
    wd = tmp_path / "wd"
    out = tmp_path / "out.fasta"

    def fake_run_cmd(cmd, debug=False, **k):
        if write_clean:
            Path(wd, "cleaned_sequences").mkdir(parents=True)
            Path(wd, "cleaned_sequences", "asm.fasta").write_text(">a\nACG\n")
            Path(wd, "fcs_adaptor_report.txt").write_text("#report\na\t4\tTRIM\n")
        return subprocess.CompletedProcess(cmd, returncode)

    with (
        patch("aaftf.fcs_screen.shutil.which", side_effect=_which(available)),
        patch("aaftf.fcs_screen.run_cmd", side_effect=fake_run_cmd) as rc,
        patch("aaftf.fcs_screen.require_databases", side_effect=lambda names, hint=None: [f"/db/{n}" for n in names]) as req,
    ):
        run(infile=str(infile), outfile=str(out), container_engine=engine, workdir=str(wd), **kw)
    return rc.call_args.args[0], req, out


@pytest.mark.parametrize(("prok", "tax"), [(False, "--euk"), (True, "--prok")])
def test_docker_default_image_and_tax(tmp_path, prok, tax):
    cmd, req, out = _run(tmp_path, available=("docker", "run_fcsadaptor.sh"), prok=prok)
    image = FCSADAPTOR["DOCKERIMAGE"] % FCSADAPTOR["VERSION"]
    assert cmd == ["/usr/bin/run_fcsadaptor.sh", "--fasta-input", str(tmp_path / "asm.fasta"), "--output-dir", str(tmp_path / "wd"), tax, "--container-engine", "docker", "--image", image]
    assert req.call_args.args[0] == []


def test_singularity_uses_database_script_and_image(tmp_path):
    cmd, req, _ = _run(tmp_path, engine="singularity", available=("apptainer",))
    assert req.call_args.args[0] == ["fcs_script", "fcs_image"]
    assert cmd[0] == "/db/fcs_script"
    assert cmd[-1] == "/db/fcs_image"


def test_singularity_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="singularity"):
        _run(tmp_path, engine="singularity", available=(), fcs_script="x.sh", image="i.sif")


def test_docker_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="docker"):
        _run(tmp_path, available=(), fcs_script="x.sh")


def test_unknown_engine_raises(tmp_path):
    with pytest.raises(ValueError, match="podman"):
        _run(tmp_path, engine="podman", fcs_script="x.sh")


def test_nonzero_exit_raises(tmp_path):
    with pytest.raises(RuntimeError, match="exit 1"):
        _run(tmp_path, returncode=1, fcs_script="x.sh")


def test_missing_clean_file_raises(tmp_path):
    with pytest.raises(RuntimeError, match="cleaned_sequences"):
        _run(tmp_path, write_clean=False, fcs_script="x.sh")


def test_success_writes_outputs(tmp_path, capsys):
    cmd, _, out = _run(tmp_path, fcs_script="x.sh", image="my/image:1", debug=True)
    assert cmd[0] == "x.sh" and cmd[-1] == "my/image:1"
    assert out.read_text() == ">a\nACG\n"
    assert Path(str(out) + ".fcs_adaptor_report.txt").read_text() == "#report\na\t4\tTRIM\n"
    assert "a\t4\tTRIM" in capsys.readouterr().out
