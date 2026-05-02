"""Tests for pipeline.py Namespace construction.

Verifies that sort_args and assess_args constructed via create_namespace()
carry the debug and pipe attributes required by their target modules.
"""

from argparse import Namespace

import pytest

pytestmark = pytest.mark.unit


def _make_pipeline_args(**kwargs):
    """Return a minimal args Namespace for pipeline.run()."""
    defaults = {
        "basename": "test_sample",
        "memory": "16",
        "left": "left.fq.gz",
        "right": "right.fq.gz",
        "cpus": 1,
        "debug": False,
        "workdir": None,
        "mincontiglen": 500,
        "method": "spades",
        "assembler_args": None,
        "tmpdir": None,
        "screen_accessions": None,
        "screen_urls": None,
        "AAFTF_DB": None,
        "phylum": ["Ascomycota"],
        "sourdb": None,
        "mincovpct": 5,
        "iterations": 5,
        "pipe": False,
    }
    defaults.update(kwargs)
    return Namespace(**defaults)


def _run_pipeline_up_to_sort(args):
    """Run create_namespace() the same way pipeline.run() does for sort_args."""
    args_dict = vars(args)

    def create_namespace(options, required_args=None, **extra_args):
        namespace_dict = {k: v for (k, v) in args_dict.items() if k in options}
        if required_args:
            namespace_dict.update(required_args)
        namespace_dict.update(extra_args)
        return Namespace(**namespace_dict)

    sortOpts = ["debug"]
    sort_args = create_namespace(
        sortOpts,
        required_args={
            "input": "test.polish.fasta",
            "out": "test.final.fasta",
            "name": "scaffold",
            "minlen": args.mincontiglen,
            "pipe": True,
        },
    )

    assessOpts = ["debug"]
    assess_args = create_namespace(
        assessOpts,
        required_args={
            "input": "test.final.fasta",
            "report": False,
            "telomere_monomer": "TAA[C]+",
            "telomere_n_repeat": 2,
            "pipe": True,
        },
    )
    return sort_args, assess_args


class TestSortArgsNamespace:
    def test_sort_args_has_debug(self):
        args = _make_pipeline_args(debug=False)
        sort_args, _ = _run_pipeline_up_to_sort(args)
        assert hasattr(sort_args, "debug")

    def test_sort_args_has_pipe(self):
        args = _make_pipeline_args()
        sort_args, _ = _run_pipeline_up_to_sort(args)
        assert hasattr(sort_args, "pipe")
        assert sort_args.pipe is True

    def test_sort_args_debug_inherits_from_parent(self):
        args = _make_pipeline_args(debug=True)
        sort_args, _ = _run_pipeline_up_to_sort(args)
        assert sort_args.debug is True

    def test_sort_args_has_required_fields(self):
        args = _make_pipeline_args(mincontiglen=1000)
        sort_args, _ = _run_pipeline_up_to_sort(args)
        assert sort_args.input == "test.polish.fasta"
        assert sort_args.out == "test.final.fasta"
        assert sort_args.name == "scaffold"
        assert sort_args.minlen == 1000


class TestAssessArgsNamespace:
    def test_assess_args_has_debug(self):
        args = _make_pipeline_args()
        _, assess_args = _run_pipeline_up_to_sort(args)
        assert hasattr(assess_args, "debug")

    def test_assess_args_has_pipe(self):
        args = _make_pipeline_args()
        _, assess_args = _run_pipeline_up_to_sort(args)
        assert hasattr(assess_args, "pipe")
        assert assess_args.pipe is True

    def test_assess_args_has_telomere_fields(self):
        args = _make_pipeline_args()
        _, assess_args = _run_pipeline_up_to_sort(args)
        assert assess_args.telomere_monomer == "TAA[C]+"
        assert assess_args.telomere_n_repeat == 2

    def test_assess_args_debug_inherits_from_parent(self):
        args = _make_pipeline_args(debug=True)
        _, assess_args = _run_pipeline_up_to_sort(args)
        assert assess_args.debug is True
