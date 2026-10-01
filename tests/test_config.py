import numpy as np
import pytest

from pleiofdr.cli import exclusion_regions, options_from_config
from pleiofdr.config import ConfigError, TextConfig, parse_cell, parse_matrix

from .conftest import ROOT


def test_parse_matrix_scalars_and_matrices():
    assert parse_matrix("0.05").tolist() == [[0.05]]
    assert parse_matrix("true").tolist() == [[1.0]]
    assert np.isnan(parse_matrix("nan")[0, 0])
    assert parse_matrix("[6 25119106 33854733; 8 7200000 12500000]").tolist() == [
        [6, 25119106, 33854733],
        [8, 7200000, 12500000],
    ]
    assert parse_matrix("[1 0 0; 1 0.5 0 ; 0 0.75 0.75]").shape == (3, 3)
    assert parse_matrix("[1, 2, 3]").tolist() == [[1, 2, 3]]
    assert parse_matrix("[]").shape == (0, 0)
    assert parse_matrix("").shape == (0, 0)


def test_parse_matrix_rejects_code():
    with pytest.raises(ConfigError):
        parse_matrix("[__import__('os')]")
    with pytest.raises(ConfigError):
        parse_matrix("[1 2; 3]")


def test_parse_cell():
    assert parse_cell("{'SSGAC_EDU_2016.mat'}") == ["SSGAC_EDU_2016.mat"]
    assert parse_cell("{'a.mat', 'b.mat' 'c.mat'}") == ["a.mat", "b.mat", "c.mat"]
    assert parse_cell("{'it''s'}") == ["it's"]
    assert parse_cell("{}") == []
    with pytest.raises(ConfigError):
        parse_cell("'a.mat'")
    with pytest.raises(ConfigError):
        parse_cell("{system('ls')}")


def test_load_file_behaviour(tmp_path):
    path = tmp_path / "config.txt"
    path.write_text(
        "# comment\n\n  mafthresh = 0.01  \nstattype=condfdr=ignored\nrandprune=false\n"
    )
    cfg = TextConfig()
    cfg.load_file(path)
    assert cfg.get_num("mafthresh") == 0.01
    assert cfg.get_str("stattype") == "condfdr"  # text after a second '=' is dropped
    assert cfg.get_bool("randprune") is False
    assert cfg.get_bool("perform_gc") is True  # runme.m default

    path.write_text("not_a_key=1\n")
    with pytest.raises(ConfigError, match="no member"):
        TextConfig().load_file(path)


def test_get_num_requires_scalar():
    cfg = TextConfig()
    with pytest.raises(ConfigError):
        cfg.get_num("manh_colorlist")


@pytest.mark.parametrize("name", ["config_default.txt", "config_template.txt"])
def test_shipped_configs_parse(name):
    cfg = TextConfig()
    text = (ROOT / name).read_text()
    if name == "config_template.txt":  # string.Template placeholders used by run_batch.py
        text = text.replace("${exclude_chr_pos}", "6 25119106 33854733")
        text = text.replace("${fdr_thresh}", "0.05").replace("${randprune_n}", "20")
        text = text.replace("${stat_type}", "conjfdr")
    path = ROOT / "tests" / f".{name}"
    path.write_text(text)
    try:
        cfg.load_file(path)
    finally:
        path.unlink()
    opts = options_from_config(cfg)
    assert opts.stattype == "conjfdr"
    assert opts.exclude_ambiguous_snps
    assert cfg.get_mat("exclude_chr_pos").shape == (1, 3)
    np.testing.assert_allclose(opts.manh_colorlist[0], [0.8, 0, 0])


def test_exclusion_regions():
    chrnumvec = np.array([1, 1, 1, 2, 2, 2])
    posvec = np.array([10, 20, 30, 10, 20, 30])
    regions = exclusion_regions(np.array([[2, 15, 30], [1, 0, 20]]), chrnumvec, posvec)
    assert regions == [(4, 5), (0, 1)]
    assert exclusion_regions(np.zeros((0, 0)), chrnumvec, posvec) == []
    with pytest.raises(ValueError, match="matches no SNPs"):
        exclusion_regions(np.array([[3, 0, 1]]), chrnumvec, posvec)
