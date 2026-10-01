from pleiofdr.selftest import main


def test_selftest_passes(tmp_path, capsys):
    assert main(["--keep", str(tmp_path)]) == 0
    assert "self-test passed" in capsys.readouterr().out
    assert (tmp_path / "out_conjfdr" / "result.mat").exists()
