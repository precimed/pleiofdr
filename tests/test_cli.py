import scipy.io

from pleiofdr.cli import main

from .conftest import ROOT, requires_demo


@requires_demo
def test_cli_end_to_end_with_seed(tmp_path):
    config = tmp_path / "config.txt"
    config.write_text(
        "\n".join(
            [
                f"reffile={ROOT}/ref_1kgPhase3eur_LDr2p1_DEMO.mat",
                f"traitfolder={ROOT}",
                "traitfile1=CTG_COG_2018_DEMO.mat",
                "traitname1=COG",
                "traitfiles={'SSGAC_EDU_2016_DEMO.mat'}",
                "traitnames={'EDU'}",
                "stattype=conjfdr",
                "fdrthresh=0.05",
                "randprune_n=5",
                "exclude_chr_pos=[]",
                "onscreen=false",
                f"outputdir={tmp_path / 'results'}",
            ]
        )
        + "\n"
    )
    assert main(["--config", str(config), "--seed", "1"]) == 0
    out = tmp_path / "results"
    expected = [
        "COG_EDU_conjfdr_0.05_loci.csv",
        "COG_EDU_conjfdr_0.05_all.csv",
        "COG_EDU_zscore_conjfdr_0.05_loci.csv",
        "COG_EDU_zscore_conjfdr_0.05_all.csv",
        "COG_EDU_conjfdr_0.05_manhattan.png",
        "COG_EDU_conjfdr_0.05_manhattan.svg",
        "COG_EDU_conjfdr_0.05.log",
        "COG_vs_EDU_qq.png",
        "EDU_vs_COG_qq.svg",
        "COG_vs_EDU_tdr.png",
        "COG_vs_EDU_enrich.png",
        "EDU_vs_COG_lookup.png",
        "result.mat",
    ]
    missing = [name for name in expected if not (out / name).exists()]
    assert not missing
    first = scipy.io.loadmat(out / "result.mat")["fdrmat"]

    # the same seed reproduces the same prune indices and results
    assert main(["--config", str(config), "--seed", "1"]) == 0
    assert (scipy.io.loadmat(out / "result.mat")["fdrmat"] == first).sum() == first.size - (
        first != first
    ).sum()
