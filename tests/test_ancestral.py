# tests/test_ancestral.py

from pathlib import Path
import shutil
import tarfile

import pysam
import pytest
import vcfpy

from annotate_dafs.ancestral import (
    _get_ancestral_fasta_url,
    _normalize_ensembl_fasta_header,
    download_ancestral_fasta,
    normalize_chromosome_name,
    retrieve_ancestral_allele_from_fasta,
    retrieve_ancestral_allele_from_info,
)


DATA_DIR = Path(__file__).parent / "data"


def test_get_ancestral_fasta_url_grch37():
    """Return the fixed Ensembl release 75 URL for GRCh37."""
    url = _get_ancestral_fasta_url("GRCh37")

    assert url == (
        "https://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/"
        "homo_sapiens_ancestor_GRCh37_e71.tar.bz2"
    )


def test_get_ancestral_fasta_url_grch37_rejects_release():
    """Reject an explicit Ensembl release for GRCh37."""
    with pytest.raises(
        ValueError,
        match="release cannot be specified for GRCh37",
    ):
        _get_ancestral_fasta_url("GRCh37", release=75)


def test_get_ancestral_fasta_url_grch38():
    """Construct the Ensembl ancestral sequence URL for GRCh38."""
    url = _get_ancestral_fasta_url("GRCh38", release=115)

    assert url == (
        "https://ftp.ensembl.org/pub/release-115/fasta/ancestral_alleles/"
        "homo_sapiens_ancestor_GRCh38.tar.gz"
    )


def test_get_ancestral_fasta_url_grch38_requires_release():
    """Require an Ensembl release for GRCh38."""
    with pytest.raises(
        ValueError,
        match="release must be specified for GRCh38",
    ):
        _get_ancestral_fasta_url("GRCh38")


def test_get_ancestral_fasta_url_rejects_unsupported_assembly():
    """Reject unsupported reference genome assemblies."""
    with pytest.raises(
        ValueError,
        match="Unsupported assembly",
    ):
        _get_ancestral_fasta_url("GRCh36")


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        (
            ">ANCESTOR_for_chromosome:GRCh37:1:1:249250621:1",
            "1",
        ),
        (
            ">ANCESTOR_for_chromosome:GRCh37:22:1:51304566:1",
            "22",
        ),
        (
            ">ANCESTOR_for_supercontig:GRCh37:GL000191.1:1:106433:1",
            "GL000191.1",
        ),
        (
            ">ANCESTOR_for_chromosome:GRCh38:X:1:156040895:1",
            "X",
        ),
    ],
)
def test_normalize_ensembl_fasta_header(header, expected):
    """Extract sequence names from Ensembl ancestral FASTA headers."""
    assert _normalize_ensembl_fasta_header(header) == expected


def test_normalize_ensembl_fasta_header_rejects_unexpected_header():
    """Reject an unexpected ancestral FASTA header."""
    with pytest.raises(
        ValueError,
        match="Unexpected Ensembl ancestral FASTA header",
    ):
        _normalize_ensembl_fasta_header(">chr1")


def test_download_ancestral_fasta(tmp_path, monkeypatch):
    """Download, combine, normalize, compress, and index ancestral FASTAs."""
    source_directory = tmp_path / "source"
    source_directory.mkdir()

    chromosome_fasta = source_directory / "chr1.fa"
    chromosome_fasta.write_text(
        ">ANCESTOR_for_chromosome:GRCh37:1:1:10:1\n"
        "ACGTACGTAA\n"
    )

    supercontig_fasta = source_directory / "supercontig.fa"
    supercontig_fasta.write_text(
        ">ANCESTOR_for_supercontig:GRCh37:GL000191.1:1:8:1\n"
        "TTGCAACC\n"
    )

    archive = tmp_path / "ancestral.tar.bz2"

    with tarfile.open(archive, "w:bz2") as tar:
        tar.add(
            chromosome_fasta,
            arcname="ancestral/chr1.fa",
        )
        tar.add(
            supercontig_fasta,
            arcname="ancestral/supercontig.fa",
        )

    def mock_download_file(url, output):
        shutil.copyfile(archive, output)

    monkeypatch.setattr(
        "annotate_dafs.ancestral._download_file",
        mock_download_file,
    )

    output_directory = tmp_path / "output"

    output = download_ancestral_fasta(
        assembly="GRCh37",
        output_directory=output_directory,
    )

    expected_fasta = (
        output_directory / "homo_sapiens_ancestor_GRCh37.fa.gz"
    )

    assert output == expected_fasta
    assert expected_fasta.exists()
    assert Path(f"{expected_fasta}.fai").exists()
    assert Path(f"{expected_fasta}.gzi").exists()

    with pysam.FastaFile(str(expected_fasta)) as fasta:
        assert fasta.references == ["1", "GL000191.1"]
        assert fasta.fetch("1") == "ACGTACGTAA"
        assert fasta.fetch("GL000191.1") == "TTGCAACC"


def make_variant(
    chromosome: str = "chr22",
    position: int = 10,
    reference: str = "C",
    alternate: str = "A",
    info: dict | None = None,
) -> vcfpy.Record:
    """Create a VCF record for testing."""
    return vcfpy.Record(
        CHROM=chromosome,
        POS=position,
        ID=["."],
        REF=reference,
        ALT=[
            vcfpy.Substitution(
                type_="SNV",
                value=alternate,
            )
        ],
        QUAL=30,
        FILTER=["PASS"],
        INFO=info or {},
    )


def test_normalize_chromosome_name():
    """Leading chr prefixes are removed case-insensitively."""
    assert normalize_chromosome_name("chr22") == "22"
    assert normalize_chromosome_name("CHR22") == "22"
    assert normalize_chromosome_name("22") == "22"


def test_retrieve_ancestral_allele_from_info():
    """Ancestral alleles are retrieved and normalized from INFO."""
    variant = make_variant(info={"AA": "c"})

    ancestral_allele = retrieve_ancestral_allele_from_info(
        variant,
        "AA",
    )

    assert ancestral_allele == "C"


def test_retrieve_ancestral_allele_from_info_list():
    """The first value is used when the INFO field contains a list."""
    variant = make_variant(info={"AA": ["a"]})

    ancestral_allele = retrieve_ancestral_allele_from_info(
        variant,
        "AA",
    )

    assert ancestral_allele == "A"


@pytest.mark.parametrize("ancestral_state", ["N", "n", "-", "."])
def test_retrieve_ancestral_allele_from_info_unresolved(ancestral_state):
    """Unresolved ancestral states return None."""
    variant = make_variant(info={"AA": ancestral_state})

    ancestral_allele = retrieve_ancestral_allele_from_info(
        variant,
        "AA",
    )

    assert ancestral_allele is None


def test_retrieve_ancestral_allele_from_info_exclude_low_confidence():
    """Lowercase ancestral alleles can be treated as unresolved."""
    variant = make_variant(info={"AA": "a"})

    ancestral_allele = retrieve_ancestral_allele_from_info(
        variant,
        "AA",
        exclude_low_confidence=True,
    )

    assert ancestral_allele is None


def test_retrieve_ancestral_allele_from_info_missing():
    """A missing ancestral allele raises an informative error."""
    variant = make_variant()

    with pytest.raises(
        ValueError,
        match="Missing ancestral allele field",
    ):
        retrieve_ancestral_allele_from_info(
            variant,
            "AA",
        )


def test_retrieve_ancestral_allele_from_info_empty():
    """An empty ancestral allele raises an informative error."""
    variant = make_variant(info={"AA": ""})

    with pytest.raises(
        ValueError,
        match="Empty ancestral allele field",
    ):
        retrieve_ancestral_allele_from_info(
            variant,
            "AA",
        )


def test_retrieve_ancestral_allele_from_fasta():
    """The ancestral allele is retrieved from the FASTA."""
    fasta = pysam.FastaFile(DATA_DIR / "ancestral.fa")

    try:
        ancestral_allele = retrieve_ancestral_allele_from_fasta(
            fasta,
            "chr22",
            10,
            "C",
        )
    finally:
        fasta.close()

    assert ancestral_allele == "C"


def test_retrieve_ancestral_allele_from_fasta_without_chr_prefix():
    """FASTA retrieval handles chromosome names without a chr prefix."""
    fasta = pysam.FastaFile(DATA_DIR / "ancestral.fa")

    try:
        ancestral_allele = retrieve_ancestral_allele_from_fasta(
            fasta,
            "22",
            10,
            "C",
        )
    finally:
        fasta.close()

    assert ancestral_allele == "C"


def test_retrieve_ancestral_allele_from_fasta_non_snv():
    """FASTA retrieval rejects non-SNV reference alleles."""
    fasta = pysam.FastaFile(DATA_DIR / "ancestral.fa")

    try:
        with pytest.raises(
            ValueError,
            match="only supports SNVs",
        ):
            retrieve_ancestral_allele_from_fasta(
                fasta,
                "chr22",
                10,
                "AT",
            )
    finally:
        fasta.close()

def test_retrieve_ancestral_allele_from_fasta_exclude_low_confidence(
    tmp_path,
):
    """Lowercase FASTA alleles can be treated as unresolved."""
    fasta_path = tmp_path / "ancestral.fa"
    fasta_path.write_text(">chr22\naAAAAAAAA\n")
    pysam.faidx(str(fasta_path))

    fasta = pysam.FastaFile(fasta_path)

    try:
        ancestral_allele = retrieve_ancestral_allele_from_fasta(
            fasta,
            "chr22",
            1,
            "A",
            exclude_low_confidence=True,
        )
    finally:
        fasta.close()

    assert ancestral_allele is None


@pytest.mark.parametrize("ancestral_state", ["N", "n", "-", "."])
def test_retrieve_ancestral_allele_from_fasta_unresolved(
    tmp_path,
    ancestral_state,
):
    """Unresolved FASTA states return None."""
    fasta_path = tmp_path / "ancestral.fa"
    fasta_path.write_text(f">chr22\n{ancestral_state}\n")
    pysam.faidx(str(fasta_path))

    fasta = pysam.FastaFile(fasta_path)

    try:
        ancestral_allele = retrieve_ancestral_allele_from_fasta(
            fasta,
            "chr22",
            1,
            "A",
        )
    finally:
        fasta.close()

    assert ancestral_allele is None