from __future__ import annotations

import itertools
import pickle
from pathlib import Path
from typing import NamedTuple, Set

import hmmer_reader
from fasta_reader import open_fasta

from iseq.domtblout import DomTBLData
from iseq.gff import read as read_gff

from ._confusion import ConfusionMatrix

__all__ = ["ProfMark"]

Sample = NamedTuple("Sample", [("prof_acc", str), ("target_id", str)])


class ProfMark:
    def __init__(
        self,
        hmmer_file: Path,
        target_file: Path,
        domtblout_file: Path,
        output_file: Path,
    ):
        from numpy import zeros

        sample_space: Set[Sample] = generate_sample_space(hmmer_file, target_file)
        true_samples = get_domtblout_samples(domtblout_file)
        sample_hits = get_output_samples(output_file)

        sample_space_id = {s: i for i, s in enumerate(sorted(sample_space))}
        true_sample_ids = [sample_space_id[k] for k in true_samples]

        P = len(true_sample_ids)
        N = len(sample_space_id) - P
        sorted_samples = zeros(N + P, int)
        for i, sample in enumerate(sample_hits):
            sorted_samples[i] = sample_space_id[sample]

        for i, sample in enumerate(sample_space - sample_hits):
            sorted_samples[i + len(sample_hits)] = sample_space_id[sample]

        self._nhits = len(sample_hits)
        self._confusion_matrix = ConfusionMatrix(true_sample_ids, N, sorted_samples)

    @property
    def nhits(self) -> int:
        return self._nhits

    @property
    def confusion_matrix(self) -> ConfusionMatrix:
        return self._confusion_matrix

    def write_pickle(self, filepath: Path):
        with open(filepath, "wb") as file:
            pickle.dump(self, file)

    @staticmethod
    def read_pickle(filepath: Path):
        with open(filepath, "rb") as file:
            return pickle.load(file)


def generate_sample_space(hmmer_file, target_file) -> Set[Sample]:
    df = hmmer_reader.fetch_metadata(hmmer_file)
    prof_accs = df["ACC"].tolist()
    target_ids = []
    for target in open_fasta(target_file):
        target_ids.append(target.id.split("|")[0])

    samples = [Sample(a, i) for a, i in itertools.product(prof_accs, target_ids)]
    return set(samples)


def get_domtblout_samples(domtblout_file) -> Set[Sample]:
    samples = []
    for row in iter(DomTBLData(domtblout_file)):
        profile_acc = row.target.accession
        target_id = row.query.name.split("|")[0]
        samples.append(Sample(profile_acc, target_id))
    return set(samples)


def get_output_samples(output_file) -> Set[Sample]:
    samples = []
    for item in read_gff(output_file).items():
        profile_acc = dict(item.attributes_astuple())["Profile_acc"]
        target_id = item.seqid.split("|")[0]
        samples.append(Sample(profile_acc, target_id))
    return set(samples)
