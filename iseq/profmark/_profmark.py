from __future__ import annotations

from .confusion import ConfusionMatrix
import itertools
import hmmer_reader
from iseq.domtblout import DomTBLData
from iseq.gff import read as read_gff
from fasta_reader import open_fasta
from typing import NamedTuple, Set


__all__ = ["Sample", "ProfMark"]

Sample = NamedTuple("Sample", [("prof_acc", str), ("target_id", str)])


class ProfMark:
    def __init__(self, hmmer_file, target_file, domtblout_file, output_file):
        from numpy import zeros

        sample_space: Set[Sample] = generate_sample_space(hmmer_file, target_file)
        true_samples = get_domtblout_samples(domtblout_file)
        samples = get_output_samples(output_file)

        sample_space_id = {s: i for i, s in enumerate(sorted(sample_space))}
        true_samples = [sample_space_id[k] for k in true_samples]

        P = len(true_samples)
        N = len(sample_space) - P
        sorted_samples = zeros(N + P, int)
        for i, sample in enumerate(samples):
            sorted_samples[i] = sample_space_id[sample]

        for i, sample in enumerate(sample_space - samples):
            sorted_samples[i + len(samples)] = sample_space_id[sample]

        self._num_output_samples = len(sorted_samples)
        self._confusion_matrix = ConfusionMatrix(true_samples, P, N, sorted_samples)

    def num_output_samples(self) -> int:
        return self._num_output_samples

    def confusion_matrix(self) -> ConfusionMatrix:
        return self._confusion_matrix


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
