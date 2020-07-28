from __future__ import annotations

import itertools
import pickle
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Set, Tuple

import hmmer_reader
from fasta_reader import open_fasta
from hmmer import read_domtbl

from iseq.gff import read as read_gff

from ._confusion import ConfusionMatrix

__all__ = ["ProfMark"]

Sample = NamedTuple("Sample", [("prof_acc", str), ("target_id", str), ("idx", int)])


class ProfMark:
    def __init__(
        self,
        hmmer_file: Path,
        target_file: Path,
        domtblout_file: Path,
        output_file: Path,
    ):
        sample_space: Set[Sample] = generate_sample_space(hmmer_file, target_file)
        true_samples = get_domtblout_samples(domtblout_file)
        ordered_sample_hits = get_ordered_output_samples(output_file)
        sample_space |= true_samples | set(ordered_sample_hits)

        self._sample_space: Set[Sample] = sample_space
        self._true_samples: Set[Sample] = true_samples
        self._ordered_sample_hits: List[Sample] = ordered_sample_hits

    def confusion_matrix(
        self, solut_space="prof-target", solut_space_idx=True
    ) -> ConfusionMatrix:
        from numpy import zeros

        if solut_space == "prof-target":
            sample_space, true_samples, ordered_sample_hits = self._prof_target_space()
        elif solut_space == "prof":
            sample_space, true_samples, ordered_sample_hits = self._prof_space()
        else:
            assert solut_space == "target"
            sample_space, true_samples, ordered_sample_hits = self._target_space()

        if not solut_space_idx:
            sample_space = set(s for s in sample_space if s.idx == 0)
            true_samples = set(s for s in true_samples if s.idx == 0)
            ordered_sample_hits = [s for s in ordered_sample_hits if s.idx == 0]

        sample_space_id = {s: i for i, s in enumerate(sorted(sample_space))}
        true_sample_ids = [sample_space_id[k] for k in true_samples]

        P = len(true_sample_ids)
        N = len(sample_space_id) - P
        sorted_samples = zeros(N + P, int)
        for i, sample in enumerate(ordered_sample_hits):
            sorted_samples[i] = sample_space_id[sample]

        for i, sample in enumerate(sample_space - set(ordered_sample_hits)):
            sorted_samples[i + len(ordered_sample_hits)] = sample_space_id[sample]

        return ConfusionMatrix(true_sample_ids, N, sorted_samples)

    def _prof_target_space(self):
        return self._sample_space, self._true_samples, self._ordered_sample_hits

    def _prof_space(self):
        sample_space = set()
        for k, n in prof_count(self._sample_space).items():
            for i in range(n):
                sample_space.add(Sample(k, "", i))

        true_samples = set()
        for k, n in prof_count(self._true_samples).items():
            for i in range(n):
                true_samples.add(Sample(k, "", i))

        ordered_sample_hits = []
        count: Dict[str, int] = {}
        for sample in self._ordered_sample_hits:
            acc = sample.prof_acc
            count[acc] = count.get(acc, -1) + 1
            ordered_sample_hits.append(Sample(acc, "", count[acc]))

        return sample_space, true_samples, ordered_sample_hits

    def _target_space(self):
        sample_space = set()
        for k, n in target_count(self._sample_space).items():
            for i in range(n):
                sample_space.add(Sample("", k, i))

        true_samples = set()
        for k, n in target_count(self._true_samples).items():
            for i in range(n):
                true_samples.add(Sample("", k, i))

        ordered_sample_hits = []
        count: Dict[str, int] = {}
        for sample in self._ordered_sample_hits:
            tgt = sample.target_id
            count[tgt] = count.get(tgt, -1) + 1
            ordered_sample_hits.append(Sample("", tgt, count[tgt]))

        return sample_space, true_samples, ordered_sample_hits

    def write_pickle(self, filepath: Path):
        with open(filepath, "wb") as file:
            pickle.dump(self, file)

    @staticmethod
    def read_pickle(filepath: Path):
        with open(filepath, "rb") as file:
            return pickle.load(file)


def prof_count(sample_space: Iterable[Sample]) -> Dict[str, int]:
    prof_count: Dict[str, int] = {}
    for sample in sample_space:
        acc = sample.prof_acc
        prof_count[acc] = prof_count.get(acc, 0) + 1
    return prof_count


def target_count(sample_space: Iterable[Sample]) -> Dict[str, int]:
    target_count: Dict[str, int] = {}
    for sample in sample_space:
        acc = sample.target_id
        target_count[acc] = target_count.get(acc, 0) + 1
    return target_count


def generate_sample_space(hmmer_file, target_file) -> Set[Sample]:
    df = hmmer_reader.fetch_metadata(hmmer_file)
    prof_accs = df["ACC"].tolist()
    target_ids = []
    for target in open_fasta(target_file):
        target_ids.append(target.id.split("|")[0])

    samples = [Sample(a, i, 0) for a, i in itertools.product(prof_accs, target_ids)]
    return set(samples)


def get_domtblout_samples(domtblout_file) -> Set[Sample]:
    samples = []
    sample_idx: Dict[Tuple[str, str], int] = {}
    for row in iter(read_domtbl(domtblout_file)):
        profile_acc = row.target.accession
        target_id = row.query.name.split("|")[0]

        if (profile_acc, target_id) not in sample_idx:
            sample_idx[(profile_acc, target_id)] = 0

        idx = sample_idx[(profile_acc, target_id)]
        samples.append(Sample(profile_acc, target_id, idx))
        sample_idx[(profile_acc, target_id)] += 1

    return set(samples)


def get_ordered_output_samples(output_file) -> List[Sample]:
    samples: List[Tuple[Sample, float]] = []
    sample_idx: Dict[Tuple[str, str], int] = {}
    for item in read_gff(output_file).items():
        profile_acc = dict(item.attributes_astuple())["Profile_acc"]
        evalue = float(dict(item.attributes_astuple())["E-value"])
        target_id = item.seqid.split("|")[0]

        if (profile_acc, target_id) not in sample_idx:
            sample_idx[(profile_acc, target_id)] = 0

        idx = sample_idx[(profile_acc, target_id)]
        samples.append((Sample(profile_acc, target_id, idx), evalue))
        sample_idx[(profile_acc, target_id)] += 1

    return [val[0] for val in sorted(samples, key=lambda i: i[1])]
