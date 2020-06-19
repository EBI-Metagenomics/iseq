import os
from pathlib import Path

import click
import ray
from fasta_reader import FASTAWriter, open_fasta
from nmm import CanonicalAminoAlphabet, DNAAlphabet, GeneticCode, Input, RNAAlphabet

from iseq.hmmsearch import HMMSearch
from iseq.profile import ProfileID
from iseq.protein import ProteinProfile

from .debug_writer import DebugWriter
from .output_writer import OutputWriter
from .pscan import update_gff_file, infer_target_alphabet


@ray.remote
class Worker:
    def __init__(self, alt_filepath, null_filepath, target_abc_name: str, debug: bool):
        self._afile = Input.create(alt_filepath)
        self._nfile = Input.create(null_filepath)
        if target_abc_name == "dna":
            target_abc = DNAAlphabet()
        elif target_abc_name == "rna":
            target_abc = RNAAlphabet()
        else:
            raise RuntimeError()
        self._gcode = GeneticCode(target_abc, CanonicalAminoAlphabet())
        self._output_items = []
        self._codon_seqs = []
        self._amino_seqs = []
        self._debug = debug
        self._debug_table = []
        self._seqids = []

    def set_offset(self, alt_offset, null_offset):
        self._afile.fseek(alt_offset)
        self._nfile.fseek(null_offset)

    def search(self, profids, targets, window):

        for profid in (ProfileID(*i.split("\t")) for i in profids):
            alt = self._afile.read()
            null = self._nfile.read()
            prof = ProteinProfile.create_from_binary(profid, null, alt)
            epsilon = prof.epsilon
            prof.window_length = window

            for tgt in targets:
                seq = prof.create_sequence(tgt.sequence.encode())
                search_results = prof.search(seq)
                ifragments = search_results.ifragments()
                seqid = f"{tgt.defline.split()[0]}"

                for ifrag in ifragments:
                    start = ifrag.interval.start
                    stop = ifrag.interval.stop
                    self._output_items.append(
                        (
                            seqid,
                            prof.profid,
                            start,
                            stop,
                            prof.window_length,
                            {"Epsilon": epsilon},
                        )
                    )
                    codon_frag = ifrag.fragment.decode()
                    self._codon_seqs.append(str(codon_frag.sequence))
                    amino_frag = codon_frag.decode(self._gcode)
                    self._amino_seqs.append(str(amino_frag.sequence))

                if self._debug:
                    self._debug_table += search_results.debug_table()
                self._seqids.append(seqid)

        self._afile.close()
        self._nfile.close()

    def output_items(self):
        return self._output_items

    def codon_seqs(self):
        return self._codon_seqs

    def amino_seqs(self):
        return self._amino_seqs

    def debug_table(self):
        return self._debug_table

    def seqids(self):
        return self._seqids


@click.command()
@click.argument(
    "profile",
    type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True),
)
@click.argument("target", type=click.File("r"))
@click.option(
    "--output",
    type=click.Path(exists=False, dir_okay=False, writable=True, resolve_path=True),
    help="Save results to OUTPUT (GFF format).",
    default="output.gff",
)
@click.option(
    "--ocodon",
    type=click.Path(exists=False, dir_okay=False, writable=True, resolve_path=True),
    help="Save codon sequences to OCODON (FASTA format).",
    default="ocodon.fasta",
)
@click.option(
    "--oamino",
    type=click.Path(exists=False, dir_okay=False, writable=True, resolve_path=True),
    help="Save amino acid sequences to OAMINO (FASTA format).",
    default="oamino.fasta",
)
@click.option(
    "--quiet/--no-quiet", "-q/-nq", help="Disable standard output.", default=False,
)
@click.option(
    "--window",
    type=int,
    help="Window length. Defaults to zero, which means no window.",
    default=0,
)
@click.option(
    "--odebug",
    type=click.File("w"),
    help="Save debug info into a tab-separated values file.",
    default=os.devnull,
)
@click.option(
    "--e-value/--no-e-value",
    help="Enable E-value computation. Defaults to True.",
    default=True,
)
@click.option(
    "--ncpus",
    help="Number of CPUs the user wishes to assign to each worker. Defaults to auto.",
    default="auto",
    type=str,
)
def bscan(
    profile,
    target,
    output,
    ocodon,
    oamino,
    quiet,
    window: int,
    odebug,
    e_value: bool,
    ncpus: str,
):
    """
    Binary scan.
    """

    if ncpus == "auto":
        num_cpus = os.cpu_count()
    else:
        num_cpus = int(ncpus)

    alt_offsets = [int(line.strip()) for line in open(profile + ".alt.idx", "r")]
    null_offsets = [int(line.strip()) for line in open(profile + ".null.idx", "r")]

    num_cpus = min(num_cpus, len(alt_offsets))

    ray.init(
        include_webui=False,
        include_java=False,
        ignore_reinit_error=True,
        num_cpus=num_cpus,
    )

    owriter = OutputWriter(output)
    cwriter = FASTAWriter(ocodon)
    awriter = FASTAWriter(oamino)
    dwriter = DebugWriter(odebug)

    alt_filepath = (profile + ".alt").encode()
    null_filepath = (profile + ".null").encode()
    meta_filepath = (profile + ".meta").encode()

    profids_list = list(
        split([line.strip() for line in open(meta_filepath, "r")], num_cpus)
    )
    alt_offsets = [i[0] for i in split(alt_offsets, num_cpus)]
    null_offsets = [i[0] for i in split(null_offsets, num_cpus)]

    target_abc = infer_target_alphabet(target)
    if isinstance(target_abc, DNAAlphabet):
        tgt_abc_id = "dna"
    elif isinstance(target_abc, RNAAlphabet):
        tgt_abc_id = "rna"
    else:
        raise RuntimeError()

    with open_fasta(target) as fasta:
        targets = list(fasta)

    workers = []
    for aoffset, noffset, profids in zip(alt_offsets, null_offsets, profids_list):
        debug = odebug is not os.devnull
        w = Worker.remote(alt_filepath, null_filepath, tgt_abc_id, debug)
        w.set_offset.remote(aoffset, noffset)
        w.search.remote(profids, targets, window)
        workers.append(w)

    for w in workers:
        output_items = ray.get(w.output_items.remote())
        codon_seqs = ray.get(w.codon_seqs.remote())
        amino_seqs = ray.get(w.amino_seqs.remote())
        for item, cseq, aseq in zip(output_items, codon_seqs, amino_seqs):
            item_id = owriter.write_item(*item)
            cwriter.write_item(item_id, cseq)
            awriter.write_item(item_id, aseq)

    if odebug is not os.devnull:
        for w in workers:
            seqids = ray.get(w.seqids.remote())
            debug_table = ray.get(w.debug_table.remote())
            for seqid, debug_row in zip(seqids, debug_table):
                dwriter.write_row(seqid, debug_row)

    owriter.close()
    cwriter.close()
    awriter.close()
    odebug.close_intelligently()

    if e_value:
        hmmsearch = HMMSearch()
        tbldata = hmmsearch.search(Path(profile), Path(oamino))
        update_gff_file(output, tbldata)

    ray.shutdown()


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m) : (i + 1) * k + min(i + 1, m)] for i in range(n))
