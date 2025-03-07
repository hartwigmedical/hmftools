# Reverse Pave
This is a library to calculate the possible coding variants that give rise to a protein variant.

It also has a stand-alone mode for batch operations and testing.

## Inputs
A `ReversePave` instance can be created from an Ensembl cache and genome:
```java
ReversePave reversePave = new ReversePave(mConfig.mEnsemblCache, mConfig.mRefGenome);
```
The instance can then be used to calculate the coding variants that give rise to a protein-level change.
A protein transcript can be specified:
```java
BaseSequenceVariants calculatedVariants = reversePave.calculateVariant(gene, transcript, variant);
```
Alternatively, just the gene and variant can be supplied. 
In this case the transcripts in the ensembl data are filtered for those that satisfy the requirements
of the variant. For example, with `(BRAF, V600E)`, we are only interested in `BRAF` transcripts that have a `V` at
position 600.
If:
- no transcript found, a warning is logged and `null` is returned
- the canonical transcript is amongst the filtered transcripts, it is selected
- a single non-canonical transcript is found, it is selected
- multiple non-canonical transcripts are found, but not the canonical transcript, then the 
set of base level changes for each is calculated, and if these sets match, they are returned, 
otherwise a warning is logged and `null` is returned.

## Outputs
The output is a `BaseSequenceVariants` object, which consists of:
- `TranscriptData` for the transcript found
- `Chromosome`
- a collection of `BaseSequenceChange` objects

In turn, a `BaseSequenceChange` has position, ref, alt and chromosome (which is redundant) fields.

## Stand-alone operation
Reverse Pave has three stand-alone modes of operation:
- process an input TSV file of gene,transcript,hgvs_protein records;
- process a `serve.json` file containing annotations and their possible underlying variants as calculated by the Transvar software;
- do round-trip processing of a Pave TSV results file.
These are explained in more detail below.

### Batch processing of HGVS annotations
The input is a TSV file containing HGVS annotations:

| GENE | TRANSCRIPT      | HGVS_PROTEIN |
|------|-----------------|--------------|
| MTOR | ENST00000361445 | L2230V       |
| BRAF |                 | V600E        |
| BRAF |                 | V600E        |

The output is the input plus a new column containing a tab-separated list of the `BaseSequenceChange`s that were
calculated by `ReversePave`.

### Processing a Serve json file
The `serve.json` file contains protein annotations and the underlying possible changes (called 'hotspots')
as calculated by the Transvar program. This file can be fed into `ReversePave` to produced an output
file similar to that produced when in batch mode.

### Round-trip testing mode
The TSV files produced by Pave contain HGVS annotations for missense mutations.
The combination of the mutation data and the annotations can be used as a way of testing
`ReversePave`: the changes calculated by `ReversePave` should include the original mutation.
In round-trip testing mode, `ReversePave` logs those missense mutations in a Pave output
file for which the calculated changes do not include the actual change. 

Note that there can be legitimate reasons for a match not being found:
- the mutation might be a mixture of a missense and a synonymous mutation in adjacent codons,
- `ReversePave` only returns single-base deletions that give rise to frameshifts,
- various other reasons.

### Arguments

| Argument           |                                  | Description                                    |
|--------------------|----------------------------------|------------------------------------------------|
| ref_genome         | Mandatory                        | Reference genome fasta file                    |
| ensembl_data_dir   | Mandatory                        | Path to Ensembl data cache directory           |
| ref_genome_version | Mandatory                        | 37 (default) or 38                             |
| mode               | batch, round_trip, or serve_json | Operational mode, default is batch             |
| vcf_input          | Optional                         | Input variant VCF, for batch mode              |
| tsv_input          | Optional                         | Input Pave TSV, for round_trip mode            |
| serve_json_input   | Optional                         | Input serve json, for serve_json mode          |
| tsv_output         | Optional                         | Output TSV file, for batch and serve_json mode |

```commandline

```
