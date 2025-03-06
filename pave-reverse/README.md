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

