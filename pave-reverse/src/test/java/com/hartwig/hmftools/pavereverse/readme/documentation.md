# ReversePave

**This was written for a code walk-through. It has not been kept up-to-date.**

This is a tool for evaluating the possible DNA changes that give rise to a given protein variant.

## Inputs
Right now, the main class is `ReversePave`, which is created with a path to the ensembl data, and a genome.
The main method is `calculateVariant(String gene, String proteinVariant)` and a couple of 
slightly different options.

We can change things to add a `main` function, batch processing, etc, as required.

## Outputs
The output format was decided after looking at `TransvarConverter` in the `serve` project. 
We can change it to provide more or less information, if required.
Currently, the output is a `BaseSequenceVariants` object, which consists of:
- `TranscriptData` for the transcript found
- `Chromosome`
- a collection of `BaseSequenceChange` objects

In turn, a `BaseSequenceChange` has position, ref, alt and chromosome (which is redundant, perhaps) fields.

Here's a code example:
```java
ReversePave reversePave = new ReversePave(ensemblDataDir, genome);
BaseSequenceVariants variant = reversePave.calculateVariant("MTOR", "L2230V");

assertEquals("ENST00000361445", variant.transcriptId());
assertEquals("1", variant.Chromosome);
checkChanges(variant,
        basesChange("A", "C", "chr1", 11_122_101),
        basesChange("TAA", "GAC", "chr1", 11_122_099),
        basesChange("TAA", "CAC", "chr1", 11_122_099),
        basesChange("TAA", "AAC", "chr1", 11_122_099)
);
```

## Transcript selection
Given a `(gene,variant)` pair, a protein transcript needs to be selected.
This is done as follows: the transcripts in the ensembl data are filtered for those that satisfy the requirements
of the variant. For example, with `(BRAF, V600E)`, we are only interested in `BRAF` transcripts that have a `V` at
position 600.
If:
- no transcript found, an exception is thrown
- the canonical transcript is amongst the filtered transcripts, it is selected
- a single non-canonical transcript is found, it is selected
- multiple non-canonical transcripts are found, but not the canonical transcript, then the hotspots
for each are calculated, and if they match, are returned, otherwise an exception is thrown.

We can easily add an option to specify the transcript to be used, if that is required.

## Only in-exon variants are considered
Only base sequence changes within a single exon are considered. If the first codon in an change region
begins in the previous exon, only changes that conserve the portion of the codon that are
in the previous exon are considered. Similarly, if the last codon in a change regions
ends in the exon following the current exon, then only changes that preserve the portion of the codon 
that is within the next exon are considered.

For example, consider `TET2:p.E1268D`. The amino acids and bases around this change are:
` ...N E E/E R T...  ...AAT GAA GA/G AGA ACT...  exon6/exon7, G at start of exon7 is at 105,259,619`
The codons `CAA` and `CAG` code for `Q`. However, to change `GA/G` to `CA/C` would require changes
in both exons 6 and 7, which is considered very unlikely. Therefore, we choose just the 
`G>C` variant.

## Left-most variants sought
An initial change position is evaluated based on the position of the first amino acid changed in 
the protein variant. We look at all potential nucleotide changes at this position that
could give rise to the required variant. We then step one base to the left, and
consider all possible changes at that position. This process continues for up to 32
bases to the left (depending on the type of variant and its location).

A nucleotide change is `compatible` if it gives rise to the required amino acid sequence.
Two nucleotide changes are `equivalent` if they have the same result when applied to the
exon in which the change occurs. Our algorithm returns the left-most nucleotide change
in each equivalence class of compatible changes.

Consider the variant `D464_Y467del` of the gene`PIK3R1`:
```text
E   Y   D   R   L   Y   E
GAA TAT GAT AGA TTA TAT GAA
Exon starts at 68_293_709
GATCAAGTTGTCAAAGAAGATAATATTGAAGCTGTAGGGAAAAAATTACATGAATATAACACTCAGTTTCAAGAAAAAAGTCGAGAATAT precedes codon for D
So the codon for D is at 68_293_709 + 90 = 68_293_799

Possible changes:
GAA TAT GAT AGA TTA TAT GAA -> GAA TAT GA                A
GAA TAT GAT AGA TTA TAT GAA -> GAA TAT G                AA
GAA TAT GAT AGA TTA TAT GAA -> GAA TAT                 GAA
GAA TAT GAT AGA TTA TAT GAA -> GAA TA                T GAA
GAA TAT GAT AGA TTA TAT GAA -> GAA T                AT GAA
GAA TAT GAT AGA TTA TAT GAA -> GAA                 TAT GAA 
GAA TAT GAT AGA TTA TAT GAA -> GA                A TAT GAA 
The last is the left-most one. It is {ref=ATATGATAGATT, Alt=, chr=5, position=68_293_795}
All variants have the same residual:GAA TAT GAA, so only this one is returned.
```

Here's an example with three different equivalence classes of compatible variants:
```text
ARID1A P20_P21dup
N   P   P   P   P   P   P   S 
AAC CCG CCG CCG CCG CCG CCC TCG
N starts at 26_696_446

The calculated hotspots are:
hotspot("A", "ACCCGCC", "chr1", 26_696_447),
hotspot("C", "CCCGCCG", "chr1", 26_696_448),
hotspot("G", "GCCGCCC", "chr1", 26_696_460)
```

## Testing
### Unit tests
Most of the classes have reasonably comprehensive unit tests.

### Integration-style test
By using a reduced genome and reduced ensembl data we can run quick tests that
use real data. These are in `ReversePaveTest`. Many of these were written using
Transvar output as a starting point and a convenient way of comparing the new
program with Transvar.

### Automated comparison with Transvar output
The Transvar-generated JSON files from the data-vm machine have been used to compare the results of the
new codebase with Transvar. (See below.)

## Comparison with Transvar
The serve.json file on the data vm contains 38,786 variants for the version 38 genome.
The hotspots for these variants were calculated using the new code, and comparisons
were made. For each variant, the hotspots from Transvar were considered the same
as those generated by the new code if there was at least one hotspot in common.
It's very hard to make more specific comparisons because for different kinds of variants
we report different hotspots.

Results:
```text
OVERALL STATS
*************
Number of genes: 1243
Number of annotations not parsed: 0
Number with no transcript found: 87
Number with multiple non-canonical transcripts giving different results: 237
Number processed without error: 38796
Number for which no variant could be calculated: 7
Number with processing error: 0
Number of annotations with similar hotspots: 38034
Difference types: {
    StartLost=1, 
    Frameshift=290, 
    Duplication=76, 
    SingleAminoAcidVariant=152, 
    Deletion=9, 
    Insertion=88, 
    DeletionInsertion=133, 
    StopGained=13
}
```
The differences are discussed on a per-type basis below, but:
- we generally report further-left frameshifts than Transvar, and it's hard to automatically reconcile these, and there are lots of them
- with duplications and insertions, there are generally many options for the inserted bases, and it's hard to 
reconcile these automatically
- many of the other differences are due to Transvar using a different transcript from what the new code selected (often canonical). There's
no information in the json file as to which transcript Transvar used, so it's hard to reconcile these automatically.

11 of the annotations resulted in processing errors. These are bugs that are being worked on.

### Frameshifts 
We only report single base deletions, and find the left-most one that will result
in the required protein change.

```text
Difference for: ARAF T253fs, calculated:
[BaseSequenceChange{Ref='GA', Alt='G', Chromosome='chrX', Position=47567013}]
from transvar:
[BaseSequenceChange{Ref='A', Alt='AG', Chromosome='chrX', Position=47567014},
BaseSequenceChange{Ref='A', Alt='AG', Chromosome='chrX', Position=47567015},
BaseSequenceChange{Ref='AA', Alt='A', Chromosome='chrX', Position=47567014},
BaseSequenceChange{Ref='A', Alt='AA', Chromosome='chrX', Position=47567015},
BaseSequenceChange{Ref='A', Alt='AC', Chromosome='chrX', Position=47567014},
BaseSequenceChange{Ref='A', Alt='AT', Chromosome='chrX', Position=47567014},
BaseSequenceChange{Ref='A', Alt='AT', Chromosome='chrX', Position=47567015},
BaseSequenceChange{Ref='AAC', Alt='A', Chromosome='chrX', Position=47567014},
BaseSequenceChange{Ref='A', Alt='AA', Chromosome='chrX', Position=47567014}]
```

Explanation:
```text
G   T   P...
GGA ACC CCC...the 2nd G of the G is at 47567014
```

### Deletions

```text
Difference for: ARAF Q349_A350del,
calculated:
[BaseSequenceChange{Ref='CCAGGCT', Alt='C', Chromosome='chrX', Position=47567400},
BaseSequenceChange{Ref='GCAGGCC', Alt='G', Chromosome='chrX', Position=47567394}]
from transvar:
[BaseSequenceChange{Ref='CCAGGCT', Alt='C', Chromosome='chrX', Position=47567400}]
```

```text
Difference for: CASP8 T272del,
calculated ENST00000673742:
[BaseSequenceChange{Ref='GACC', Alt='G', Chromosome='chr2', Position=201284823},
BaseSequenceChange{Ref='CACG', Alt='C', Chromosome='chr2', Position=201284826}]
from transvar:
[BaseSequenceChange{Ref='ACTG', Alt='A', Chromosome='chr2', Position=201284870},
BaseSequenceChange{Ref='CTGC', Alt='C', Chromosome='chr2', Position=201284871}]
```
Explanation:  Transvar uses a different transcript


### Duplications
```text
Duplication difference for: HRAS G13dup,
calculated: ENST00000311189 canonical: true
[BaseSequenceChange{Ref='A', Alt='ACAC', Chromosome='chr11', Position=534282}]
from transvar:
[BaseSequenceChange{Ref='C', Alt='CACC', Chromosome='chr11', Position=534283}]
```

Explanation: there are many options for base sequences to produce G and the new code
happened to pick a different one from Transvar.
```text
<- 14 13 12 11... V   G   G   A...
<-                GTG TGG CGG CCG....
->                CAC ACC GCC GGC...A of V is at 534282
reversePave: CACACC ACC GCC = GTG TGG TGG CGG = V G G G A..
```

### Start lost changes
StartLost difference for: PPARG M1?. Transvar used a different transcript. Not sure why.

### Insertions
There can be enormously large numbers of options for base sequences giving rise to a 
string of amino acids. To handle this, the new code just chooses one such option and
finds the left most place where it can be inserted and still produce the overall required
protein change. This almost always produces a different result from that found by Transvar.
Note: an alternative would have been to choose the left-most insertion point for each
alternative and then choose an alternative that had the left most point amongst all found,
but this is not really feasible and would be of questionable value.

```text
Difference for: AKT2 L78_Q79insHANTFVIRCL,
calculated: ENST00000392038 canonical: true
[BaseSequenceChange{Ref='G', Alt='GTAAACACCTGATTACGAAGGTGTTTGCGTG', Chromosome='chr19', Position=40255210}]
from transvar:
[BaseSequenceChange{Ref='G', Alt='GAAGACACCTGATTACAAATGTGTTTGCATG', Chromosome='chr19', Position=40255210}]
>> different choices of codons for the required amino acids
```

Sometimes the new code finds an insertion point further to the left than that returned
by Transvar:
```text
Difference for: ERBB2 E698_P699insLL,
calculated: ENST00000269571 canonical: true
[BaseSequenceChange{Ref='G', Alt='GTTATTA', Chromosome='chr17', Position=39723546}]
from transvar:
[BaseSequenceChange{Ref='C', Alt='CTTCTTC', Chromosome='chr17', Position=39723547}]

V E P... GTG GAG CCG...second G of E is at 39723546
transvar: GTG GAG CTTCTTCCG...V E L L P...
reversePave: GTG GAGTTATTA CCG...V E L L P...
```

### Single Amino Acid Differences
Whenever a difference has been found, it has been because Transvar uses a different
transcript.

### Stop Gained
Whenever a difference has been found, it has been because Transvar uses a different
transcript.

### Deletion insertion
No differences found.

## Outstanding issues
No special treatment for deletion-insertion changes of length two, which could potentially be 
due to small number of base changes.

## Discussion points
Are there any algorithmic changes required?

Do we need to change the API?

Is any kind of stand-alone option required?

Is batch processing required?

Is there a way of assessing how the new code would affect the pipeline? Peter has mentioned a database
containg variants and their hotspots.