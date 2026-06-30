# Ultima and SBX specific routines in SAGE
Due to the nature of Ultima and SBX sequencing data, some routines in SAGE need to be added or significantly overhauled to support these technologies faithfully.

## Ultima
SAGE applies additional processing to Ultima data to account for homopolymer-associated sequencing errors and the differing ways errors are encoded compared to Illumina data.

### Read-context extension

SAGE performs Ultima-specific extension of read-context cores prior to variant quality calculation.

This extension is required because Ultima local realignment (described later) operates by comparing the read and reference sequence across the variant core. This requires the read and reference to be anchored to each other on both sides of the candidate variant.

Thus, SAGE expands read-context cores until homopolymers are completely covered in both the read bases and associated ref bases, with an matching base outside the homopolymer region(s) on each side of the core. Read contexts that cannot be anchored due to clipping or insufficient flanking sequence fail to produce a read context. In append mode, truncation of flanking sequence (to a minimum of 1bp) is permitted to allow previously established contexts to be extended.

### Ultima variant quality modelling

SAGE applies an Ultima-specific variant quality model to estimate the probability that a read supports a true variant rather than a sequencing artefact.

Unlike standard base-quality driven models, the Ultima quality model is defined relative to homopolymers and uses additional per-base BAM tags. Thus, variants are decomposed into their constituent raw Ultima variants (insertions, deletions, complete hompolymer deletions) that collectively synthesise the variant, and the per-read quality is calculated from this decomposition in accordance with the Ultima sequencer design specifications.

Variants that cannot be interpreted as a plausible synthesis of raw Ultima variants (e.g. long novel indels or many MNVs) are instead assigned a fallback fixed quality for each supporting read (28)

### Ultima local realignment

Unlike Illumina data, where alignments generally provide a suitable representation of the underlying sequencing event, Ultima reads frequently contain homopolymer-associated insertions, deletions and alignment ambiguities. As a result, nearby SNVs, MNVs and indels may be represented differently from the event(s) that most likely generated the read.

Before calculating variant quality, SAGE performs Ultima-specific local realignment of the extended variant core. This process evaluates alternative representations of nearby sequence differences, with particular emphasis on homopolymer length changes and complete homopolymer deletions. The resulting representation more closely reflects likely Ultima sequencing error mechanisms than the original aligner representation alone. Each constituent raw Ultima variant is then passed into the Ultima-specific variant quality model as discussed previously.

## SBX
SAGE applies bespoke handling to medium quality (simplex) bases  to treat them in accordance with their empirical error rates as 'second class' citizens.

### Jitter and MSI determination

Jitter determination is assessed solely on the basis of locally duplex evidence.

Reads containing discordant duplex bases or medium-quality bases within the variant core are excluded from jitter matching. Similarly, only these reads are used by Redux to calculate MSI jitter rates.

This allows for more precision in MSI indel calling at low VAFs where low error rates are required for precise calling.

### Tumor quality scoring

For SBX, SAGE classifies qualities 26-29 as medium-quality evidence - this includes all simplex bases (which are marked as qual 27 by Redux).

Medium-quality reads contribute to AD, but do not contribute directly to the Tumor Quality Probability (TQP) model which is the basis of the minTumorQual filter. Instead, they act as a boost to the amount of supporting reads going into the probabilistic TQP calculation after the rest of the model (DP, per-read quals) has been computed, based on the amount of quality-weighted evidence present. For example, if each medium-quality read has half the empirical phred score of a duplex read, the AD in the TQP calculation is boosted by 1 for every 2 medium-quality reads.

This approach allows simplex reads to influence variant support while preventing the TQP calculation from being highly sensitive to per-read qual variance induced by local simplex depth ratios. The resulting  calculation remains driven primarily by high-confidence evidence while retaining sensitivity to variants supported primarily by simplex observations with a duplex anchor.
