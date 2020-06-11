package com.hartwig.hmftools.serve.transvar;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarAnnotation;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarComplexInsertDelete;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarDeletion;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarDuplication;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarInsertion;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarRecord;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarSnvMnv;
import com.hartwig.hmftools.serve.util.AminoAcidFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

class TransvarInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(TransvarInterpreter.class);

    @NotNull
    private final IndexedFastaSequenceFile refGenome;

    @NotNull
    static TransvarInterpreter fromRefGenomeFastaFile(@NotNull String refGenomeFastaFile) throws FileNotFoundException {
        return new TransvarInterpreter(new IndexedFastaSequenceFile(new File(refGenomeFastaFile)));
    }

    private TransvarInterpreter(@NotNull IndexedFastaSequenceFile refGenome) {
        this.refGenome = refGenome;
    }

    @NotNull
    List<VariantHotspot> convertRecordToHotspots(@NotNull TransvarRecord record, @NotNull Strand strand) {
        TransvarAnnotation annotation = record.annotation();
        if (annotation instanceof TransvarSnvMnv) {
            return convertSnvMnvToHotspots(record, (TransvarSnvMnv) annotation, strand);
        } else if (annotation instanceof TransvarDuplication) {
            return convertDuplicationToHotspots(record, (TransvarDuplication) annotation);
        } else if (annotation instanceof TransvarInsertion) {
            return convertInsertionToHotspots(record, (TransvarInsertion) annotation, strand);
        } else if (annotation instanceof TransvarDeletion) {
            return convertDeletionToHotspots(record, (TransvarDeletion) annotation);
        } else if (annotation instanceof TransvarComplexInsertDelete) {
            return convertComplexInsertDeleteToHotspots(record, (TransvarComplexInsertDelete) annotation);
        } else {
            LOGGER.warn("Unrecognized annotation type in transvar record: '{}'. Skipping interpretation.",
                    annotation.getClass().toString());
            return Lists.newArrayList();
        }
    }

    @NotNull
    private List<VariantHotspot> convertSnvMnvToHotspots(@NotNull TransvarRecord record, @NotNull TransvarSnvMnv snvMnv,
            @NotNull Strand strand) {
        List<VariantHotspot> hotspots = Lists.newArrayList();

        if (record.variantSpanMultipleExons()) {
            // In this case we only generate hotspots for a simple SNV.
            if (snvMnv.gdnaRef().length() == 1) {
                hotspots.add(ImmutableVariantHotspotImpl.builder()
                        .chromosome(record.chromosome())
                        .position(record.gdnaPosition())
                        .ref(snvMnv.gdnaRef())
                        .alt(snvMnv.gdnaAlt())
                        .build());
            }
        } else {
            // We need to look up which index of the ref codon is changed (0, 1 or 2) in case of SNV/MNV.
            int gdnaCodonIndex = findIndexInRefCodonForGdnaMatch(snvMnv, strand);

            for (String candidateCodon : snvMnv.candidateCodons()) {
                hotspots.add(fromCandidateCodon(record, snvMnv.referenceCodon(), candidateCodon, gdnaCodonIndex, strand));
            }
        }

        return hotspots;
    }

    @NotNull
    private List<VariantHotspot> convertDuplicationToHotspots(@NotNull TransvarRecord record, @NotNull TransvarDuplication dup) {
        List<VariantHotspot> hotspots = Lists.newArrayList();
        if (!record.variantSpanMultipleExons()) {
            // Dups don't have ref and alt information so need to look it up in ref genome.
            long position = record.gdnaPosition() - 1;
            String preMutatedSequence = refGenome.getSubsequenceAt(record.chromosome(), position, position).getBaseString();

            String dupBases =
                    refGenome.getSubsequenceAt(record.chromosome(), position + 1, position + dup.duplicatedBaseCount()).getBaseString();

            hotspots.add(ImmutableVariantHotspotImpl.builder()
                    .chromosome(record.chromosome())
                    .position(position)
                    .ref(preMutatedSequence)
                    .alt(preMutatedSequence + dupBases)
                    .build());
        } else {
            LOGGER.debug("Duplication spanning multiple exons. Ignoring {}", record);
        }

        return hotspots;
    }

    @NotNull
    private List<VariantHotspot> convertInsertionToHotspots(@NotNull TransvarRecord record, @NotNull TransvarInsertion insertion,
            @NotNull Strand strand) {
        List<VariantHotspot> hotspots = Lists.newArrayList();
        if (!record.variantSpanMultipleExons()) {
            long position = record.gdnaPosition();
            String preMutatedSequence = refGenome.getSubsequenceAt(record.chromosome(), position, position).getBaseString();

            ImmutableVariantHotspotImpl.Builder hotspotBuilder =
                    ImmutableVariantHotspotImpl.builder().chromosome(record.chromosome()).position(position).ref(preMutatedSequence);

            // We assume inserts of length 3 are always (inframe) amino acid inserts, we expand those hotspots.
            if (insertion.insertedBases().length() == 3) {
                for (String trinucleotide : AminoAcidFunctions.allTrinucleotidesForSameAminoAcid(insertion.insertedBases(), strand)) {
                    hotspots.add(hotspotBuilder.alt(preMutatedSequence + trinucleotide).build());
                }
            } else {
                hotspots.add(hotspotBuilder.alt(preMutatedSequence + insertion.insertedBases()).build());
            }
        } else {
            LOGGER.debug("Insertion spanning multiple exons. Ignoring {}", record);
        }
        return hotspots;
    }

    @NotNull
    private List<VariantHotspot> convertDeletionToHotspots(@NotNull TransvarRecord record, @NotNull TransvarDeletion deletion) {
        List<VariantHotspot> hotspots = Lists.newArrayList();
        if (!record.variantSpanMultipleExons()) {
            if (deletion.unalignedGDNAPosition() <= record.gdnaPosition()) {
                ImmutableVariantHotspotImpl.Builder hotspotBuilder = ImmutableVariantHotspotImpl.builder().chromosome(record.chromosome());
                for (long start = deletion.unalignedGDNAPosition(); start <= record.gdnaPosition(); start++) {
                    long adjustedPosition = start - 1;
                    String preMutatedSequence =
                            refGenome.getSubsequenceAt(record.chromosome(), adjustedPosition, adjustedPosition).getBaseString();
                    String deletedSequence = refGenome.getSubsequenceAt(record.chromosome(),
                            adjustedPosition + 1,
                            adjustedPosition + deletion.deletedBaseCount()).getBaseString();

                    hotspots.add(hotspotBuilder.position(adjustedPosition)
                            .ref(preMutatedSequence + deletedSequence)
                            .alt(preMutatedSequence)
                            .build());
                }
            } else {
                LOGGER.warn("Unaligned GDNA > position. Unsure why! Record={}", record);
            }
        } else {
            LOGGER.debug("Deletion spanning multiple exons. Ignoring {}", record);
        }

        return hotspots;
    }

    @NotNull
    private List<VariantHotspot> convertComplexInsertDeleteToHotspots(@NotNull TransvarRecord record,
            @NotNull TransvarComplexInsertDelete insDel) {
        List<VariantHotspot> hotspots = Lists.newArrayList();
        if (!record.variantSpanMultipleExons()) {
            long position = record.gdnaPosition() - 1;
            String preMutatedSequence = refGenome.getSubsequenceAt(record.chromosome(), position, position).getBaseString();

            String deletedBases =
                    refGenome.getSubsequenceAt(record.chromosome(), position + 1, position + insDel.deletedBaseCount()).getBaseString();
            ImmutableVariantHotspotImpl.Builder hotspotBuilder = ImmutableVariantHotspotImpl.builder()
                    .chromosome(record.chromosome())
                    .position(position)
                    .ref(preMutatedSequence + deletedBases);

            hotspots.add(hotspotBuilder.alt(preMutatedSequence + insDel.insertedSequence()).build());
            // For now we only add alternative sequences for insertions of one AA
            if (insDel.insertedSequence().length() == 3) {
                for (String candidateAlternativeSequence : insDel.candidateAlternativeSequences()) {
                    if (!candidateAlternativeSequence.equals(insDel.insertedSequence())) {
                        hotspots.add(hotspotBuilder.alt(preMutatedSequence + candidateAlternativeSequence).build());
                    }
                }
            }
        } else {
            LOGGER.debug("Complex insert/delete spanning multiple exons. Ignoring {}", record);
        }
        return hotspots;
    }

    private static int findIndexInRefCodonForGdnaMatch(@NotNull TransvarSnvMnv snvMnv, @NotNull Strand strand) {
        String codonCompatibleRef = strand.equals(Strand.FORWARD) ? snvMnv.gdnaRef() : reverseAndFlip(snvMnv.gdnaRef());
        String codonCompatibleAlt = strand.equals(Strand.FORWARD) ? snvMnv.gdnaAlt() : reverseAndFlip(snvMnv.gdnaAlt());

        // Function only supports SNV and MNV
        assert codonCompatibleRef.length() == codonCompatibleAlt.length();

        // We look for the candidate codon where the mutation is exclusively the mutation implied by the ref>alt
        int mutLength = codonCompatibleRef.length();
        for (String candidateCodon : snvMnv.candidateCodons()) {
            for (int i = 0; i < 4 - mutLength; i++) {
                if (snvMnv.referenceCodon().substring(i, i + mutLength).equals(codonCompatibleRef) && candidateCodon.substring(i,
                        i + mutLength).equals(codonCompatibleAlt)) {
                    boolean match = true;
                    for (int j = 0; j < 3; j++) {
                        // No match if we find another mismatch outside of the ref->alt range.
                        if ((j - i >= mutLength || j - i < 0) && !snvMnv.referenceCodon()
                                .substring(j, j + 1)
                                .equals(candidateCodon.substring(j, j + 1))) {
                            match = false;
                        }
                    }
                    if (match) {
                        // For reverse strand searches we need to correct for the flipping of MNVs.
                        return strand == Strand.FORWARD ? i : i + mutLength - 1;
                    }
                }
            }
        }

        throw new IllegalStateException("Could not find codon index for gDNA match for " + snvMnv);
    }

    @NotNull
    private static VariantHotspot fromCandidateCodon(@NotNull TransvarRecord record, @NotNull String referenceCodon,
            @NotNull String candidateCodon, int gdnaCodonIndex, @NotNull Strand strand) {
        String strandAdjustedRefCodon = strand == Strand.FORWARD ? referenceCodon : reverseAndFlip(referenceCodon);
        String strandAdjustedCandidateCodon = strand == Strand.FORWARD ? candidateCodon : reverseAndFlip(candidateCodon);
        int strandAdjustedGdnaCodingIndex = strand == Strand.FORWARD ? gdnaCodonIndex : 2 - gdnaCodonIndex;

        int firstMutatedPosition = -1;
        int lastMutatedPosition = -1;
        for (int i = 0; i < 3; i++) {
            if (!strandAdjustedRefCodon.substring(i, i + 1).equals(strandAdjustedCandidateCodon.substring(i, i + 1))) {
                if (firstMutatedPosition == -1) {
                    firstMutatedPosition = i;
                }
                lastMutatedPosition = i;
            }
        }

        String ref = strandAdjustedRefCodon.substring(firstMutatedPosition, lastMutatedPosition + 1);
        String alt = strandAdjustedCandidateCodon.substring(firstMutatedPosition, lastMutatedPosition + 1);

        return ImmutableVariantHotspotImpl.builder()
                .chromosome(record.chromosome())
                .position(record.gdnaPosition() - strandAdjustedGdnaCodingIndex + firstMutatedPosition)
                .ref(ref)
                .alt(alt)
                .build();
    }

    @NotNull
    private static String reverseAndFlip(@NotNull String string) {
        StringBuilder stringBuilder = new StringBuilder();
        for (int i = string.length() - 1; i >= 0; i--) {
            stringBuilder.append(flipBase(string.charAt(i)));
        }
        return stringBuilder.toString();
    }

    private static char flipBase(char base) {
        switch (base) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'G':
                return 'C';
            case 'C':
                return 'G';
        }

        throw new IllegalArgumentException("Cannot flip invalid base: " + base);
    }
}
