package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class TransvarInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(TransvarInterpreter.class);

    private TransvarInterpreter() {
    }

    @NotNull
    static List<VariantHotspot> extractHotspotsFromTransvarRecord(@NotNull TransvarRecord record, @NotNull HmfTranscriptRegion transcript) {
        if (record.transcript().equals(transcript.transcriptID())) {
            return convertRecordToHotspots(record, transcript.strand());
        } else {
            LOGGER.debug(" Skipped interpretation as transvar transcript '{}' does not match canonical transcript '{}'",
                    record.transcript(),
                    transcript.transcriptID());
            return Lists.newArrayList();
        }
    }

    @NotNull
    @VisibleForTesting
    static List<VariantHotspot> convertRecordToHotspots(@NotNull TransvarRecord record, @NotNull Strand strand) {
        int gdnaCodonIndex = findIndexInRefCodonForGdnaMatch(record, strand);

        List<VariantHotspot> hotspots = Lists.newArrayList();
        for (String candidateCodon : record.candidateCodons()) {
            hotspots.add(fromCandidateCodon(record, candidateCodon, gdnaCodonIndex, strand));
        }

        return hotspots;
    }

    private static int findIndexInRefCodonForGdnaMatch(@NotNull TransvarRecord record, @NotNull Strand strand) {
        String codonCompatibleRef = strand.equals(Strand.FORWARD) ? record.gdnaRef() : flipBase(record.gdnaRef());
        String codonCompatibleAlt = strand.equals(Strand.FORWARD) ? record.gdnaAlt() : flipBase(record.gdnaAlt());

        // We look for the reference codon and candidate codon where the mutation is exclusively the mutation implied by the ref>alt
        for (String candidateCodon : record.candidateCodons()) {
            for (int i = 0; i < 3; i++) {
                if (record.referenceCodon().substring(i, i + 1).equals(codonCompatibleRef) && candidateCodon.substring(i, i + 1)
                        .equals(codonCompatibleAlt)) {
                    boolean match = true;
                    for (int j = 0; j < 3; j++) {
                        if (j != i && !record.referenceCodon().substring(j, j + 1).equals(candidateCodon.substring(j, j + 1))) {
                            match = false;
                        }
                    }
                    if (match) {
                        return i;
                    }
                }
            }
        }

        throw new IllegalStateException("Could not find codon index for GDNA match for " + record);
    }

    @NotNull
    private static VariantHotspot fromCandidateCodon(@NotNull TransvarRecord record, @NotNull String candidateCodon, int gdnaCodonIndex,
            @NotNull Strand strand) {
        String strandAdjustedRefCodon = strand == Strand.FORWARD ? record.referenceCodon() : reverseAndFlip(record.referenceCodon());
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
            stringBuilder.append(flipBase(string.substring(i, i + 1)));
        }
        return stringBuilder.toString();
    }

    @NotNull
    private static String flipBase(@NotNull String base) {
        assert base.length() == 1;

        switch (base) {
            case "A":
                return "T";
            case "T":
                return "A";
            case "G":
                return "C";
            case "C":
                return "G";
        }

        throw new IllegalArgumentException("Cannot flip base: " + base);
    }
}
