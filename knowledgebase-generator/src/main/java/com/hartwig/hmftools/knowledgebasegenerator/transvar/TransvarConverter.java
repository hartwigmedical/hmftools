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

final class TransvarConverter {

    private static final Logger LOGGER = LogManager.getLogger(TransvarConverter.class);

    private static final String FIELD_DELIMITER = "\t";

    private TransvarConverter() {
    }

    @NotNull
    static List<VariantHotspot> transvarToHotpots(@NotNull String transvarLine, @NotNull HmfTranscriptRegion transcript) {
        TransvarRecord record = TransvarConverter.toTransvarRecord(transvarLine);

        if (record.transcript().equals(transcript.transcriptID())) {
            return convertRecordToHotspots(record, transcript.strand());
        } else {
            LOGGER.warn("Skipped conversion for record as transcript '{}' does not match canonical transcript '{}'",
                    record.transcript(),
                    transcript.transcriptID());
            return Lists.newArrayList();
        }
    }

    @NotNull
    @VisibleForTesting
    static TransvarRecord toTransvarRecord(@NotNull String transvarLine) {
        ImmutableTransvarRecord.Builder builder = ImmutableTransvarRecord.builder();

        String[] fields = transvarLine.split(FIELD_DELIMITER);

        populateTranscript(builder, fields[1]);
        populateCoordinates(builder, fields[4]);
        populateCodonInfo(builder, fields[6]);

        return builder.build();
    }

    private static void populateTranscript(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String field) {
        // Field looks like "${transcript} (protein_coding)"
        String[] parts = field.trim().split(" ");
        builder.transcript(parts[0]);
    }

    private static void populateCoordinates(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String field) {
        // Field looks like "chr${chr}:g.${gdnaPos}${gdnaRef}>${dnaAlt}/c.${cdnaPos}${cdnaRef}>${cdnaAlt}/p.${aaRef}${aaPos}{aaAlt}"
        String[] chromosomeAndGNDA = (field.split("/")[0]).split(":");

        // Remove "chr" from the chromosome
        builder.chromosome(chromosomeAndGNDA[0].substring(3));

        StringBuilder gdnaPos = new StringBuilder();
        StringBuilder gdnaRef = new StringBuilder();
        StringBuilder gdnaAlt = new StringBuilder();

        boolean foundNonInteger = false;
        boolean foundRefToAltChar = false;
        for (int i = 2; i < chromosomeAndGNDA[1].length(); i++) {
            char charToEvaluate = chromosomeAndGNDA[1].charAt(i);

            if (!foundNonInteger) {
                if (Character.isDigit(charToEvaluate)) {
                    gdnaPos.append(charToEvaluate);
                } else {
                    foundNonInteger = true;
                }
            }

            if (foundNonInteger) {
                if (foundRefToAltChar) {
                    gdnaAlt.append(charToEvaluate);
                } else if (String.valueOf(charToEvaluate).equals(">")) {
                    foundRefToAltChar = true;
                } else {
                    gdnaRef.append(charToEvaluate);
                }
            }
        }

        builder.gdnaPosition(Integer.parseInt(gdnaPos.toString()));
        builder.gdnaRef(gdnaRef.toString());
        builder.gdnaAlt(gdnaAlt.toString());
    }

    private static void populateCodonInfo(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String field) {
        // Field is semicolon-separated. Relevant fields look like "reference_codon=${ref};candidate_codons=${alt1},${alt2},...,${altN};"
        String[] infoFields = field.split(";");

        for (String infoField : infoFields) {
            if (infoField.contains("reference_codon")) {
                builder.referenceCodon(infoField.split("=")[1]);
            } else if (infoField.contains("candidate_codons")) {
                String candidates = infoField.split("=")[1];
                builder.addCandidateCodons(candidates.split(","));
            }
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
