package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import org.jetbrains.annotations.NotNull;

final class TransvarConverter {

    private static final String FIELD_DELIMITER = "\t";

    private TransvarConverter() {
    }

    @NotNull
    static TransvarRecord toTransvarRecord(@NotNull String transvarLine) {
        ImmutableTransvarRecord.Builder builder = ImmutableTransvarRecord.builder();

        String[] fields = transvarLine.split(FIELD_DELIMITER);

        populateTranscript(builder, fields[1]);
        populateCoordinates(builder, fields[4]);
        populateCodonInfo(builder, fields[5]);

        return builder.build();
    }

    private static void populateTranscript(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String field) {
        // Field looks like "${transcript} (protein_coding)"
        String[] parts = field.trim().split(" ");
        builder.transcript(parts[0]);
    }

    private static void populateCoordinates(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String field) {
        // Field looks like "chr${chr}:g.${gdnaPos}${gdnaRef}>${dnaAlt}/c.${cdnaPos}${cdnaRef}>${cdnaAlt}/p.${aaRef}${aaPos}{aaAlt}"
        String[] chromosomeAndGNDA = (field.split("///")[0]).split(":");

        builder.chromosome(chromosomeAndGNDA[0].substring(3));
    }

    private static void populateCodonInfo(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String field) {
        // Field is semicolon-separated. Relevant fields look like "reference_codon=${ref};candidate_codons=${alt1},${alt2},...,${altN};"

    }
}
