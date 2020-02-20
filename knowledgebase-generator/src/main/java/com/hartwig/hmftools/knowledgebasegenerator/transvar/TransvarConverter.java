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
            }
            else if (infoField.contains("candidate_codons")) {
                String candidates = infoField.split("=")[1];
                builder.addCandidateCodons(candidates.split(","));
            }
        }
    }
}
