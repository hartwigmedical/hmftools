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
        // General case: "chr${chr}:g.${gdnaPos}${gdnaRef}>${dnaAlt}/c.${cdnaPos}${cdnaRef}>${cdnaAlt}/p.${aaRef}${aaPos}{aaAlt}"
        //  For MNV the g. part looks like ${gdnaPosStart}_${gdnaPosEnd}del${ref}ins${alt}
        String[] chromosomeAndGDNA = (field.split("/")[0]).split(":");

        // Remove "chr" from the chromosome
        builder.chromosome(chromosomeAndGDNA[0].substring(3));

        // Remove "g." from the gdna annotation
        String gdna = chromosomeAndGDNA[1].substring(2);

        if (gdna.contains("_")) {
            populateForMNV(builder, gdna);
        } else {
            populateForSNV(builder, gdna);
        }
    }

    private static void populateForMNV(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String gdna) {
        String[] gdnaParts = gdna.split("_");
        builder.gdnaPosition(Integer.parseInt(gdnaParts[0]));

        int delStart = gdnaParts[1].indexOf("del");
        int insStart = gdnaParts[1].indexOf("ins");
        builder.gdnaRef(gdnaParts[1].substring(delStart+3, insStart));
        builder.gdnaAlt(gdnaParts[1].substring(insStart+3));
    }

    private static void populateForSNV(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String gdna) {
        StringBuilder gdnaPos = new StringBuilder();
        StringBuilder gdnaRef = new StringBuilder();
        StringBuilder gdnaAlt = new StringBuilder();

        boolean foundNonInteger = false;
        boolean foundRefToAltChar = false;
        for (int i = 0; i < gdna.length(); i++) {
            char charToEvaluate = gdna.charAt(i);

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
}
