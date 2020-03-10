package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class TransvarConverter {

    private static final String FIELD_DELIMITER = "\t";

    private static final int TRANSCRIPT_COLUMN = 1;
    private static final int COORDINATES_COLUMN = 4;
    private static final int MESSAGE_COLUMN = 6;

    private static final String MSG_NO_VALID_TRANSCRIPT_FOUND = "no_valid_transcript_found";
    private static final String MSG_INVALID_MUTATION_STRING = "Error_invalid_mutation_string";

    private static final String RANGE_INDICATOR = "_";
    private static final String DELETION = "del";
    private static final String INSERTION = "ins";

    private TransvarConverter() {
    }

    @Nullable
    static TransvarRecord toTransvarRecord(@NotNull String transvarLine) {
        String[] fields = transvarLine.split(FIELD_DELIMITER);

        if (fields[MESSAGE_COLUMN].contains(MSG_NO_VALID_TRANSCRIPT_FOUND)
                || fields[MESSAGE_COLUMN].contains(MSG_INVALID_MUTATION_STRING)) {
            return null;
        }

        ImmutableTransvarRecord.Builder builder = ImmutableTransvarRecord.builder();

        populateTranscript(builder, fields[TRANSCRIPT_COLUMN]);
        populateCoordinatesRefAlt(builder, fields[COORDINATES_COLUMN]);
        populateCodonInfo(builder, fields[MESSAGE_COLUMN]);

        TransvarRecord record = builder.build();

        if (isLong(record.gdnaRef()) || isLong(record.gdnaAlt())) {
            // For long indels, transvar gives the length of the indel rather than the exact bases.
            return null;
        }

        return record;
    }

    private static void populateTranscript(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String field) {
        // Field looks like "${transcript} (protein_coding)"
        String[] parts = field.trim().split(" ");
        builder.transcript(parts[0]);
    }

    private static void populateCoordinatesRefAlt(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String field) {
        // General case: "chr${chr}:g.${gdnaPos}${gdnaRef}>${dnaAlt}/c.${cdnaPos}${cdnaRef}>${cdnaAlt}/p.${aaRef}${aaPos}{aaAlt}"
        //  For MNV the g. part looks like ${gdnaPosStart}_${gdnaPosEnd}del${ref}ins${alt}
        String[] chromosomeAndGDNA = (field.split("/")[0]).split(":");

        // Remove "chr" from the chromosome
        builder.chromosome(chromosomeAndGDNA[0].substring(3));

        // Remove "g." from the gdna annotation
        String gdna = chromosomeAndGDNA[1].substring(2);

        if (gdna.contains(RANGE_INDICATOR)) {
            if (gdna.contains(DELETION) || gdna.contains(INSERTION)) {
                populateForInsertionDeletion(builder, gdna);
            } else {
                populateForDuplication(builder, gdna);
            }
        } else {
            populateForSNV(builder, gdna);
        }
    }

    private static void populateForInsertionDeletion(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String gdna) {
        String[] gdnaParts = gdna.split(RANGE_INDICATOR);
        builder.gdnaPosition(Long.parseLong(gdnaParts[0]));

        String delInsPart = gdnaParts[1];
        if (delInsPart.contains(DELETION) && delInsPart.contains(INSERTION)) {
            // This should look like 'delCinsG'
            int delStart = delInsPart.indexOf(DELETION);
            int insStart = delInsPart.indexOf(INSERTION);
            builder.gdnaRef(delInsPart.substring(delStart + DELETION.length(), insStart));
            builder.gdnaAlt(delInsPart.substring(insStart + INSERTION.length()));
        } else if (delInsPart.contains(DELETION)) {
            // This should look like 'delC'
            builder.gdnaRef(delInsPart.substring(delInsPart.indexOf(DELETION) + DELETION.length()));
            builder.gdnaAlt(Strings.EMPTY);
        } else if (delInsPart.contains(INSERTION)) {
            // This should look like 'insTTGT'
            builder.gdnaRef(Strings.EMPTY);
            builder.gdnaAlt(delInsPart.substring(delInsPart.indexOf(INSERTION) + INSERTION.length()));
        } else {
            throw new IllegalStateException("Cannot process range gDNA as no '" + DELETION + "' or  '" + INSERTION + "' found: " + gdna);
        }
    }

    private static void populateForDuplication(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String gdna) {
        String[] gdnaParts = gdna.split(RANGE_INDICATOR);
        builder.gdnaPosition(Long.parseLong(gdnaParts[0]));

        if (isLong(gdnaParts[1])) {
            // Assume the variant is a dup with format 'start_end'
            builder.gdnaRef(Strings.EMPTY);
            builder.gdnaAlt(Strings.EMPTY);

            long diff = Long.parseLong(gdnaParts[1]) - Long.parseLong(gdnaParts[0]);
            builder.dupLength(1 + (int) diff);
        } else {
            throw new IllegalStateException("Cannot process duplication for gDNA: " + gdna);
        }
    }

    private static boolean isLong(@NotNull String value) {
        try {
            Long.parseLong(value);
            return true;
        } catch (NumberFormatException exp) {
            return false;
        }
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

        builder.gdnaPosition(Long.parseLong(gdnaPos.toString()));
        builder.gdnaRef(gdnaRef.toString());
        builder.gdnaAlt(gdnaAlt.toString());
    }

    private static void populateCodonInfo(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String field) {
        // Field is semicolon-separated.
        //  SNV relevant fields look like "reference_codon=${ref};candidate_codons=${alt1},${alt2},...,${altN};"
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
