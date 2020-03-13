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
    private static final String MSG_ERROR_INDICATION_PREFIX = "Error_";

    private static final String RANGE_INDICATOR = "_";
    private static final String DELETION = "del";
    private static final String INSERTION = "ins";
    private static final String DUPLICATION = "dup";

    private TransvarConverter() {
    }

    @Nullable
    static TransvarRecord toTransvarRecord(@NotNull String transvarLine) {
        String[] fields = transvarLine.split(FIELD_DELIMITER);

        String message = fields[MESSAGE_COLUMN];
        if (message.contains(MSG_NO_VALID_TRANSCRIPT_FOUND) || message.trim().startsWith(MSG_ERROR_INDICATION_PREFIX)) {
            return null;
        }

        ImmutableTransvarRecord.Builder builder = ImmutableTransvarRecord.builder();

        populateTranscript(builder, fields[TRANSCRIPT_COLUMN]);
        populateCoordinatesRefAlt(builder, fields[COORDINATES_COLUMN]);
        populateCodonInfo(builder, message);

        TransvarRecord record = builder.build();

        if (isLong(record.gdnaRef()) || isLong(record.gdnaAlt())) {
            // For long indels, transvar gives the length of the indel rather than the exact bases. We ignore such records.
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
            if (gdna.contains(INSERTION) || gdna.contains(DELETION)) {
                populateForInsertionDeletion(builder, gdna);
            } else {
                populateForDuplication(builder, gdna);
            }
        } else {
            populateForSNV(builder, gdna);
        }
    }

    private static void populateForInsertionDeletion(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String gdna) {
        assert gdna.contains(RANGE_INDICATOR);

        String[] gdnaParts = gdna.split(RANGE_INDICATOR);
        long start = Long.parseLong(gdnaParts[0]);
        builder.gdnaPosition(start);

        String delInsPart = gdnaParts[1];
        int delStart = delInsPart.indexOf(DELETION);
        int insStart = delInsPart.indexOf(INSERTION);

        if (delStart < 0 && insStart < 0) {
            throw new IllegalStateException("Cannot process range gDNA as no '" + DELETION + "' or  '" + INSERTION + "' found: " + gdna);
        }

        String gdnaRef = Strings.EMPTY;
        String gdnaAlt = Strings.EMPTY;
        if (insStart >= 0) {
            gdnaAlt = delInsPart.substring(insStart + INSERTION.length());

            if (delStart >= 0) {
                // This should look like '123delCinsG'
                gdnaRef = delInsPart.substring(delStart + DELETION.length(), insStart);
            }
        } else {
            // This should look like '123delC'
            gdnaRef = delInsPart.substring(delStart + DELETION.length());
        }

        // Fill in the indel length in case of an indel
        if (!gdnaRef.isEmpty() && gdnaAlt.isEmpty()) {
            builder.indelLength(gdnaRef.length());
        } else if (gdnaRef.isEmpty() && !gdnaAlt.isEmpty()) {
            builder.indelLength(gdnaAlt.length());
        }

        builder.gdnaRef(gdnaRef);
        builder.gdnaAlt(gdnaAlt);
    }

    private static void populateForDuplication(@NotNull ImmutableTransvarRecord.Builder builder, @NotNull String gdna) {
        assert gdna.contains(RANGE_INDICATOR);

        String[] gdnaParts = gdna.split(RANGE_INDICATOR);
        long start = Long.parseLong(gdnaParts[0]);
        builder.gdnaPosition(start);

        // DUPs simply look like 'start_end' and come with no ref/alt information, but some come with something "dupTTT" appended to it.
        // In both cases we ignore the "dup" part.
        builder.gdnaRef(Strings.EMPTY);
        builder.gdnaAlt(Strings.EMPTY);

        String dupPart = gdnaParts[1];
        if (dupPart.contains(DUPLICATION)) {
            builder.indelLength(1 + (int) (Long.parseLong(dupPart.substring(0, dupPart.indexOf(DUPLICATION))) - start));
        } else if (isLong(dupPart)) {
            builder.indelLength(1 + (int) (Long.parseLong(dupPart) - start));
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
        // SNVs look like 1234T>C
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
