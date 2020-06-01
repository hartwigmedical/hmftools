package com.hartwig.hmftools.serve.transvar;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarDeletion;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarDuplication;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarInsertion;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarRecord;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarSnvMnv;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarAnnotation;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarDeletion;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarInsertion;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarRecord;

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

        TransvarRecord record = createRecord(fields[TRANSCRIPT_COLUMN], fields[COORDINATES_COLUMN], message);

        // (Very) long insertions and deletions report a deleted count rather than a list of bases. We ignore these.
        if (record.annotation() instanceof TransvarDeletion && isLong(((TransvarDeletion) record.annotation()).deletedBases())) {
            return null;
        } else if (record.annotation() instanceof TransvarInsertion && isLong(((TransvarInsertion) record.annotation()).insertedBases())) {
            return null;
        }

        return record;
    }

    @NotNull
    private static TransvarRecord createRecord(@NotNull String transcriptField, @NotNull String coordinateField,
            @NotNull String messageField) {
        ImmutableTransvarRecord.Builder recordBuilder = ImmutableTransvarRecord.builder();

        // Field looks like "${transcript} (protein_coding)"
        recordBuilder.transcript(transcriptField.trim().split(" ")[0]);

        // General case: "chr${chr}:g.${gdnaPos}${gdnaRef}>${dnaAlt}/c.${cdnaPos}${cdnaRef}>${cdnaAlt}/p.${aaRef}${aaPos}{aaAlt}"
        //  For MNV the g. part looks like ${gdnaPosStart}_${gdnaPosEnd}del${ref}ins${alt}
        //  For any type of indel there will be an _ in the range of the gdna coordinates.
        String[] chromosomeAndGDNA = (coordinateField.split("/")[0]).split(":");

        // Remove "chr" from the chromosome
        recordBuilder.chromosome(chromosomeAndGDNA[0].substring(3));

        // Remove "g." from the gdna annotation
        String gdna = chromosomeAndGDNA[1].substring(2);

        if (gdna.contains(RANGE_INDICATOR)) {
            recordBuilder.gdnaPosition(Long.parseLong(gdna.split(RANGE_INDICATOR)[0]));
            TransvarAnnotation annotation;
            if (gdna.contains(INSERTION) || gdna.contains(DELETION)) {
                annotation = annotationForInsertionDeletion(gdna, messageField);
            } else {
                annotation = annotationForDuplication(gdna);
            }

            return recordBuilder.annotation(annotation).build();
        } else {
            return createForSNV(recordBuilder, gdna, messageField);
        }
    }

    @NotNull
    private static TransvarAnnotation annotationForInsertionDeletion(@NotNull String gdna, @NotNull String messageField) {
        assert gdna.contains(RANGE_INDICATOR);

        String delInsPart = gdna.split(RANGE_INDICATOR)[1];
        int delStart = delInsPart.indexOf(DELETION);
        int insStart = delInsPart.indexOf(INSERTION);

        if (delStart < 0 && insStart < 0) {
            throw new IllegalStateException("Cannot process range gDNA as no '" + DELETION + "' or  '" + INSERTION + "' found: " + gdna);
        }

        if (insStart >= 0) {
            String insertedBases = delInsPart.substring(insStart + INSERTION.length());

            if (delStart >= 0) {
                if (delStart + DELETION.length() == insStart) {
                    // This should look like 123delinsGGT
                    // TODO support complex deletions+insertions
                    return new TransvarAnnotation() {};
                } else {
                    // This should look like '123delCinsG'
                    String deletedBases = delInsPart.substring(delStart + DELETION.length(), insStart);
                    assert deletedBases.length() == insertedBases.length();

                    return ImmutableTransvarSnvMnv.builder()
                            .gdnaRef(deletedBases)
                            .gdnaAlt(insertedBases)
                            .referenceCodon(extractReferenceCodonFromMessageField(messageField))
                            .candidateCodons(extractCandidateCodonsFromMessageField(messageField))
                            .build();
                }
            } else {
                return ImmutableTransvarInsertion.builder().insertedBases(insertedBases).build();
            }
        } else {
            // This should look like '123delC'
            String deletedBases = delInsPart.substring(delStart + DELETION.length());
            return ImmutableTransvarDeletion.builder().deletedBases(deletedBases).build();
        }
    }

    @NotNull
    private static TransvarAnnotation annotationForDuplication(@NotNull String gdna) {
        assert gdna.contains(RANGE_INDICATOR);

        String[] gdnaParts = gdna.split(RANGE_INDICATOR);
        long position = Long.parseLong(gdnaParts[0]);

        String dupPart = gdnaParts[1];
        int duplicatedBaseCount;
        if (dupPart.contains(DUPLICATION)) {
            duplicatedBaseCount = 1 + (int) (Long.parseLong(dupPart.substring(0, dupPart.indexOf(DUPLICATION))) - position);
        } else if (isLong(dupPart)) {
            duplicatedBaseCount = 1 + (int) (Long.parseLong(dupPart) - position);
        } else {
            throw new IllegalStateException("Cannot process duplication for gDNA: " + gdna);
        }
        return ImmutableTransvarDuplication.builder().duplicatedBaseCount(duplicatedBaseCount).build();
    }

    @NotNull
    private static TransvarRecord createForSNV(@NotNull ImmutableTransvarRecord.Builder recordBuilder, @NotNull String gdna,
            @NotNull String messageField) {
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

        TransvarAnnotation snvAnnotation = ImmutableTransvarSnvMnv.builder()
                .gdnaRef(gdnaRef.toString())
                .gdnaAlt(gdnaAlt.toString())
                .referenceCodon(extractReferenceCodonFromMessageField(messageField))
                .candidateCodons(extractCandidateCodonsFromMessageField(messageField))
                .build();

        return recordBuilder.gdnaPosition(Long.parseLong(gdnaPos.toString())).annotation(snvAnnotation).build();
    }

    @NotNull
    private static String extractReferenceCodonFromMessageField(@NotNull String messageField) {
        // Field is semicolon-separated.
        //  Relevant fields look like "reference_codon=${ref};candidate_codons=${alt1},${alt2},...,${altN};"
        String[] infoFields = messageField.split(";");

        for (String infoField : infoFields) {
            if (infoField.contains("reference_codon")) {
                return infoField.split("=")[1];
            }
        }

        throw new IllegalStateException("No reference_codon found in message field: " + messageField);
    }

    @NotNull
    private static List<String> extractCandidateCodonsFromMessageField(@NotNull String messageField) {
        // Field is semicolon-separated.
        //  Relevant fields look like "reference_codon=${ref};candidate_codons=${alt1},${alt2},...,${altN};"
        String[] infoFields = messageField.split(";");

        for (String infoField : infoFields) {
            if (infoField.contains("candidate_codons")) {
                String candidates = infoField.split("=")[1];
                return Arrays.asList(candidates.split(","));
            }
        }

        throw new IllegalStateException("No candidate_codons found in message field: " + messageField);
    }

    private static boolean isLong(@NotNull String value) {
        try {
            Long.parseLong(value);
            return true;
        } catch (NumberFormatException exp) {
            return false;
        }
    }
}
