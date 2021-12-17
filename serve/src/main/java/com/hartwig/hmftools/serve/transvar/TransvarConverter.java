package com.hartwig.hmftools.serve.transvar;

import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_DELETION;
import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_DUPLICATION;
import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_INSERTION;
import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_RANGE_INDICATOR;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarComplexInsertDelete;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarDeletion;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarDuplication;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarFrameshift;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarInsertion;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarRecord;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarSnvMnv;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarAnnotation;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarDuplication;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarInsertion;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarRecord;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class TransvarConverter {

    private static final String FIELD_DELIMITER = "\t";

    private static final int TRANSCRIPT_COLUMN = 1;
    private static final int COORDINATES_COLUMN = 4;
    private static final int LOCATION_COLUMN = 5;
    private static final int MESSAGE_COLUMN = 6;

    private static final String CSQN_FRAMESHIFT = "Frameshift";

    private static final String MSG_NO_VALID_TRANSCRIPT_FOUND = "no_valid_transcript_found";
    private static final String MSG_ERROR_INDICATION_PREFIX = "Error_";

    private TransvarConverter() {
    }

    @Nullable
    static TransvarRecord toTransvarRecord(@NotNull String transvarLine) {
        String[] fields = transvarLine.split(FIELD_DELIMITER);

        String messageField = fields[MESSAGE_COLUMN];
        if (messageField.contains(MSG_NO_VALID_TRANSCRIPT_FOUND) || messageField.trim().startsWith(MSG_ERROR_INDICATION_PREFIX)) {
            return null;
        }

        TransvarRecord record = createRecord(fields[TRANSCRIPT_COLUMN], fields[COORDINATES_COLUMN], fields[LOCATION_COLUMN], messageField);

        // (Very) long insertions report an inserted count rather than a list of bases. We ignore these.
        if (record.annotation() instanceof TransvarInsertion && isInteger(((TransvarInsertion) record.annotation()).insertedBases())) {
            return null;
        }

        // Duplications are somewhat educated guesses, so need to be sure they refer to a dup at least once in the raw output.
        if (record.annotation() instanceof TransvarDuplication && !transvarLine.contains(HGVS_DUPLICATION)) {
            return null;
        }

        return record;
    }

    @NotNull
    private static TransvarRecord createRecord(@NotNull String transcriptField, @NotNull String coordinateField,
            @NotNull String locationField, @NotNull String messageField) {
        ImmutableTransvarRecord.Builder recordBuilder = ImmutableTransvarRecord.builder();

        // Field looks like "${transcript} (protein_coding)"
        recordBuilder.transcript(transcriptField.trim().split(" ")[0]);
        recordBuilder.variantSpanMultipleExons(variantSpanMultipleExons(locationField));

        // General case: "chr${chr}:g.${gdnaPos}${gdnaRef}>${dnaAlt}/c.${cdnaPos}${cdnaRef}>${cdnaAlt}/p.${aaRef}${aaPos}{aaAlt}"
        //  For MNV the g. part looks like ${gdnaPosStart}_${gdnaPosEnd}del${ref}ins${alt}
        //  For any type of indel there will be an _ in the range of the gdna coordinates.
        String[] chromosomeAndGDNA = (coordinateField.split("/")[0]).split(":");

        // Remove "chr" from the chromosome
        recordBuilder.chromosome(chromosomeAndGDNA[0].substring(3));

        // Remove "g." from the gdna annotation, and also remove potential parentheses
        String gdna = chromosomeAndGDNA[1].substring(2).replaceAll("[()]", "");

        if (gdna.contains(HGVS_RANGE_INDICATOR)) {
            String[] gdnaParts = gdna.split(HGVS_RANGE_INDICATOR);
            int position = Integer.parseInt(gdnaParts[0]);
            recordBuilder.gdnaPosition(position);

            TransvarAnnotation annotation;
            if (gdna.contains(HGVS_INSERTION) || gdna.contains(HGVS_DELETION)) {
                annotation = annotationForInsertionDeletion(position, gdnaParts[1], messageField);
            } else if (isFrameshift(messageField)) {
                annotation = annotationForFrameshift(coordinateField);
            } else {
                annotation = annotationForDuplication(position, gdnaParts[1]);
            }

            return recordBuilder.annotation(annotation).build();
        } else {
            return createForSNV(recordBuilder, gdna, messageField);
        }
    }

    @NotNull
    private static TransvarAnnotation annotationForFrameshift(@NotNull String coordinateField) {
        return ImmutableTransvarFrameshift.builder().isFrameshiftInsideStartCodon(coordinateField.contains("p.M1fs")).build();
    }

    @NotNull
    private static TransvarAnnotation annotationForInsertionDeletion(long position, @NotNull String delInsPart,
            @NotNull String messageField) {
        int delStart = delInsPart.indexOf(HGVS_DELETION);
        int insStart = delInsPart.indexOf(HGVS_INSERTION);

        assert delStart >= 0 || insStart >= 0;

        if (insStart >= 0) {
            // This should end in something like 'insG'
            String insertedBases = delInsPart.substring(insStart + HGVS_INSERTION.length());

            if (delStart >= 0) {
                if (delStart + HGVS_DELETION.length() == insStart) {
                    // This should look like '123delinsGGT' and is a complex insertion + deletion rather than a simple MNV
                    int deletedBaseCount = 1 + (int) (Long.parseLong(delInsPart.substring(0, delStart)) - position);
                    return ImmutableTransvarComplexInsertDelete.builder()
                            .deletedBaseCount(deletedBaseCount)
                            .insertedSequence(insertedBases)
                            .candidateAlternativeCodons(extractCandidateAlternativeCodonsFromMessageField(messageField))
                            .build();
                } else {
                    // This should look like '123delCinsG'
                    String deletedBases = delInsPart.substring(delStart + HGVS_DELETION.length(), insStart);
                    assert deletedBases.length() == insertedBases.length();

                    return ImmutableTransvarSnvMnv.builder()
                            .gdnaRef(deletedBases)
                            .gdnaAlt(insertedBases)
                            .referenceCodon(extractReferenceCodonFromMessageField(messageField))
                            .candidateCodons(extractCandidateCodonsFromMessageField(messageField))
                            .build();
                }
            } else {
                // This should look like '123insC'
                return ImmutableTransvarInsertion.builder()
                        .insertedBases(insertedBases)
                        .leftAlignedGDNAPosition(extractLeftAlignedGDNAPositionFromMessageField(messageField))
                        .build();
            }
        } else {
            // This should look like '123delC' or '123del97' in case of very long dels.
            String deletedBases = delInsPart.substring(delStart + HGVS_DELETION.length());
            int deletedBaseCount = isInteger(deletedBases) ? Integer.parseInt(deletedBases) : deletedBases.length();
            return ImmutableTransvarDeletion.builder()
                    .deletedBaseCount(deletedBaseCount)
                    .leftAlignedGDNAPosition(extractLeftAlignedGDNAPositionFromMessageField(messageField))
                    .build();
        }
    }

    @NotNull
    private static TransvarAnnotation annotationForDuplication(long position, @NotNull String dupPart) {
        int duplicatedBaseCount;
        if (dupPart.contains(HGVS_DUPLICATION)) {
            duplicatedBaseCount = 1 + (int) (Long.parseLong(dupPart.substring(0, dupPart.indexOf(HGVS_DUPLICATION))) - position);
        } else if (isLong(dupPart)) {
            duplicatedBaseCount = 1 + (int) (Long.parseLong(dupPart) - position);
        } else {
            throw new IllegalStateException("Cannot process duplication for gDNA: " + dupPart);
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

        return recordBuilder.gdnaPosition(Integer.parseInt(gdnaPos.toString())).annotation(snvAnnotation).build();
    }

    private static boolean isFrameshift(@NotNull String messageField) {
        String csqn = extractOptionalValueFromMessageField(messageField, "CSQN");
        return csqn != null && csqn.equals(CSQN_FRAMESHIFT);
    }

    private static int extractLeftAlignedGDNAPositionFromMessageField(@NotNull String messageField) {
        // Looks like g.139399409_139399411delCAC
        String leftAlignedGDNA = extractValueFromMessageField(messageField, "left_align_gDNA");
        return Integer.parseInt(leftAlignedGDNA.substring(2).split(HGVS_RANGE_INDICATOR)[0]);
    }

    @NotNull
    private static String extractReferenceCodonFromMessageField(@NotNull String messageField) {
        return extractValueFromMessageField(messageField, "reference_codon");
    }

    @NotNull
    private static List<String> extractCandidateCodonsFromMessageField(@NotNull String messageField) {
        String fieldValue = extractValueFromMessageField(messageField, "candidate_codons");

        return Arrays.asList(fieldValue.split(","));
    }

    @NotNull
    private static List<String> extractCandidateAlternativeCodonsFromMessageField(@NotNull String messageField) {
        String fieldValue = extractOptionalValueFromMessageField(messageField, "candidate_alternative_sequence");

        return fieldValue != null ? Arrays.asList(fieldValue.split("/")) : Lists.newArrayList();
    }

    private static boolean variantSpanMultipleExons(@NotNull String locationField) {
        // This looks like:
        //  - "inside_[cds_in_exons_[1,2]]" for SNV
        //  - "from_[cds_in_exon_6]_to_[cds_in_exon_7]" for inframes.
        return locationField.contains("cds_in_exons") || (locationField.contains("from") && locationField.contains("to"));
    }

    @NotNull
    private static String extractValueFromMessageField(@NotNull String messageField, @NotNull String fieldName) {
        String value = extractOptionalValueFromMessageField(messageField, fieldName);

        if (value == null) {
            throw new IllegalStateException("No '" + fieldName + "' found in message field: " + messageField);
        }

        return value;
    }

    @Nullable
    private static String extractOptionalValueFromMessageField(@NotNull String messageField, @NotNull String fieldName) {
        String[] infoFields = messageField.split(";");

        for (String infoField : infoFields) {
            if (infoField.contains(fieldName)) {
                return infoField.split("=")[1];
            }
        }

        return null;
    }

    private static boolean isLong(@NotNull String value) {
        try {
            Long.parseLong(value);
            return true;
        } catch (NumberFormatException exp) {
            return false;
        }
    }

    private static boolean isInteger(@NotNull String value) {
        try {
            Integer.parseInt(value);
            return true;
        } catch (NumberFormatException exp) {
            return false;
        }
    }
}
