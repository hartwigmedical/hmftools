package com.hartwig.hmftools.common.variant;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SomaticVariantFactory {

    private static final Logger LOGGER = LogManager.getLogger(SomaticVariantFactory.class);
    private static final String VCF_COLUMN_SEPARATOR = "\t";

    private static final int CHROMOSOME_COLUMN = 0;
    private static final int POSITION_COLUMN = 1;

    private static final int ID_COLUMN = 2;
    private static final String DBSNP_IDENTIFIER = "rs";
    private static final String COSMIC_IDENTIFIER = "COSM";
    private static final String ID_SEPARATOR = ";";

    private static final int REF_COLUMN = 3;
    private static final int ALT_COLUMN = 4;

    private static final int FILTER_COLUMN = 6;

    private static final int INFO_COLUMN = 7;
    private static final String INFO_FIELD_SEPARATOR = ";";
    private static final String CALLER_ALGO_IDENTIFIER = "set=";
    private static final String CALLER_ALGO_START = "=";
    private static final String CALLER_ALGO_SEPARATOR = "-";
    private static final String CALLER_FILTERED_IDENTIFIER = "filterIn";
    private static final String CALLER_INTERSECTION_IDENTIFIER = "Intersection";

    private static final int SAMPLE_COLUMN = 9;
    private static final String SAMPLE_FIELD_SEPARATOR = ":";
    private static final int AF_COLUMN_INDEX = 1;
    private static final String AF_FIELD_SEPARATOR = ",";

    private SomaticVariantFactory() {
    }

    @NotNull
    public static String sampleFromHeaderLine(@NotNull final String headerLine) {
        final String[] values = headerLine.split(VCF_COLUMN_SEPARATOR);
        final String sample = values[SAMPLE_COLUMN];
        // KODU: In v1.7, the sample would contain the whole path of the VCF.
        if (sample.contains(File.separator)) {
            final String[] parts = sample.split(File.separator);
            final String[] subParts = parts[parts.length - 1].split("_");
            // KODU: Assume last part starts with "R_T"
            return subParts[1];
        } else {
            return sample;
        }
    }

    @NotNull
    public static String toVCFLine(@NotNull final SomaticVariant variant) {
        return toVCFLine(variant, null);
    }

    @NotNull
    public static String toVCFLine(@NotNull final SomaticVariant variant, @Nullable String newInfoValue) {
        final String originalVCFLine = variant.originalVCFLine();
        if (newInfoValue != null) {
            final String[] values = originalVCFLine.split(VCF_COLUMN_SEPARATOR);
            String reconstructedVCFLine = values[0];
            for (int i = 1; i < values.length; i++) {
                final String value = i == INFO_COLUMN ? newInfoValue : values[i];
                reconstructedVCFLine += (VCF_COLUMN_SEPARATOR + value);
            }

            return reconstructedVCFLine;
        } else {
            return originalVCFLine;
        }
    }

    @NotNull
    public static SomaticVariant fromVCFLine(@NotNull final String line) {
        final String[] values = line.split(VCF_COLUMN_SEPARATOR);

        final VariantType type = VariantExtractorFunctions.determineVariantType(values[REF_COLUMN].trim(),
                values[ALT_COLUMN].trim());
        final SomaticVariant.Builder builder = SomaticVariant.Builder.fromVCF(line, type);

        builder.chromosome(values[CHROMOSOME_COLUMN].trim());
        builder.position(Long.valueOf(values[POSITION_COLUMN].trim()));
        builder.filter(values[FILTER_COLUMN].trim());

        final String idValue = values[ID_COLUMN].trim();
        if (!idValue.isEmpty()) {
            final String[] ids = idValue.split(ID_SEPARATOR);
            for (final String id : ids) {
                if (id.contains(DBSNP_IDENTIFIER)) {
                    builder.dnsnpID(id);
                } else if (id.contains(COSMIC_IDENTIFIER)) {
                    builder.cosmicID(id);
                }
            }
        }

        final String info = values[INFO_COLUMN].trim();
        builder.callers(extractCallers(info));
        builder.annotations(VariantAnnotationFactory.fromVCFInfoField(info));

        // KODU: For testing
        final List<VariantAnnotation> annotations = VariantAnnotationFactory.fromVCFInfoField(info);
        boolean annotationHasMissense = false;
        for (final VariantAnnotation annotation : annotations) {
            annotationHasMissense =
                    annotationHasMissense || annotation.consequences().contains(VariantConsequence.MISSENSE_VARIANT);
        }
        final boolean infoHasMissense = info.contains("missense");
        if (annotationHasMissense != infoHasMissense) {
            LOGGER.warn("Mismatch between annotated missense and normal missense: " + info);
        }
        // KODU: End testing

        final String sampleInfo = values[SAMPLE_COLUMN].trim();
        final ReadCount readCounts = extractReadCounts(sampleInfo);
        if (readCounts == null) {
            LOGGER.warn("Could not parse read counts from " + line);
        } else {
            builder.alleleFrequency(readCounts.alleleFrequency());
            builder.readDepth(readCounts.readDepth());
        }

        return builder.build();
    }

    @NotNull
    private static List<String> extractCallers(@NotNull final String info) {
        final Optional<String> setValue = Arrays.stream(info.split(INFO_FIELD_SEPARATOR)).filter(
                infoLine -> infoLine.contains(CALLER_ALGO_IDENTIFIER)).map(
                infoLine -> infoLine.substring(infoLine.indexOf(CALLER_ALGO_START) + 1,
                        infoLine.length())).findFirst();
        if (!setValue.isPresent()) {
            LOGGER.warn("No caller info found in info field: " + info);
            return Lists.newArrayList();
        }

        final String[] allCallers = setValue.get().split(CALLER_ALGO_SEPARATOR);
        final List<String> finalCallers = Lists.newArrayList();
        if (allCallers.length > 0 && allCallers[0].equals(CALLER_INTERSECTION_IDENTIFIER)) {
            finalCallers.addAll(SomaticVariantConstants.ALL_CALLERS);
        } else {
            finalCallers.addAll(
                    Arrays.stream(allCallers).filter(caller -> !caller.startsWith(CALLER_FILTERED_IDENTIFIER)).collect(
                            Collectors.toList()));
        }
        return finalCallers;
    }

    @Nullable
    private static ReadCount extractReadCounts(@NotNull final String sampleData) {
        final String[] sampleFields = sampleData.split(SAMPLE_FIELD_SEPARATOR);
        if (sampleFields.length < 2) {
            return null;
        }

        final String[] afFields = sampleFields[AF_COLUMN_INDEX].split(AF_FIELD_SEPARATOR);

        if (afFields.length < 2) {
            return null;
        }

        final int readCount = Integer.valueOf(afFields[1]);
        int totalReadCount = 0;
        for (final String afField : afFields) {
            totalReadCount += Integer.valueOf(afField);
        }

        return new ReadCount(readCount, totalReadCount);
    }

    private static class ReadCount {
        private final int alleleReadCount;
        private final int totalReadCount;

        private ReadCount(final int alleleReadCount, final int totalReadCount) {
            this.alleleReadCount = alleleReadCount;
            this.totalReadCount = totalReadCount;
        }

        private double alleleFrequency() {
            return (double) alleleReadCount / totalReadCount;
        }

        private int readDepth() {
            return totalReadCount;
        }
    }
}
