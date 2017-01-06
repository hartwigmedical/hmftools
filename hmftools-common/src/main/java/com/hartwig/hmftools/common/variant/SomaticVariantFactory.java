package com.hartwig.hmftools.common.variant;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class SomaticVariantFactory {

    private static final Logger LOGGER = LogManager.getLogger(SomaticVariantFactory.class);
    private static final String VCF_COLUMN_SEPARATOR = "\t";

    private static final int CHROMOSOME_COLUMN = 0;
    private static final int POSITION_COLUMN = 1;

    private static final int ID_COLUMN = 2;
    private static final String DBSNP_IDENTIFIER = "rs";
    private static final String COSMIC_IDENTIFIER = "COSM";

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
    private static final String MISSENSE_IDENTIFIER = "missense_variant";

    private static final int SAMPLE_COLUMN = 9;
    private static final String SAMPLE_FIELD_SEPARATOR = ":";
    private static final int AF_COLUMN_INDEX = 1;
    private static final String AF_FIELD_SEPARATOR = ",";

    private SomaticVariantFactory() {
    }

    @NotNull
    public static String sampleFromHeaderLine(@NotNull final String headerLine) {
        final String[] values = headerLine.split(VCF_COLUMN_SEPARATOR);
        return values[SAMPLE_COLUMN];
    }

    @NotNull
    public static SomaticVariant fromVCFLine(@NotNull final String line) {
        final String[] values = line.split(VCF_COLUMN_SEPARATOR);

        final VariantType type = VariantExtractorFunctions.determineVariantType(values[REF_COLUMN].trim(),
                values[ALT_COLUMN].trim());
        final SomaticVariant.Builder builder = new SomaticVariant.Builder(type);

        builder.filter(values[FILTER_COLUMN].trim());

        final String info = values[INFO_COLUMN].trim();
        builder.callers(extractCallers(info));

        final String sampleInfo = values[SAMPLE_COLUMN].trim();
        final double alleleFrequency = calcAlleleFrequency(sampleInfo);
        if (Double.isNaN(alleleFrequency)) {
            LOGGER.warn("Could not parse alleleFrequency from " + line);
        }
        builder.alleleFrequency(alleleFrequency);

        builder.chromosome(values[CHROMOSOME_COLUMN].trim());
        builder.position(Long.valueOf(values[POSITION_COLUMN].trim()));

        final String id = values[ID_COLUMN];
        builder.isDBSNP(id.contains(DBSNP_IDENTIFIER));
        builder.isCOSMIC(id.contains(COSMIC_IDENTIFIER));
        builder.isMissense(info.contains(MISSENSE_IDENTIFIER));

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
        List<String> finalCallers = Lists.newArrayList();
        if (allCallers.length > 0 && allCallers[0].equals(CALLER_INTERSECTION_IDENTIFIER)) {
            finalCallers.addAll(SomaticVariantConstants.ALL_CALLERS);
        } else {
            finalCallers.addAll(
                    Arrays.stream(allCallers).filter(caller -> !caller.startsWith(CALLER_FILTERED_IDENTIFIER)).collect(
                            Collectors.toList()));
        }
        return finalCallers;
    }

    private static double calcAlleleFrequency(@NotNull final String sampleInfo) {
        final String[] sampleFields = sampleInfo.split(SAMPLE_FIELD_SEPARATOR);
        if (sampleFields.length < 2) {
            return Double.NaN;
        }

        final String[] afFields = sampleFields[AF_COLUMN_INDEX].split(AF_FIELD_SEPARATOR);

        if (afFields.length < 2) {
            return Double.NaN;
        }

        final int readCount = Integer.valueOf(afFields[1]);
        int totalReadCount = 0;
        for (String afField : afFields) {
            totalReadCount += Integer.valueOf(afField);
        }

        return (double) readCount / totalReadCount;
    }
}
