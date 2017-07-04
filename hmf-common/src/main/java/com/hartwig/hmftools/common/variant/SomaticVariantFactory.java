package com.hartwig.hmftools.common.variant;

import com.google.common.collect.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.common.variant.VariantFactory.*;

public final class SomaticVariantFactory {

    private static final Logger LOGGER = LogManager.getLogger(SomaticVariantFactory.class);

    private static final String DBSNP_IDENTIFIER = "rs";
    private static final String COSMIC_IDENTIFIER = "COSM";
    private static final String ID_SEPARATOR = ";";

    private static final String INFO_FIELD_SEPARATOR = ";";
    private static final String CALLER_ALGO_IDENTIFIER = "set=";
    private static final String CALLER_ALGO_START = "=";
    private static final String CALLER_ALGO_SEPARATOR = "-";
    private static final String CALLER_FILTERED_IDENTIFIER = "filterIn";
    private static final String CALLER_INTERSECTION_IDENTIFIER = "Intersection";

    private static final int SAMPLE_DATA_COLUMN = 9;
    private static final int SAMPLE_DATA_ALLELE_FREQUENCY_COLUMN = 1;

    private SomaticVariantFactory() {
    }

    @NotNull
    public static String sampleFromHeaderLine(@NotNull final String headerLine) {
        return VariantFactory.sampleFromHeaderLine(headerLine, SAMPLE_DATA_COLUMN);
    }

    @NotNull
    public static String toVCFLine(@NotNull final SomaticVariant variant) {
        final String[] values = variant.originalVCFLine().split(VCF_COLUMN_SEPARATOR);
        // KODU: Currently assume only the filter value could change, but this is potentially wrong!
        if (!values[FILTER_COLUMN].trim().equals(variant.filter())) {
            values[FILTER_COLUMN] = variant.filter();
        }
        String reassembledLine = values[0];
        for (int i = 1; i < values.length; i++) {
            reassembledLine += (VCF_COLUMN_SEPARATOR + values[i]);
        }
        return reassembledLine;
    }

    @NotNull
    public static SomaticVariant fromVCFLine(@NotNull final String line) {
        final SomaticVariant.Builder builder = SomaticVariant.Builder.fromVCF(line);
        final String[] values = line.split(VCF_COLUMN_SEPARATOR);
        VariantFactory.withLine(builder, values);

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

        final String sampleData = values[SAMPLE_DATA_COLUMN].trim();
        final AlleleFrequencyData alleleFrequencyData = extractAlleleFrequencyData(sampleData);
        
        if (null == alleleFrequencyData) {
            LOGGER.warn("Could not parse allele frequencies from " + sampleData);
        } else {
            builder.totalReadCount(alleleFrequencyData.totalReadCount());
            builder.alleleReadCount(alleleFrequencyData.alleleReadCount());
        }

        return builder.build();
    }
    
    @NotNull
    public static SomaticTruthSetVariant fromTruthSetVCFLine(@NotNull final String line) {
        final SomaticTruthSetVariant.Builder builder = SomaticTruthSetVariant.Builder.fromVCF(line);
        final String[] values = line.split(VCF_COLUMN_SEPARATOR);
        VariantFactory.withLine(builder, values);
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
    private static AlleleFrequencyData extractAlleleFrequencyData(@NotNull final String sampleData) {
        final String[] sampleFields = VariantFactoryFunctions.splitSampleDataFields(sampleData);
        if (sampleFields.length < 2) {
            return null;
        }

        return VariantFactoryFunctions.determineAlleleFrequencies(sampleFields[SAMPLE_DATA_ALLELE_FREQUENCY_COLUMN]);
    }
}
