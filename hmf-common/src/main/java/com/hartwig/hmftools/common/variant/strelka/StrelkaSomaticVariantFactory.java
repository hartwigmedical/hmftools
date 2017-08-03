package com.hartwig.hmftools.common.variant.strelka;

import static com.hartwig.hmftools.common.variant.VariantFactory.ALT_COLUMN;
import static com.hartwig.hmftools.common.variant.VariantFactory.CHROMOSOME_COLUMN;
import static com.hartwig.hmftools.common.variant.VariantFactory.FILTER_COLUMN;
import static com.hartwig.hmftools.common.variant.VariantFactory.POSITION_COLUMN;
import static com.hartwig.hmftools.common.variant.VariantFactory.REF_COLUMN;
import static com.hartwig.hmftools.common.variant.VariantFactory.VCF_COLUMN_SEPARATOR;

import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class StrelkaSomaticVariantFactory {

    private static final Logger LOGGER = LogManager.getLogger(SomaticVariantFactory.class);

    private static final String INFO_FIELD_SEPARATOR = ";";
    public static final String FORMAT_SEPARATOR = ":";
    public static final String FORMAT_VALUES_SEPARATOR = ",";

    private static final int INFO_COLUMN = 7;
    private static final int FORMAT_COLUMN = 8;
    private static final int TUMOR_SAMPLE_DATA_COLUMN = 10;

    @NotNull
    public static StrelkaSomaticVariant fromVCFLine(@NotNull final String line) {
        final String[] values = line.split(VCF_COLUMN_SEPARATOR);
        final String ref = values[REF_COLUMN].trim();
        final String alt = values[ALT_COLUMN].trim();
        final VariantType variantType = VariantType.fromRefAlt(ref, alt);
        final String chromosome = values[CHROMOSOME_COLUMN].trim();
        final long position = Long.valueOf(values[POSITION_COLUMN].trim());
        final String filter = values[FILTER_COLUMN].trim();
        final ArrayListMultimap<String, String> tumorData = readTumorData(chromosome, position, values);
        final Map<String, String> infoData = readInfoFields(values);
        return ImmutableStrelkaSomaticVariant.of(line, variantType, chromosome, position, ref, alt, filter, infoData, tumorData);
    }

    @NotNull
    private static ArrayListMultimap<String, String> readTumorData(@NotNull final String chromosome, final long position,
            @NotNull final String[] values) {
        final ArrayListMultimap<String, String> tumorData = ArrayListMultimap.create();
        final String tumorDataString = values[TUMOR_SAMPLE_DATA_COLUMN].trim();
        final String formatString = values[FORMAT_COLUMN].trim();
        final String[] formatKeys = formatString.split(FORMAT_SEPARATOR);
        final String[] tumorDataValues = tumorDataString.split(FORMAT_SEPARATOR);
        if (formatKeys.length != tumorDataValues.length) {
            LOGGER.warn("variant at chromosome " + chromosome + "(" + position + ") has different lengths for format and tumor data");
        }
        for (int index = 0; index < formatKeys.length; index++) {
            final String key = formatKeys[index];
            final String[] keyValues = tumorDataValues[index].split(FORMAT_VALUES_SEPARATOR);
            tumorData.putAll(key, Arrays.stream(keyValues).collect(Collectors.toList()));
        }
        return tumorData;
    }

    @NotNull
    private static Map<String, String> readInfoFields(@NotNull final String[] values) {
        final Map<String, String> infoData = Maps.newHashMap();
        final String[] infoDataStrings = values[INFO_COLUMN].split(INFO_FIELD_SEPARATOR);
        for (final String infoDataString : infoDataStrings) {
            final String[] infoDataFields = infoDataString.split("=", 2);
            if (infoDataFields.length > 1) {
                infoData.put(infoDataFields[0], infoDataFields[1]);
            } else {
                infoData.put(infoDataFields[0], "");
            }
        }
        return infoData;
    }
}
