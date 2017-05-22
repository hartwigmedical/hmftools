package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

import java.io.File;

public class VariantFactory {
    static final String VCF_COLUMN_SEPARATOR = "\t";
    static final int CHROMOSOME_COLUMN = 0;
    static final int POSITION_COLUMN = 1;
    static final int ID_COLUMN = 2;
    static final int REF_COLUMN = 3;
    static final int ALT_COLUMN = 4;
    static final int FILTER_COLUMN = 6;
    static final int INFO_COLUMN = 7;


    static <T extends Variant> VariantBuilder<T> withLine(VariantBuilder<T> builder, @NotNull String[] values) {

        final String ref = values[REF_COLUMN].trim();
        final String alt = values[ALT_COLUMN].trim();

        builder.ref(ref);
        builder.alt(alt);
        builder.type(VariantType.fromRefAlt(ref, alt));

        builder.chromosome(values[CHROMOSOME_COLUMN].trim());
        builder.position(Long.valueOf(values[POSITION_COLUMN].trim()));
        builder.filter(values[FILTER_COLUMN].trim());

        return builder;
    }

    @NotNull
    static String sampleFromHeaderLine(@NotNull final String headerLine, int sampleColumn) {
        final String[] values = headerLine.split(VCF_COLUMN_SEPARATOR);
        final String sample = values[sampleColumn];
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
}
