package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

class VariantFactory {
    public static final String VCF_COLUMN_SEPARATOR = "\t";
    public static final int CHROMOSOME_COLUMN = 0;
    public static final int POSITION_COLUMN = 1;
    public static final int ID_COLUMN = 2;
    public static final int REF_COLUMN = 3;
    public static final int ALT_COLUMN = 4;
    public static final int FILTER_COLUMN = 6;
    public static final int INFO_COLUMN = 7;

    public static <T extends Variant> VariantBuilder<T> withLine(VariantBuilder<T> builder, @NotNull String[] values) {

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
}
