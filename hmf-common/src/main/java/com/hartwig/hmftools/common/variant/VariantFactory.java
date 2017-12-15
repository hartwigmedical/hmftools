package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Deprecated
final class VariantFactory {
    static final String VCF_COLUMN_SEPARATOR = "\t";
    private static final int CHROMOSOME_COLUMN = 0;
    private static final int POSITION_COLUMN = 1;
    static final int ID_COLUMN = 2;
    private static final int REF_COLUMN = 3;
    private static final int ALT_COLUMN = 4;
    static final int FILTER_COLUMN = 6;
    static final int INFO_COLUMN = 7;

    static <T extends Variant> void withLine(@NotNull final VariantBuilder<T> builder, @NotNull final String[] values) {
        final String ref = values[REF_COLUMN].trim();
        final String alt = values[ALT_COLUMN].trim();

        builder.ref(ref);
        builder.alt(alt);
        builder.type(VariantType.fromRefAlt(ref, alt));

        builder.chromosome(values[CHROMOSOME_COLUMN].trim());
        builder.position(Long.valueOf(values[POSITION_COLUMN].trim()));
        builder.filter(values[FILTER_COLUMN].trim());
    }

    @Nullable
    static String sampleFromHeaderLine(@NotNull final String headerLine, int sampleColumn) {
        final String[] values = headerLine.split(VCF_COLUMN_SEPARATOR);
        if (values.length <= sampleColumn) {
            return null;
        } else {
            return values[sampleColumn];
        }
    }
}
