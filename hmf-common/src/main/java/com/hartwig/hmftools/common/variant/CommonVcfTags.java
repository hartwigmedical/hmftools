package com.hartwig.hmftools.common.variant;

import htsjdk.variant.variantcontext.Genotype;

public final class CommonVcfTags
{
    // common
    public static final String PASS = "PASS";

    public static final String QUAL = "QUAL";
    public static final String QUAL_DESC = "Variant quality";

    public static final String REPORTED_FLAG = "REPORTED";
    public static final String REPORTED_DESC = "Variant is reported in the driver catalog";

    public static int getGenotypeAttributeAsInt(final Genotype genotype, final String attribute, int defaultValue)
    {
        Object value = genotype.getExtendedAttribute(attribute);
        return value == null ? defaultValue : Integer.parseInt(value.toString());
    }

    public static double getGenotypeAttributeAsDouble(final Genotype genotype, final String attribute, double defaultValue)
    {
        Object value = genotype.getExtendedAttribute(attribute);
        return value == null ? defaultValue : Double.parseDouble(value.toString());
    }
}
