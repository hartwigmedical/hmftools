package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.SageVcfTags.LIST_SEPARATOR;

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
        if(genotype == null)
            return defaultValue;

        Object value = genotype.getExtendedAttribute(attribute);
        return value == null ? defaultValue : Integer.parseInt(value.toString());
    }

    public static double getGenotypeAttributeAsDouble(final Genotype genotype, final String attribute, double defaultValue)
    {
        if(genotype == null)
            return defaultValue;

        Object value = genotype.getExtendedAttribute(attribute);
        return value == null ? defaultValue : Double.parseDouble(value.toString());
    }

    public static int[] parseIntegerList(final Genotype genotype, final String vcfTag)
    {
        final String[] stringValues = genotype.getExtendedAttribute(vcfTag, 0).toString().split(LIST_SEPARATOR, -1);
        int[] values = new int[stringValues.length];

        for(int i = 0; i < stringValues.length; ++i)
        {
            values[i] = Integer.parseInt(stringValues[i]);
        }

        return values;
    }
}
