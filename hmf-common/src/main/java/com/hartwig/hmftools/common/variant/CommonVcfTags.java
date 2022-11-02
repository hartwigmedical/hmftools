package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public final class CommonVcfTags
{
    // common
    public static final String PASS = "PASS";

    public static final String REPORTED_FLAG = "REPORTED";
    public static final String REPORTED_DESC = "Variant is reported in the driver catalog";

    private CommonVcfTags() { }

    public static int getGenotypeAttributeAsInt(final Genotype genotype, final String attribute, int defaultVaue)
    {
        Object value = genotype.getExtendedAttribute(attribute);
        return value == null ? defaultVaue : Integer.parseInt(value.toString());
    }

    public static double getGenotypeAttributeAsDouble(final Genotype genotype, final String attribute, double defaultVaue)
    {
        Object value = genotype.getExtendedAttribute(attribute);
        return value == null ? defaultVaue : Double.parseDouble(value.toString());
    }



}
