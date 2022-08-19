package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public final class VariantVcfTags
{
    // common
    public static final String PASS = "PASS";

    public static final String REPORTED_FLAG = "REPORTED";
    public static final String REPORTED_DESC = "Variant is reported in the driver catalog";

    // set by Sage
    public static final String TIER = "TIER";
    private static final String TIER_DESCRIPTION = "Tier: [HOTSPOT, PANEL, HIGH_CONFIDENCE, LOW_CONFIDENCE]";

    public static final String LOCAL_PHASE_SET = "LPS";
    private static final String PHASE_DESCRIPTION = "Local Phase Set";

    // set by Pave
    public static final String GNOMAD_FREQ = "GND_FREQ";
    public static final String GNOMAD_FREQ_DESC = "Gnomad variant frequency";

    // set by Purple
    public static final String PURPLE_AF_INFO = "PURPLE_AF";
    private static final String PURPLE_AF_DESC = "Purity adjusted variant allelic frequency";

    public static final String PURPLE_CN_INFO = "PURPLE_CN";
    private static final String PURPLE_CN_DESC = "Purity adjusted copy number surrounding variant location";

    public static final String PURPLE_BIALLELIC_FLAG = "BIALLELIC";
    private static final String PURPLE_BIALLELIC_DESC = "Variant is biallelic";

    public static final String PURPLE_VARIANT_CN_INFO = "PURPLE_VCN";
    private static final String PURPLE_VARIANT_CN_DESC = "Purity adjusted variant copy number";

    public static final String PURPLE_GERMLINE_INFO = "PURPLE_GERMLINE";
    private static final String PURPLE_GERMLINE_DESC = "Germline classification surrounding variant location";

    public static final String PURPLE_MINOR_ALLELE_CN_INFO = "PURPLE_MACN";
    private static final String PURPLE_MINOR_ALLELE_PLOIDY_DESC = "Purity adjusted minor allele ploidy surrounding variant location";

    private VariantVcfTags() { }

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

    @NotNull
    public static VCFHeader addSageMetaData(@NotNull final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(TIER, 1, VCFHeaderLineType.String, TIER_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(LOCAL_PHASE_SET, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, PHASE_DESCRIPTION));
        return header;
    }

    @NotNull
    public static VCFHeader purpleGermlineHeader(@NotNull final String purpleVersion, @NotNull final VCFHeader template)
    {
        template.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF_INFO, 1, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_VARIANT_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_VARIANT_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_MINOR_ALLELE_CN_INFO, 1,  VCFHeaderLineType.Float, PURPLE_MINOR_ALLELE_PLOIDY_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_BIALLELIC_FLAG, 0, VCFHeaderLineType.Flag, PURPLE_BIALLELIC_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));

        return template;
    }

    @NotNull
    public static VCFHeader purpleSomaticHeader(@NotNull final String purpleVersion, @NotNull final VCFHeader template)
    {
        template.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF_INFO, 1, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_VARIANT_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_VARIANT_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_MINOR_ALLELE_CN_INFO, 1,  VCFHeaderLineType.Float, PURPLE_MINOR_ALLELE_PLOIDY_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_GERMLINE_INFO, 1, VCFHeaderLineType.String, PURPLE_GERMLINE_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_BIALLELIC_FLAG, 0, VCFHeaderLineType.Flag, PURPLE_BIALLELIC_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));

        return template;
    }

}
