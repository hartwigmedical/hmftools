package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_DESC;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;
import static com.hartwig.hmftools.common.variant.Hotspot.HOTSPOT;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public final class PurpleVcfTags
{
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

    public static VCFHeader addGermlineHeader(@NotNull final String purpleVersion, @NotNull final VCFHeader template)
    {
        template.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF_INFO, 1, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_VARIANT_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_VARIANT_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_MINOR_ALLELE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_MINOR_ALLELE_PLOIDY_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_BIALLELIC_FLAG, 0, VCFHeaderLineType.Flag, PURPLE_BIALLELIC_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));

        return template;
    }

    @NotNull
    public static VCFHeader addSomaticHeader(@NotNull final String purpleVersion, @NotNull final VCFHeader template)
    {
        template.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF_INFO, 1, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_VARIANT_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_VARIANT_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_MINOR_ALLELE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_MINOR_ALLELE_PLOIDY_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_GERMLINE_INFO, 1, VCFHeaderLineType.String, PURPLE_GERMLINE_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_BIALLELIC_FLAG, 0, VCFHeaderLineType.Flag, PURPLE_BIALLELIC_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));

        return template;
    }
}
