package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VariantHeader {

    public static final String PASS = "PASS";
    public static final String REPORTED_FLAG = "REPORTED";
    public static final String PURPLE_AF_INFO = "PURPLE_AF";
    public static final String PURPLE_CN_INFO = "PURPLE_CN";
    public static final String PURPLE_BIALLELIC_FLAG = "BIALLELIC";
    public static final String PURPLE_VARIANT_CN_INFO = "PURPLE_VCN";
    public static final String PURPLE_GERMLINE_INFO = "PURPLE_GERMLINE";
    public static final String PURPLE_MINOR_ALLELE_CN_INFO = "PURPLE_MACN";

    // KEEP FOR BACKWARDS COMPATIBILITY
    public static final String PURPLE_VARIANT_PLOIDY_INFO = "PURPLE_PLOIDY";
    public static final String PURPLE_MINOR_ALLELE_PLOIDY_INFO = "PURPLE_MAP";

    private static final String PURPLE_CN_DESC = "Purity adjusted copy number surrounding variant location";
    private static final String PURPLE_MINOR_ALLELE_PLOIDY_DESC = "Purity adjusted minor allele ploidy surrounding variant location";
    private static final String PURPLE_GERMLINE_DESC = "Germline classification surrounding variant location";
    private static final String PURPLE_AF_DESC = "Purity adjusted variant allelic frequency";
    private static final String PURPLE_PLOIDY_DESC = "Purity adjusted variant copy number";
    private static final String PURPLE_BIALLELIC_DESC = "Variant is biallelic";
    public static final String REPORTED_DESC = "Variant is reported in the driver catalog";

    @NotNull
    public static VCFHeader germlineHeader(@NotNull final String purpleVersion, @NotNull final VCFHeader template) {
        template.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF_INFO, 1, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_VARIANT_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_PLOIDY_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_MINOR_ALLELE_CN_INFO, 1,  VCFHeaderLineType.Float, PURPLE_MINOR_ALLELE_PLOIDY_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_BIALLELIC_FLAG, 0, VCFHeaderLineType.Flag, PURPLE_BIALLELIC_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));

        return template;
    }

    @NotNull
    public static VCFHeader somaticHeader(@NotNull final String purpleVersion, @NotNull final VCFHeader template) {
        template.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF_INFO, 1, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_VARIANT_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_PLOIDY_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_MINOR_ALLELE_CN_INFO, 1,  VCFHeaderLineType.Float, PURPLE_MINOR_ALLELE_PLOIDY_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_GERMLINE_INFO, 1, VCFHeaderLineType.String, PURPLE_GERMLINE_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_BIALLELIC_FLAG, 0, VCFHeaderLineType.Flag, PURPLE_BIALLELIC_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));

        return template;
    }

}
