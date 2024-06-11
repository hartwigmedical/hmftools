package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_DESC;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public final class PurpleVcfTags
{
    public static final String PURPLE_JUNCTION_COPY_NUMBER = "PURPLE_JCN";
    public static final String PURPLE_JUNCTION_COPY_NUMBER_DESC = "Purity adjusted copy number of variant junction";

    public static final String PURPLE_PLOIDY_INFO = "PURPLE_PLOIDY";

    public static final String PURPLE_CN_CHANGE = "PURPLE_CN_CHANGE";
    public static final String PURPLE_CN_CHANGE_DESC = "Purity adjusted change in copy number at each breakend";

    public static final String PURPLE_AF = "PURPLE_AF";
    public static final String PURPLE_AF_DESC = "Purity adjusted variant allelic frequency";

    public static final String PURPLE_CN = "PURPLE_CN";
    public static final String PURPLE_CN_DESC = "Purity adjusted copy number surrounding variant location";

    public static final String PURPLE_BIALLELIC_PROB = "BIALLELIC_PROB";
    public static final String PURPLE_BIALLELIC_PROB_DESC = "Biallelic probability";
    
    public static final String PURPLE_BIALLELIC_FLAG = "BIALLELIC";
    public static final String PURPLE_BIALLELIC_DESC = "If the probability this variant is biallelic equals or is greater than 0.50";

    public static final String PURPLE_VARIANT_CN = "PURPLE_VCN";
    private static final String PURPLE_VARIANT_CN_DESC = "Purity adjusted variant copy number";

    public static final String PURPLE_GERMLINE_INFO = "PURPLE_GERMLINE";
    private static final String PURPLE_GERMLINE_DESC = "Germline classification surrounding variant location";

    public static final String PURPLE_MINOR_ALLELE_CN_INFO = "PURPLE_MACN";
    private static final String PURPLE_MINOR_ALLELE_PLOIDY_DESC = "Purity adjusted minor allele ploidy surrounding variant location";

    public static final String SUBCLONAL_LIKELIHOOD_FLAG = "SUBCL";
    public static final String SUBCLONAL_LIKELIHOOD_FLAG_DESCRIPTION = "Non-zero subclonal likelihood";

    public static final String KATAEGIS_FLAG = "KT";
    public static final String KATAEGIS_FLAG_DESCRIPTION = "Forward/reverse kataegis id";

    public static final String PANEL_GERMLINE_VAF_DISTANCE = "VAF_DIS_MIN";
    public static final String PANEL_GERMLINE_VAF_DISTANCE_DESC = "Panel germline VAF distance, somatic VAF distance";

    public static final String PANEL_SOMATIC_LIKELIHOOD = "SOM_LH";
    public static final String PANEL_SOMATIC_LIKELIHOOD_DESC = "Panel somatic likelihood [HIGH,MEDIUM,LOW]";

    public static final String REPORTABLE_TRANSCRIPTS = "REPORTABLE_TRANSCRIPTS";
    public static final String REPORTABLE_TRANSCRIPTS_DESC = "List of reportable transcript when non-canonical are reportable";
    public static final String REPORTABLE_TRANSCRIPTS_DELIM = "|";


    public static void addGermlineHeader(@NotNull final String purpleVersion, @NotNull final VCFHeader header)
    {
        header.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));
        header.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF, 1, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN, 1, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_VARIANT_CN, 1, VCFHeaderLineType.Float, PURPLE_VARIANT_CN_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_MINOR_ALLELE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_MINOR_ALLELE_PLOIDY_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_BIALLELIC_FLAG, 0, VCFHeaderLineType.Flag, PURPLE_BIALLELIC_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(REPORTABLE_TRANSCRIPTS, 1, VCFHeaderLineType.String, REPORTABLE_TRANSCRIPTS_DESC));
    }

    @NotNull
    public static void addSomaticHeader(@NotNull final String purpleVersion, @NotNull final VCFHeader header)
    {
        header.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));
        header.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF, 1, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN, 1, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_VARIANT_CN, 1, VCFHeaderLineType.Float, PURPLE_VARIANT_CN_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_MINOR_ALLELE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_MINOR_ALLELE_PLOIDY_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_GERMLINE_INFO, 1, VCFHeaderLineType.String, PURPLE_GERMLINE_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_BIALLELIC_FLAG, 0, VCFHeaderLineType.Flag, PURPLE_BIALLELIC_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_BIALLELIC_PROB, 1, VCFHeaderLineType.Float, PURPLE_BIALLELIC_PROB_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(REPORTED_FLAG, 0, VCFHeaderLineType.Flag, REPORTED_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(REPORTABLE_TRANSCRIPTS, 1, VCFHeaderLineType.String, REPORTABLE_TRANSCRIPTS_DESC));
    }
}
