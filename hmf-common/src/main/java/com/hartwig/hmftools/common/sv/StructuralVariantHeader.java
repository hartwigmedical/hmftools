package com.hartwig.hmftools.common.sv;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF_DESC;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN_CHANGE;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN_CHANGE_DESC;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN_DESC;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_JUNCTION_COPY_NUMBER;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_JUNCTION_COPY_NUMBER_DESC;

import static htsjdk.variant.vcf.VCFHeaderLineCount.UNBOUNDED;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class StructuralVariantHeader
{
    private static final String RECOVERED_DESC = "Entry has been recovered";
    private static final String RECOVERY_FILTER_DESC = "Filter before recovery";
    private static final String RECOVERY_METHOD_DESC =
            "Method used to recover, one of [UNBALANCED_SV_START, UNBALANCED_SV_END, UNSUPPORTED_BREAKEND_START, UNSUPPORTED_BREAKEND_END]";
    private static final String INFERRED_DESC = "Breakend inferred from copy number transition";
    private static final String IMPRECISE_DESC = "Imprecise structural variation";
    private static final String CIPOS_DESC = "Confidence interval around POS for imprecise variants";
    private static final String SVTYPE_DESC = "Type of structural variant";

    @NotNull
    public static VCFHeader generateHeader(@NotNull final String purpleVersion, @NotNull final VCFHeader template)
    {
        final VCFHeader outputVCFHeader = new VCFHeader(template.getMetaDataInInputOrder(), template.getGenotypeSamples());
        outputVCFHeader.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));

        outputVCFHeader.addMetaDataLine(VCFStandardHeaderLines.getFormatLine("GT"));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(StructuralVariantFactory.RECOVERED,
                0,
                VCFHeaderLineType.Flag,
                RECOVERED_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(StructuralVariantFactory.INFERRED, 0, VCFHeaderLineType.Flag, INFERRED_DESC));
        outputVCFHeader.addMetaDataLine(new VCFFilterHeaderLine(StructuralVariantFactory.INFERRED, INFERRED_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(StructuralVariantFactory.IMPRECISE,
                0,
                VCFHeaderLineType.Flag,
                IMPRECISE_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(StructuralVariantFactory.CIPOS, 2, VCFHeaderLineType.Integer, CIPOS_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(StructuralVariantFactory.SVTYPE, 1, VCFHeaderLineType.String, SVTYPE_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF, UNBOUNDED, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN, UNBOUNDED, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(StructuralVariantFactory.RECOVERY_METHOD, 1, VCFHeaderLineType.String, RECOVERY_METHOD_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(StructuralVariantFactory.RECOVERY_FILTER, UNBOUNDED, VCFHeaderLineType.String, RECOVERY_FILTER_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_JUNCTION_COPY_NUMBER, 1, VCFHeaderLineType.Float,
                PURPLE_JUNCTION_COPY_NUMBER_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN_CHANGE,
                UNBOUNDED,
                VCFHeaderLineType.Float,
                PURPLE_CN_CHANGE_DESC));
        return outputVCFHeader;
    }

}
