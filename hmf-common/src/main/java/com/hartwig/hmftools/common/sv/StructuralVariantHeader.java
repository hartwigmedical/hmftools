package com.hartwig.hmftools.common.sv;

import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INFERRED;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INFERRED_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.RECOVERED;
import static com.hartwig.hmftools.common.sv.SvVcfTags.RECOVERED_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.RECOVERY_FILTER;
import static com.hartwig.hmftools.common.sv.SvVcfTags.RECOVERY_FILTER_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.RECOVERY_METHOD;
import static com.hartwig.hmftools.common.sv.SvVcfTags.RECOVERY_METHOD_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SVTYPE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SVTYPE_DESC;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.IMPRECISE;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.IMPRECISE_DESC;
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
    @NotNull
    public static VCFHeader generateHeader(@NotNull final String purpleVersion, @NotNull final VCFHeader template)
    {
        final VCFHeader outputVCFHeader = new VCFHeader(template.getMetaDataInInputOrder(), template.getGenotypeSamples());
        outputVCFHeader.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));

        outputVCFHeader.addMetaDataLine(VCFStandardHeaderLines.getFormatLine("GT"));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(RECOVERED,
                0,
                VCFHeaderLineType.Flag,
                RECOVERED_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(INFERRED, 0, VCFHeaderLineType.Flag, INFERRED_DESC));
        outputVCFHeader.addMetaDataLine(new VCFFilterHeaderLine(INFERRED, INFERRED_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(IMPRECISE,
                0,
                VCFHeaderLineType.Flag,
                IMPRECISE_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(CIPOS, 2, VCFHeaderLineType.Integer, CIPOS_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(SVTYPE, 1, VCFHeaderLineType.String, SVTYPE_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF, UNBOUNDED, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN, UNBOUNDED, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(RECOVERY_METHOD, 1, VCFHeaderLineType.String, RECOVERY_METHOD_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(RECOVERY_FILTER, UNBOUNDED, VCFHeaderLineType.String, RECOVERY_FILTER_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_JUNCTION_COPY_NUMBER, 1, VCFHeaderLineType.Float,
                PURPLE_JUNCTION_COPY_NUMBER_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN_CHANGE,
                UNBOUNDED,
                VCFHeaderLineType.Float,
                PURPLE_CN_CHANGE_DESC));
        return outputVCFHeader;
    }
}
