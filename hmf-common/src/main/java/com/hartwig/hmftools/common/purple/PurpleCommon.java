package com.hartwig.hmftools.common.purple;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

import org.jetbrains.annotations.NotNull;

public final class PurpleCommon
{
    public static final String PURPLE_SOMATIC_VCF_SUFFIX = ".purple.somatic.vcf.gz";
    public static final String PURPLE_GERMLINE_VCF_SUFFIX = ".purple.germline.vcf.gz";
    public static final String PURPLE_SV_VCF_SUFFIX = ".purple.sv.vcf.gz";
    public static final String PURPLE_SV_GERMLINE_VCF_SUFFIX = ".purple.sv.germline.vcf.gz";
    public static final String PURPLE_PURITY_SUFFIX = ".purple.purity.tsv";
    public static final String PURPLE_QC_SUFFIX = ".purple.qc";

    public static String purpleSomaticVcfFile(final String purpleDir, final String sampleId)
    {
        return FileWriterUtils.checkAddDirSeparator(purpleDir) + sampleId + PURPLE_SOMATIC_VCF_SUFFIX;
    }

    public static String purpleGermlineVcfFile(final String purpleDir, final String sampleId)
    {
        return FileWriterUtils.checkAddDirSeparator(purpleDir) + sampleId + PURPLE_GERMLINE_VCF_SUFFIX;
    }

    public static String purpleSomaticSvFile(final String purpleDir, final String sampleId)
    {
        return FileWriterUtils.checkAddDirSeparator(purpleDir) + sampleId + PURPLE_SV_VCF_SUFFIX;
    }

    public static String purpleGermlineSvFile(final String purpleDir, final String sampleId)
    {
        return FileWriterUtils.checkAddDirSeparator(purpleDir) + sampleId + PURPLE_SV_GERMLINE_VCF_SUFFIX;
    }

    public static String purplePurityFile(final String purpleDir, final String sampleId)
    {
        return FileWriterUtils.checkAddDirSeparator(purpleDir) + sampleId + PURPLE_PURITY_SUFFIX;
    }

    public static String purpleQcFile(final String purpleDir, final String sampleId)
    {
        return FileWriterUtils.checkAddDirSeparator(purpleDir) + sampleId + PURPLE_QC_SUFFIX;
    }

    public static DecimalFormat decimalFormat(@NotNull String format)
    {
        return new DecimalFormat(format, DecimalFormatSymbols.getInstance(Locale.ENGLISH));
    }
}
