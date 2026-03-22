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

    public static final String PURPLE_PLOT_INPUT_CIRCOS = ".input.png";
    public static final String PURPLE_PLOT_FINAL_CIRCOS = ".circos.png";
    public static final String PURPLE_PLOT_COPY_NUMBER = ".copynumber.png";
    public static final String PURPLE_PLOT_PURITY_RANGE = ".purity.range.png";
    public static final String PURPLE_PLOT_SOMATIC_CLONALITY = ".somatic.clonality.png";
    public static final String PURPLE_PLOT_SOMATIC_RAINFALL = ".somatic.rainfall.png";
    public static final String PURPLE_PLOT_SOMATIC_CN = ".somatic.png";
    public static final String PURPLE_PLOT_MINOR_ALLELE = ".map.png";

    public static final double DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO = 3;
    public static final double DRIVER_AMPLIFICATION_CANDIDATE_PLOIDY_RATIO = 0.85;
    public static final double DEFAULT_DRIVER_HET_DELETION_THRESHOLD = 0.6;
    public static final double LOH_MINOR_ALLEL_CN = 0.5;

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

    public static String purplePlotFile(final String purpleDir, final String sampleId, final String plotSuffix)
    {
        return FileWriterUtils.checkAddDirSeparator(purpleDir) + sampleId + plotSuffix;
    }

    public static DecimalFormat decimalFormat(@NotNull String format)
    {
        return new DecimalFormat(format, DecimalFormatSymbols.getInstance(Locale.ENGLISH));
    }
}
