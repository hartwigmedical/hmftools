package com.hartwig.hmftools.compar.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CIDER_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CIDER_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CUPPA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CUPPA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LILAC_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LILAC_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PEACH_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PEACH_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TEAL_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TEAL_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.ComparConfig.NEW_SOURCE;
import static com.hartwig.hmftools.compar.ComparConfig.REF_SOURCE;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.pipeline.reference.api.PipelineOutputStructure;
import com.hartwig.pipeline.reference.api.PipelineVersion;

public record FileSources(String source, String sampleDir, String linx, String purple, String linxGermline, String cuppa, String lilac,
                          String chord, String peach, String virus, String somaticVcf, String somaticUnfilteredVcf, String tumorFlagstat,
                          String germlineFlagstat, String tumorBamMetrics, String germlineBamMetrics, String snpGenotype, String cider,
                          String teal, PipelineVersion pipelineVersion, PipelineOutputStructure pipelineOutputStructure)
{
    private static final String SAMPLE_DIR = "sample_dir";
    private static final String SOMATIC_VCF = "somatic_vcf";
    private static final String SOMATIC_UNFILTERED_VCF = "somatic_unfiltered_vcf";
    private static final String TUMOR_FLAGSTAT = "tumor_flagstat_dir";
    private static final String GERMLINE_FLAGSTAT = "germline_flagstat_dir";
    private static final String TUMOR_BAM_METRICS = "tumor_bam_metrics_dir";
    private static final String GERMLINE_BAM_METRICS = "germline_bam_metrics_dir";
    private static final String SNP_GENOTYPE = "snp_genotype_dir";

    private static final String PIPELINE_VERSION = "pipeline_version";
    private static final String PIPELINE_OUTPUT_STRUCTURE = "pipeline_output_structure";

    private static final Map<String, String> PIPELINE_VERSION_VALUE_TRANSLATIONS = Map.of("2.0", "6.0", "5.34", "5.34", "6.0", "6.0");

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        List<String> sourceNames = Lists.newArrayList(REF_SOURCE, NEW_SOURCE);

        for(String sourceName : sourceNames)
        {
            configBuilder.addPath(
                    formSourceConfig(SAMPLE_DIR, sourceName), false,
                    formSourceDescription("Sample data root directory", sourceName));

            addPathConfig(configBuilder, LINX_DIR_CFG, LINX_DIR_DESC, sourceName);
            addPathConfig(configBuilder, LINX_GERMLINE_DIR_CFG, LINX_GERMLINE_DIR_DESC, sourceName);
            addPathConfig(configBuilder, PURPLE_DIR_CFG, PURPLE_DIR_DESC, sourceName);
            addPathConfig(configBuilder, LILAC_DIR_CFG, LILAC_DIR_DESC, sourceName);
            addPathConfig(configBuilder, CHORD_DIR_CFG, CHORD_DIR_DESC, sourceName);
            addPathConfig(configBuilder, CUPPA_DIR_CFG, CUPPA_DIR_DESC, sourceName);
            addPathConfig(configBuilder, PEACH_DIR_CFG, PEACH_DIR_DESC, sourceName);
            addPathConfig(configBuilder, VIRUS_DIR_CFG, VIRUS_DIR_DESC, sourceName);
            addPathConfig(configBuilder, CIDER_DIR_CFG, CIDER_DIR_DESC, sourceName);
            addPathConfig(configBuilder, TEAL_DIR_CFG, TEAL_DIR_DESC, sourceName);
            addPathConfig(configBuilder, TUMOR_FLAGSTAT, formSourceDescription("Tumor flagstat", sourceName), sourceName);
            addPathConfig(configBuilder, GERMLINE_FLAGSTAT, formSourceDescription("Germline flagstat", sourceName), sourceName);
            addPathConfig(configBuilder, TUMOR_BAM_METRICS, formSourceDescription("Tumor BAM metrics", sourceName), sourceName);
            addPathConfig(configBuilder, GERMLINE_BAM_METRICS, formSourceDescription("Germline BAM metrics", sourceName), sourceName);
            addPathConfig(configBuilder, SNP_GENOTYPE, formSourceDescription("SNP genotype", sourceName), sourceName);

            configBuilder.addPath(
                    formSourceConfig(SOMATIC_VCF, sourceName), false,
                    formSourceDescription("Non-Purple somatic VCF (eg Pave, Sage)", sourceName));

            configBuilder.addPath(
                    formSourceConfig(SOMATIC_UNFILTERED_VCF, sourceName), false,
                    formSourceDescription("VCF to search for filtered variants", sourceName));

            configBuilder.addConfigItem(formSourceConfig(PIPELINE_VERSION, sourceName), false, formSourceDescription(
                    "Pipeline version for deriving tool directories from sample directory in non-default format. '" + PIPELINE_OUTPUT_STRUCTURE
                            + "' needs to also be set for this to work. Options: " + optionList(PIPELINE_VERSION_VALUE_TRANSLATIONS.keySet()), sourceName));
            configBuilder.addConfigItem(formSourceConfig(PIPELINE_OUTPUT_STRUCTURE, sourceName), false, formSourceDescription(
                    "Pipeline type for deriving tool directories from sample directory in non-default format. '" + PIPELINE_VERSION
                            + "' needs to also be set for this to work. Options: PIPELINE5, ONCOANALYSER, DATABASE", sourceName));
        }
    }

    private static String optionList(final Set<String> options)
    {
        return options.stream().sorted().collect(Collectors.joining(", "));
    }

    private static void addPathConfig(final ConfigBuilder configBuilder, final String toolDir, final String toolDesc,
            final String sourceName)
    {
        configBuilder.addPrefixedPath(
                formSourceConfig(toolDir, sourceName), false, formSourceDescription(toolDesc, sourceName),
                formSourceConfig(SAMPLE_DIR, sourceName));
    }

    private static String formSourceDescription(final String desc, final String sourceName)
    {
        return format("%s: source %s", desc, sourceName);
    }

    private static String formSourceConfig(final String config, final String sourceName)
    {
        return format("%s_%s", config, sourceName);
    }

    private static String getConfigValue(final ConfigBuilder configBuilder, final String config, final String sourceName)
    {
        return configBuilder.getValue(formSourceConfig(config, sourceName), "");
    }

    public static FileSources fromConfig(final String sourceName, final ConfigBuilder configBuilder)
    {
        return new FileSources(sourceName,
                checkAddDirSeparator(getConfigValue(configBuilder, SAMPLE_DIR, sourceName)),
                getConfigValue(configBuilder, LINX_DIR_CFG, sourceName),
                getConfigValue(configBuilder, PURPLE_DIR_CFG, sourceName),
                getConfigValue(configBuilder, LINX_GERMLINE_DIR_CFG, sourceName),
                getConfigValue(configBuilder, CUPPA_DIR_CFG, sourceName),
                getConfigValue(configBuilder, LILAC_DIR_CFG, sourceName),
                getConfigValue(configBuilder, CHORD_DIR_CFG, sourceName),
                getConfigValue(configBuilder, PEACH_DIR_CFG, sourceName),
                getConfigValue(configBuilder, VIRUS_DIR_CFG, sourceName),
                getConfigValue(configBuilder, SOMATIC_VCF, sourceName),
                getConfigValue(configBuilder, SOMATIC_UNFILTERED_VCF, sourceName),
                getConfigValue(configBuilder, TUMOR_FLAGSTAT, sourceName),
                getConfigValue(configBuilder, GERMLINE_FLAGSTAT, sourceName),
                getConfigValue(configBuilder, TUMOR_BAM_METRICS, sourceName),
                getConfigValue(configBuilder, GERMLINE_BAM_METRICS, sourceName),
                getConfigValue(configBuilder, SNP_GENOTYPE, sourceName),
                getConfigValue(configBuilder, CIDER_DIR_CFG, sourceName),
                getConfigValue(configBuilder, TEAL_DIR_CFG, sourceName),
                resolvePipelineVersion(getConfigValue(configBuilder, PIPELINE_VERSION, sourceName)),
                resolvePipelineOutputStructure(getConfigValue(configBuilder, PIPELINE_OUTPUT_STRUCTURE, sourceName))
        );
    }

    private static PipelineVersion resolvePipelineVersion(String pipelineVersion)
    {
        if(pipelineVersion.isEmpty())
        {
            return null;
        }
        else if(PIPELINE_VERSION_VALUE_TRANSLATIONS.containsKey(pipelineVersion))
        {
            return PipelineVersion.fromString(PIPELINE_VERSION_VALUE_TRANSLATIONS.get(pipelineVersion));
        }
        else
        {
            throw new IllegalArgumentException("Pipeline version " + pipelineVersion + " is not one of the allowed options: "
                    + optionList(PIPELINE_VERSION_VALUE_TRANSLATIONS.keySet()));
        }
    }

    private static PipelineOutputStructure resolvePipelineOutputStructure(String pipelineOutputStructure)
    {
        if(pipelineOutputStructure.isEmpty())
            return null;

        try
        {
            return PipelineOutputStructure.valueOf(pipelineOutputStructure);
        }
        catch(IllegalArgumentException e)
        {
            throw new IllegalArgumentException("Invalid pipeline output structure: " + pipelineOutputStructure);
        }
    }
}
