package com.hartwig.hmftools.compar.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
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
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TEAL_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TEAL_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.ComparConfig.NEW_SOURCE;
import static com.hartwig.hmftools.compar.ComparConfig.REF_SOURCE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.pipeline.PipelineToolDirectories;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class FileSources
{
    public final String Source;
    public final String Linx;
    public final String LinxGermline;
    public final String Purple;
    public final String Cuppa;
    public final String Lilac;
    public final String Chord;
    public final String SomaticVcf;
    public final String SomaticUnfilteredVcf;

    public final String Cider;

    public final String Teal;

    private static final String SAMPLE_DIR = "sample_dir";
    private static final String SOMATIC_VCF = "somatic_vcf";
    private static final String SOMATIC_UNFILTERED_VCF = "somatic_unfiltered_vcf";

    public FileSources(final String source, final String linx, final String purple, final String linxGermline, final String cuppa,
            final String lilac, final String chord, final String somaticVcf, final String somaticUnfilteredVcf,
            final String cider, String teal)
    {
        Source = source;
        Linx = linx;
        LinxGermline = linxGermline;
        Purple = purple;
        Cuppa = cuppa;
        Lilac = lilac;
        Chord = chord;
        SomaticVcf = somaticVcf;
        SomaticUnfilteredVcf = somaticUnfilteredVcf;
        Cider = cider;
        Teal = teal;
    }

    public static FileSources sampleInstance(final FileSources fileSources, final String sampleId)
    {
        return new FileSources(
                fileSources.Source,
                convertWildcardSamplePath(fileSources.Linx, sampleId),
                convertWildcardSamplePath(fileSources.Purple, sampleId),
                convertWildcardSamplePath(fileSources.LinxGermline, sampleId),
                convertWildcardSamplePath(fileSources.Cuppa, sampleId),
                convertWildcardSamplePath(fileSources.Lilac, sampleId),
                convertWildcardSamplePath(fileSources.Chord, sampleId),
                convertWildcardSamplePath(fileSources.SomaticVcf, sampleId),
                convertWildcardSamplePath(fileSources.SomaticUnfilteredVcf, sampleId),
                convertWildcardSamplePath(fileSources.Cider, sampleId),
                convertWildcardSamplePath(fileSources.Teal, sampleId)
                );
    }

    public static RefGenomeVersion liftoverSourceGenomeVersion(final String source)
    {
        return source.equals(REF_SOURCE) ? V37 : V38;
    }

    private static void addPathConfig(final ConfigBuilder configBuilder, final String toolDir, final String toolDesc, final String sourceName)
    {
        configBuilder.addPrefixedPath(
                formSourceConfig(toolDir, sourceName), false, formSourceDescription(toolDesc, sourceName),
                formSourceConfig(SAMPLE_DIR, sourceName));
    }

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
            addPathConfig(configBuilder, CIDER_DIR_CFG, CIDER_DIR_DESC, sourceName);
            addPathConfig(configBuilder, TEAL_DIR_CFG, TEAL_DIR_DESC, sourceName);

            configBuilder.addPath(
                    formSourceConfig(SOMATIC_VCF, sourceName), false,
                    formSourceDescription("Non-Purple somatic VCF (eg Pave, Sage)", sourceName));

            configBuilder.addPath(
                    formSourceConfig(SOMATIC_UNFILTERED_VCF, sourceName), false,
                    formSourceDescription("VCF to search for filtered variants", sourceName));
        }
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
        String sampleDir = checkAddDirSeparator(getConfigValue(configBuilder, SAMPLE_DIR, sourceName));

        String linxDir = getDirectory(configBuilder, sampleDir, PipelineToolDirectories.LINX_SOMATIC_DIR, LINX_DIR_CFG, sourceName);

        String linxGermlineDir = getDirectory(
                configBuilder, sampleDir, PipelineToolDirectories.LINX_GERMLINE_DIR, LINX_GERMLINE_DIR_CFG, sourceName);

        String purpleDir = getDirectory(configBuilder, sampleDir, PipelineToolDirectories.PURPLE_DIR, PURPLE_DIR_CFG, sourceName);
        String cuppaDir = getDirectory(configBuilder, sampleDir, PipelineToolDirectories.CUPPA_DIR, CUPPA_DIR_CFG, sourceName);
        String lilacDir = getDirectory(configBuilder, sampleDir, PipelineToolDirectories.LILAC_DIR, LILAC_DIR_CFG, sourceName);
        String chordDir = getDirectory(configBuilder, sampleDir, PipelineToolDirectories.CHORD_DIR, CHORD_DIR_CFG, sourceName);

        String somaticVcf = getConfigValue(configBuilder, SOMATIC_VCF, sourceName);
        String somaticUnfilteredVcf = getConfigValue(configBuilder, SOMATIC_UNFILTERED_VCF, sourceName);

        String ciderDir = getDirectory(configBuilder, sampleDir, PipelineToolDirectories.CIDER_DIR, CIDER_DIR_CFG, sourceName);
        String tealDir = getDirectory(configBuilder, sampleDir, PipelineToolDirectories.TEAL_DIR, TEAL_DIR_CFG, sourceName);

        return new FileSources(
                sourceName, linxDir, purpleDir, linxGermlineDir, cuppaDir, lilacDir, chordDir, somaticVcf, somaticUnfilteredVcf,
                ciderDir, tealDir);
    }

    private static String getDirectory(
            final ConfigBuilder configBuilder, final String sampleDir, final String toolDefaultDir,
            final String config, final String sourceName)
    {
        // if a tool directory is specified in config, then it overrides the default pipeline directory
        // if the root sample directory is specified, then the tool directory is relative to that, otherwise is absolute

        String configStr = formSourceConfig(config, sourceName);

        if(!configBuilder.hasValue(configStr) && sampleDir.isEmpty())
            return "";

        String toolDir = configBuilder.getValue(configStr, toolDefaultDir);

        String directory = "";

        if(sampleDir.isEmpty())
            directory = toolDir;
        else if(toolDir.isEmpty())
            directory = sampleDir;
        else
            directory = format("%s%s", sampleDir, toolDir);

        return checkAddDirSeparator(directory);
    }
}
