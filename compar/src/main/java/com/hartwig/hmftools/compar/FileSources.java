package com.hartwig.hmftools.compar;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_DESC;
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
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.CommonUtils.SUB_ITEM_DELIM;
import static com.hartwig.hmftools.compar.ComparConfig.NEW_SOURCE;
import static com.hartwig.hmftools.compar.ComparConfig.REF_SOURCE;

import java.util.List;

import com.google.common.collect.Lists;
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

    public final boolean RequiresLiftover;

    private static final String SAMPLE_DIR = "sample_dir";
    private static final String SOMATIC_VCF = "somatic_vcf";
    private static final String SOMATIC_UNFILTERED_VCF = "somatic_unfiltered_vcf";
    private static final String REQUIRES_LIFTOVER = "liftover";

    public FileSources(final String source, final String linx, final String purple, final String linxGermline, final String cuppa,
            final String lilac, final String chord, final String somaticVcf, final String somaticUnfilteredVcf, boolean requiresLiftover)
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
        RequiresLiftover = requiresLiftover;
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
                fileSources.RequiresLiftover);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        List<String> sourceNames = Lists.newArrayList(REF_SOURCE, NEW_SOURCE);

        for(String sourceName : sourceNames)
        {
            configBuilder.addPath(
                    formSourceConfig(SAMPLE_DIR, sourceName), false,
                    formSourceDescription("Sample data root directory", sourceName));

            configBuilder.addPrefixedPath(
                    formSourceConfig(LINX_DIR_CFG, sourceName), false, formSourceDescription(LINX_DIR_DESC, sourceName), SAMPLE_DIR);

            configBuilder.addPath(
                    formSourceConfig(LINX_GERMLINE_DIR_CFG, sourceName), false, formSourceDescription(LINX_GERMLINE_DIR_DESC, sourceName));

            configBuilder.addPath(
                    formSourceConfig(PURPLE_DIR_CFG, sourceName), false, formSourceDescription(PURPLE_DIR_DESC, sourceName));

            configBuilder.addPath(
                    formSourceConfig(LILAC_DIR_CFG, sourceName), false, formSourceDescription(LILAC_DIR_DESC, sourceName));

            configBuilder.addPath(
                    formSourceConfig(CHORD_DIR_CFG, sourceName), false, formSourceDescription(CHORD_DIR_DESC, sourceName));

            configBuilder.addPath(
                    formSourceConfig(CUPPA_DIR_CFG, sourceName), false, formSourceDescription(CUPPA_DIR_DESC, sourceName));

            configBuilder.addPath(
                    formSourceConfig(SOMATIC_VCF, sourceName), false,
                    formSourceDescription("Non-Purple somatic VCF (eg Pave, Sage)", sourceName));

            configBuilder.addPath(
                    formSourceConfig(SOMATIC_UNFILTERED_VCF, sourceName), false,
                    formSourceDescription("VCF to search for filtered variants", sourceName));

            configBuilder.addFlag(
                    formSourceConfig(REQUIRES_LIFTOVER, sourceName), formSourceDescription("Liftover positions", sourceName));
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
        boolean requiresLiftover = configBuilder.hasFlag(formSourceConfig(REQUIRES_LIFTOVER, sourceName));

        return new FileSources(
                sourceName, linxDir, purpleDir, linxGermlineDir, cuppaDir, lilacDir, chordDir,
                somaticVcf, somaticUnfilteredVcf, requiresLiftover);
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

    @Deprecated
    private static String getDirectory(final String sampleDir, final String typeDir)
    {
        String directory = "";

        if(sampleDir.isEmpty())
            directory = typeDir;
        else if(typeDir.isEmpty())
            directory = sampleDir;
        else
            directory = format("%s%s", sampleDir, typeDir);

        return checkAddDirSeparator(directory);
    }

    @Deprecated
    public static FileSources fromConfig(final String sourceName, final String fileSourceStr)
    {
        String[] values = fileSourceStr.split(ITEM_DELIM, -1);

        String sampleDir = "";

        int itemIndex = 0;

        if(values[itemIndex].startsWith(SAMPLE_DIR))
        {
            sampleDir = checkAddDirSeparator(values[itemIndex].split(SUB_ITEM_DELIM)[1]);
            ++itemIndex;
        }

        // by default use the sample root directory and then default pipeline directory names per tool and mode
        String linxDir = getDirectory(sampleDir, PipelineToolDirectories.LINX_SOMATIC_DIR);
        String linxGermlineDir = getDirectory(sampleDir, PipelineToolDirectories.LINX_GERMLINE_DIR);
        String purpleDir = getDirectory(sampleDir, PipelineToolDirectories.PURPLE_DIR);
        String cuppaDir = getDirectory(sampleDir, PipelineToolDirectories.CUPPA_DIR);
        String lilacDir = getDirectory(sampleDir, PipelineToolDirectories.LILAC_DIR);
        String chordDir = getDirectory(sampleDir, PipelineToolDirectories.CHORD_DIR);
        String somaticVcf = "";
        String somaticUnfilteredVcf = "";

        boolean requiresLiftover = false;

        for(int i = itemIndex; i < values.length; ++i)
        {
            if(values[i].equals(REQUIRES_LIFTOVER))
            {
                requiresLiftover = true;
                continue;
            }

            String[] itemStr = values[i].split(SUB_ITEM_DELIM);

            if(itemStr.length != 2)
                return null;

            String type = itemStr[0];
            String value = itemStr[1];

            if(type.equals(LINX_DIR_CFG))
            {
                linxDir = getDirectory(sampleDir, value);
            }
            else if(type.equals(LINX_GERMLINE_DIR_CFG))
            {
                linxGermlineDir = getDirectory(sampleDir, value);
            }
            else if(type.equals(PURPLE_DIR_CFG))
            {
                purpleDir = getDirectory(sampleDir, value);
            }
            else if(type.equals(LILAC_DIR_CFG))
            {
                lilacDir = value;
            }
            else if(type.equals(CHORD_DIR_CFG))
            {
                chordDir = value;
            }
            else if(type.equals(CUPPA_DIR_CFG))
            {
                cuppaDir = value;
            }
            else if(type.equals(SOMATIC_VCF))
            {
                somaticVcf = value;
            }
            else if(type.equals(SOMATIC_UNFILTERED_VCF))
            {
                somaticUnfilteredVcf = value;
            }
        }

        return new FileSources(
                sourceName, linxDir, purpleDir, linxGermlineDir, cuppaDir, lilacDir, chordDir,
                somaticVcf, somaticUnfilteredVcf, requiresLiftover);
    }
}
