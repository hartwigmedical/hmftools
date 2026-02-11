package com.hartwig.hmftools.common.pipeline;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;

import java.io.IOException;
import java.util.Arrays;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public record PipelineToolDirectories(
        String amberDir,
        String chordDir,
        String ciderDir,
        String cobaltDir,
        String cuppaDir,
        String esveeDir,
        String germlineFlagstatDir,
        String germlineMetricsDir,
        String isofoxDir,
        String lilacDir,
        String linxGermlineDir,
        String linxSomaticDir,
        String orangeDir,
        String paveGermlineDir,
        String paveSomaticDir,
        String peachDir,
        String purpleDir,
        String sageGermlineDir,
        String sageSomaticDir,
        String sigsDir,
        String snpGenotypeDir,
        String tealDir,
        String tumorFlagstatDir,
        String tumorMetricsDir,
        String vChordDir,
        String virusBreakendDir,
        String virusInterpreterDir
) {
    public static final String PIPELINE_FORMAT_CFG = "pipeline_format";
    public static final String PIPELINE_FORMAT_DESC =
            "Assumed directory structure for tool directories. Possible values: " + Arrays.stream(PipelineOutputStructure.values())
                    .map(Enum::name)
                    .collect(Collectors.joining(", ")) + ". Default: " + PipelineOutputStructure.PIP5_V6_0.name();
    public static final String PIPELINE_FORMAT_FILE_CFG = "pipeline_format_file";
    public static final String PIPELINE_FORMAT_FILE_DESC = "File describing expected tool directory structure.";

    public static final PipelineToolDirectories OA_V2_0_FORMAT = new PipelineToolDirectories(
            "amber",
            "chord",
            "cider",
            "cobalt",
            "cuppa",
            "esvee",
            "bamtools/$_bamtools",
            "bamtools/$_bamtools",
            "isofox",
            "lilac",
            "linx/germline_annotations",
            "linx/somatic_annotations",
            "orange",
            "pave",
            "pave",
            "peach",
            "purple",
            "sage/germline",
            "sage/somatic",
            "sigs",
            "",
            "teal",
            "bamtools/*_bamtools",
            "bamtools/*_bamtools",
            "vchord",  // not yet implemented in this version
            "virusbreakend",
            "virusinterpreter"
    );

    public static final PipelineToolDirectories OA_V2_2_FORMAT = new PipelineToolDirectories(
            "amber",
            "chord",
            "cider",
            "cobalt",
            "cuppa",
            "esvee",
            "bamtools/$_bamtools",
            "bamtools/$_bamtools",
            "isofox",
            "lilac",
            "linx/germline_annotations",
            "linx/somatic_annotations",
            "orange",
            "pave",
            "pave",
            "peach",
            "purple",
            "sage/germline",
            "sage/somatic",
            "sigs",
            "",
            "teal",
            "bamtools/*_bamtools",
            "bamtools/*_bamtools",
            "vchord",  // not yet implemented in this version
            "virusbreakend",
            "virusinterpreter"
    );

    public static final PipelineToolDirectories OA_V2_3_FORMAT = OA_V2_2_FORMAT;
    public static final PipelineToolDirectories OA_V3_0_FORMAT = OA_V2_2_FORMAT;

    public static final PipelineToolDirectories PIP5_V6_0_FORMAT = new PipelineToolDirectories(
            "amber",
            "chord",
            "cider",
            "cobalt",
            "cuppa",
            "esvee",
            "$/bam_metrics",
            "$/bam_metrics",
            "isofox",
            "lilac",
            "linx_germline",
            "linx",
            "orange",
            "pave_germline",
            "pave_somatic",
            "peach",
            "purple",
            "sage_germline",
            "sage_somatic",
            "sigs",
            "$/snp_genotype",
            "teal",
            "*/bam_metrics",
            "*/bam_metrics",
            "vchord",
            "virusbreakend",
            "virusintrprtr"
    );

    public static final PipelineToolDirectories DB_V6_0_FORMAT = new PipelineToolDirectories(
            "amber",
            "chord",
            "cider",
            "cobalt",
            "cuppa",
            "esvee",
            "bamtools/$_bamtools",
            "bamtools/$_bamtools",
            "isofox",
            "lilac",
            "linx/germline_annotations",
            "linx/somatic_annotations",
            "orange",
            "pave/germline",
            "pave/somatic",
            "peach",
            "purple",
            "sage/germline",
            "sage/somatic",
            "sigs",
            "snp_genotype/$",
            "teal",
            "bamtools/*_bamtools",
            "bamtools/*_bamtools",
            "vchord",  // not yet implemented in this version
            "virusbreakend",
            "virusinterpreter"
    );

    public static PipelineToolDirectories resolveToolDirectories(
            final ConfigBuilder configBuilder, final String pipelineFormatConfigStr,
            final String pipelineFormatFileConfigStr, final String tumorSampleId)
    {
        return resolveToolDirectories(configBuilder, pipelineFormatConfigStr, pipelineFormatFileConfigStr, tumorSampleId, null);
    }

    public static PipelineToolDirectories resolveToolDirectories(
            final ConfigBuilder configBuilder, final String pipelineFormatConfigStr,
            final String pipelineFormatFileConfigStr, final String tumorSampleId, final String normalSampleId)
    {
        PipelineToolDirectories withWildcardSampleIds =
                resolveToolDirectories(configBuilder, pipelineFormatConfigStr, pipelineFormatFileConfigStr);
        return withWildcardSampleIds.resolveSampleIds(tumorSampleId, normalSampleId);
    }

    public static PipelineToolDirectories resolveToolDirectories(
            final ConfigBuilder configBuilder, final String pipelineFormatConfigStr,
            final String pipelineFormatFileConfigStr)
    {
        if(configBuilder.hasValue(pipelineFormatFileConfigStr))
        {
            return resolveToolDirectoriesFromFile(configBuilder.getValue(pipelineFormatFileConfigStr));
        }
        else if(configBuilder.hasValue(pipelineFormatConfigStr))
        {
            PipelineOutputStructure outputStructure = PipelineOutputStructure.valueOf(configBuilder.getValue(pipelineFormatConfigStr));
            return PipelineToolDirectories.resolveToolDirectoriesFromDefault(outputStructure);
        }
        else
        {
            return PipelineToolDirectories.resolveToolDirectoriesFromDefault(PipelineOutputStructure.OA_V2_2);
        }
    }

    public static void addPipelineFormatOptions(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(PIPELINE_FORMAT_CFG, false, PIPELINE_FORMAT_DESC);
        configBuilder.addPath(PIPELINE_FORMAT_FILE_CFG, false, PIPELINE_FORMAT_FILE_DESC);
    }

    private static PipelineToolDirectories resolveToolDirectoriesFromDefault(final PipelineOutputStructure outputStructure)
    {
        return switch(outputStructure)
        {
            case OA_V2_0 -> OA_V2_0_FORMAT;
            case OA_V2_2 -> OA_V2_2_FORMAT;
            case OA_V2_3 -> OA_V2_3_FORMAT;
            case OA_V3_0 -> OA_V3_0_FORMAT;
            case PIP5_V6_0 -> PIP5_V6_0_FORMAT;
            case DB_V6_0 -> DB_V6_0_FORMAT;
        };
    }

    private static PipelineToolDirectories resolveToolDirectoriesFromFile(final String filePath)
    {
        try
        {
            return PipelineToolDirectoriesFile.read(filePath);
        }
        catch(IOException e)
        {
            throw new IllegalArgumentException("Could not load tool subdirectory config file: " + filePath, e);
        }
    }

    private PipelineToolDirectories resolveSampleIds(final String tumorSampleId, final String normalSampleId)
    {
        return new PipelineToolDirectories(
                convertWildcardSamplePath(amberDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(chordDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(ciderDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(cobaltDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(cuppaDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(esveeDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(germlineFlagstatDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(germlineMetricsDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(isofoxDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(lilacDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(linxGermlineDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(linxSomaticDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(orangeDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(paveGermlineDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(paveSomaticDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(peachDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(purpleDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(sageGermlineDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(sageSomaticDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(sigsDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(snpGenotypeDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(tealDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(tumorFlagstatDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(tumorMetricsDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(vChordDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(virusBreakendDir, tumorSampleId, normalSampleId),
                convertWildcardSamplePath(virusInterpreterDir, tumorSampleId, normalSampleId)
        );
    }
}
