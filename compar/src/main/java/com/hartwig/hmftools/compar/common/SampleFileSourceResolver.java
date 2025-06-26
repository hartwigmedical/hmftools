package com.hartwig.hmftools.compar.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.pipeline.PipelineToolDirectories;
import com.hartwig.pipeline.reference.api.DataType;
import com.hartwig.pipeline.reference.api.Pipeline;
import com.hartwig.pipeline.reference.api.PipelineFile;
import com.hartwig.pipeline.reference.api.PipelineFilePath;
import com.hartwig.pipeline.reference.api.PipelineFiles;
import com.hartwig.pipeline.reference.api.PipelineMolecule;
import com.hartwig.pipeline.reference.api.PipelineRun;
import com.hartwig.pipeline.reference.api.SampleType;
import com.hartwig.pipeline.reference.api.Tool;

public record SampleFileSourceResolver(FileSources fileSources, String sampleId, String germlineSampleId)
{
    public String resolveLinxDirectory()
    {
        return resolveDirectory(fileSources.linx(), sampleId, germlineSampleId, null,
                DataType.LINX_DRIVER_CATALOG, SampleType.TUMOR, PipelineToolDirectories.LINX_SOMATIC_DIR);
    }

    public String resolveLinxGermlineDirectory()
    {
        return resolveDirectory(fileSources.linxGermline(), sampleId, germlineSampleId, null,
                DataType.LINX_GERMLINE_DRIVER_CATALOG, SampleType.REFERENCE, PipelineToolDirectories.LINX_GERMLINE_DIR);
    }

    public String resolvePurpleDirectory()
    {
        return resolveDirectory(fileSources.purple(), sampleId, germlineSampleId, null,
                DataType.PURPLE_SOMATIC_DRIVER_CATALOG, null, PipelineToolDirectories.PURPLE_DIR);
    }

    public String resolveCuppaDirectory()
    {
        return resolveDirectory(fileSources.cuppa(), sampleId, germlineSampleId, null,
                DataType.CUPPA_PRED_SUMM, null, PipelineToolDirectories.CUPPA_DIR);
    }

    public String resolveLilacDirectory()
    {
        return resolveDirectory(fileSources.lilac(), sampleId, germlineSampleId, null,
                DataType.LILAC_OUTPUT, null, PipelineToolDirectories.LILAC_DIR);
    }

    public String resolveChordDirectory()
    {
        return resolveDirectory(fileSources.chord(), sampleId, germlineSampleId, null,
                DataType.CHORD_PREDICTION, null, PipelineToolDirectories.CHORD_DIR);
    }

    public String resolvePeachDirectory()
    {
        return resolveDirectory(fileSources.peach(), sampleId, germlineSampleId, null,
                DataType.PEACH_GENOTYPE, null, PipelineToolDirectories.PEACH_DIR);
    }

    public String resolveVirusDirectory()
    {
        return resolveDirectory(fileSources.virus(), sampleId, germlineSampleId, null,
                DataType.VIRUS_INTERPRETATION, null, PipelineToolDirectories.VIRUS_INTERPRETER_DIR);
    }

    public String resolveTumorFlagstatDirectory()
    {
        return resolveDirectory(fileSources.tumorFlagstat(), sampleId, germlineSampleId, null,
                DataType.METRICS_FLAG_COUNT, SampleType.TUMOR, "*/bam_metrics");
    }

    public String resolveGermlineFlagstatDirectory()
    {
        return resolveDirectory(fileSources.germlineFlagstat(), sampleId, germlineSampleId, null,
                DataType.METRICS_FLAG_COUNT, SampleType.REFERENCE, "$/bam_metrics");
    }

    public String resolveTumorBamMetricsDirectory()
    {
        return resolveDirectory(fileSources.tumorBamMetrics(), sampleId, germlineSampleId, null,
                DataType.METRICS_COVERAGE, SampleType.TUMOR, "*/bam_metrics");
    }

    public String resolveGermlineBamMetricsDirectory()
    {
        return resolveDirectory(fileSources.germlineBamMetrics(), sampleId, germlineSampleId, null,
                DataType.METRICS_COVERAGE, SampleType.REFERENCE, "$/bam_metrics");
    }

    public String resolveSnpGenotypeDirectory()
    {
        return resolveDirectory(fileSources.snpGenotype(), sampleId, germlineSampleId, Tool.SNP_GENOTYPE,
                null, null, "$/snp_genotype");
    }

    public String resolveCiderDirectory()
    {
        return resolveDirectory(fileSources.cider(), sampleId, germlineSampleId, null,
                DataType.CIDER_VDJ, null, PipelineToolDirectories.CIDER_DIR);
    }

    public String resolveTealDirectory()
    {
        return resolveDirectory(fileSources.teal(), sampleId, germlineSampleId, null,
                DataType.TEAL_SOMATIC_BREAKEND, null, PipelineToolDirectories.TEAL_DIR);
    }

    public String resolveSomaticVcfPath()
    {
        return convertWildcardSamplePath(fileSources.somaticVcf(), sampleId, germlineSampleId);
    }

    public String resolveSomaticUnfilteredVcfPath()
    {
        return resolvePath(fileSources.somaticUnfilteredVcf(), sampleId, germlineSampleId, DataType.SOMATIC_VARIANTS_SAGE);
    }

    private String resolveDirectory(final String toolDirectoryFromConfig, final String sampleId,
            final String germlineSampleId, final Tool tool, final DataType dataType, final SampleType sampleType,
            final String hardCodedDefaultToolDirectory)
    {
        // if a tool directory is specified in config, then it overrides the default pipeline directory
        // if the root sample directory is specified, then the tool directory is relative to that, otherwise is absolute
        if(toolDirectoryFromConfig.isEmpty() && fileSources.sampleDir().isEmpty())
        {
            return "";
        }

        String directory;
        if(fileSources.sampleDir().isEmpty())
        {
            directory = toolDirectoryFromConfig;
        }
        else if(!toolDirectoryFromConfig.isEmpty())
        {
            directory = prependSampleDirectory(toolDirectoryFromConfig);
        }
        else
        {
            String defaultToolDirectory =
                    resolveDefaultToolDirectory(sampleId, germlineSampleId, tool, dataType, sampleType, hardCodedDefaultToolDirectory);
            directory = prependSampleDirectory(defaultToolDirectory);
        }

        return checkAddDirSeparator(convertWildcardSamplePath(directory, sampleId, germlineSampleId));
    }

    private String prependSampleDirectory(final String relativePath)
    {
        return format("%s%s", checkAddDirSeparator(fileSources.sampleDir()), relativePath);
    }

    private String resolveDefaultToolDirectory(final String sampleId, final String germlineSampleId, final Tool tool,
            final DataType dataType, final SampleType sampleType, final String hardCodedDefaultToolDirectory)
    {
        if(fileSources.pipelineVersion() != null && fileSources.pipelineOutputStructure() != null)
        {
            // try to use pipeline-reference-data library to find default tool directory for specific pipeline and version
            PipelineFilePath pathOrNull = resolveImpliedPipelinePathOrNull(sampleId, germlineSampleId, tool, dataType, sampleType);
            if(pathOrNull != null)
            {
                return pathOrNull.getFull().getParent().toString();
            }
        }

        return hardCodedDefaultToolDirectory;

    }

    private String resolvePath(String pathFromConfig, String sampleId, String germlineSampleId, DataType dataType)
    {
        if(!pathFromConfig.isEmpty())
        {
            return convertWildcardSamplePath(pathFromConfig, sampleId, germlineSampleId);
        }

        if(fileSources.pipelineVersion() != null && fileSources.pipelineOutputStructure() != null)
        {
            // try to use pipeline-reference-data library to find expected path to file
            PipelineFilePath pathOrNull = resolveImpliedPipelinePathOrNull(sampleId, germlineSampleId, null, dataType, null);
            if(pathOrNull != null)
            {
                return convertWildcardSamplePath(prependSampleDirectory(pathOrNull.getFull().toString()), sampleId, germlineSampleId);
            }
        }

        return "";
    }

    private PipelineFilePath resolveImpliedPipelinePathOrNull(final String sampleId, final String germlineSampleId, final Tool tool,
            final DataType dataType, final SampleType sampleType)
    {
        Pipeline pipeline = Pipeline.fromMoleculeAndVersion(PipelineMolecule.DNA, fileSources.pipelineVersion());
        Predicate<PipelineFile>[] filters = createPipelineFileFilters(tool, dataType, sampleType);

        List<PipelineFile> matchingPipelineFiles = PipelineFiles.get(new PipelineRun(pipeline, sampleId, germlineSampleId), filters);

        if(!matchingPipelineFiles.isEmpty())
        {
            return matchingPipelineFiles.get(0).getPathOrNull(fileSources.pipelineOutputStructure());
        }
        else
        {
            return null;
        }
    }

    private static Predicate<PipelineFile>[] createPipelineFileFilters(final Tool tool, final DataType dataType,
            final SampleType sampleType)
    {
        ArrayList<Predicate<PipelineFile>> filters = new ArrayList<>();
        if(tool != null)
            filters.add(PipelineFiles.toolIsAnyOf(tool));
        if(dataType != null)
            filters.add(PipelineFiles.dataTypeIsAnyOf(dataType));
        if(sampleType != null)
            filters.add(PipelineFiles.sampleTypeIs(sampleType));
        return filters.toArray(new Predicate[0]);
    }
}
