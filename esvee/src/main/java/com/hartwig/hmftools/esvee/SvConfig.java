package com.hartwig.hmftools.esvee;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.CONFIG_FILE_DELIM;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.esvee.SvConstants.ASSEMBLY_BAM_FILE_ID;
import static com.hartwig.hmftools.esvee.SvConstants.DEFAULT_HTML_SUMMARY_DIR;
import static com.hartwig.hmftools.esvee.SvConstants.SV_PREP_JUNCTIONS_FILE_ID;
import static com.hartwig.hmftools.esvee.WriteType.BREAKEND_TSV;
import static com.hartwig.hmftools.esvee.WriteType.HTML_SUMMARY;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.ValidationStringency;

public class SvConfig
{
    public final List<String> SampleNames;
    public final List<String> BamFiles;
    public final List<String> JunctionFiles;

    public final RefGenomeVersion RefGenVersion;
    public final RefGenomeCoordinates RefGenomeCoords;
    public final String RefGenomeFile;
    public final RefGenomeInterface RefGenome;
    public final String RefGenomeImageFile;

    public final ValidationStringency BamStringency;

    public final String VcfFile;
    public final List<WriteType> WriteTypes;

    public final String OutputDir;
    public final String HtmlOutputDir;
    public final boolean PlotDiagrams;

    public final boolean PerfDebug;
    public final double PerfLogTime;
    public final boolean OtherDebug;

    public final int Threads;

    public static final String OUTPUT_VCF = "output_vcf";
    public static final String REF_GENOME_IMAGE = "ref_genome_image";
    public static final String JUNCTION_FILES = "junction_files";

    // alternative to specifically load a tumor and/or ref sample and BAM
    public static final String SAMPLE_IDS = "samples";
    public static final String SAMPLE_BAMS = "bam_files";

    public static final String WRITE_TYPES = "write_types";
    public static final String HTML_SUMMARY_DIR = "html_dir";
    public static final String PLOT_DIAGRAMS = "plot_diagrams";
    public static final String OTHER_DEBUG = "other_debug";
    public static final String PERF_LOG_TIME = "perf_log_time";

    public static final Logger SV_LOGGER = LogManager.getLogger(SvConfig.class);

    public SvConfig(final ConfigBuilder configBuilder)
    {
        SampleNames = Arrays.stream(configBuilder.getValue(SAMPLE_IDS).split(CONFIG_FILE_DELIM)).collect(Collectors.toList());
        BamFiles = Arrays.stream(configBuilder.getValue(SAMPLE_BAMS).split(CONFIG_FILE_DELIM)).collect(Collectors.toList());

        if(SampleNames.isEmpty() || SampleNames.size() != BamFiles.size())
        {
            SV_LOGGER.error("sample IDs must match bam files");
            System.exit(1);
        }

        JunctionFiles = Lists.newArrayList();

        if(configBuilder.hasValue(JUNCTION_FILES))
        {
            Arrays.stream(configBuilder.getValue(JUNCTION_FILES).split(CONFIG_FILE_DELIM)).forEach(x -> JunctionFiles.add(x));
        }
        else
        {
            String bamPath = pathFromFile(BamFiles.get(0));

            for(String sampleId : SampleNames)
            {
                String junctionFile = bamPath + sampleId + SV_PREP_JUNCTIONS_FILE_ID;

                if(Files.exists(Paths.get(junctionFile)))
                    JunctionFiles.add(junctionFile);
            }
        }

        WriteTypes = Lists.newArrayList();

        if(configBuilder.hasValue(WRITE_TYPES))
        {
            String writeTypesStr = configBuilder.getValue(WRITE_TYPES);

            if(writeTypesStr.equals(WriteType.ALL))
            {
                Arrays.stream(WriteType.values()).forEach(x -> WriteTypes.add(x));
            }
            else
            {
                String[] writeTypes = writeTypesStr.split(ITEM_DELIM, -1);
                Arrays.stream(writeTypes).forEach(x -> WriteTypes.add(WriteType.valueOf(x)));
            }
        }

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenome = loadRefGenome(RefGenomeFile);
        RefGenomeImageFile = configBuilder.getValue(REF_GENOME_IMAGE);

        BamStringency = ValidationStringency.STRICT;

        RefGenomeCoords = RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        VcfFile = configBuilder.getValue(OUTPUT_VCF);

        if(configBuilder.hasValue(OUTPUT_DIR))
        {
            OutputDir = parseOutputDir(configBuilder);
        }
        else
        {
            OutputDir = pathFromFile(VcfFile);
        }

        if(configBuilder.hasValue(HTML_SUMMARY_DIR))
        {
            HtmlOutputDir = configBuilder.getValue(HTML_SUMMARY_DIR);
        }
        else
        {
            HtmlOutputDir = OutputDir + DEFAULT_HTML_SUMMARY_DIR + File.separator;
        }

        PlotDiagrams = configBuilder.hasFlag(PLOT_DIAGRAMS);

        PerfDebug = configBuilder.hasFlag(PERF_DEBUG);
        PerfLogTime = configBuilder.getDecimal(PERF_LOG_TIME);
        OtherDebug = configBuilder.hasFlag(OTHER_DEBUG);

        Threads = parseThreads(configBuilder);
    }

    public String primaryBam() { return BamFiles.get(0); }

    // FIXME: decide on how represented
    public String tumorBam() { return BamFiles.get(0); }
    public String referenceBam() { return BamFiles.size() > 1 ? BamFiles.get(1) : null; }

    public String outputFilename(final WriteType writeType)
    {
        String filename = writeType == HTML_SUMMARY ? HtmlOutputDir : OutputDir;

        filename += SampleNames.get(0);

        switch(writeType)
        {
            case ASSEMBLY_BAM:
                filename += ASSEMBLY_BAM_FILE_ID;
                break;

            case BREAKEND_TSV:
                filename += BREAKEND_TSV;
                break;

            case HTML_SUMMARY:
                filename += "to_be_determined";
                break;
        }

        return filename;
    }


    public boolean writeHtmlFiles() { return WriteTypes.contains(HTML_SUMMARY); }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(TUMOR_BAM, false, TUMOR_BAM_DESC);
        configBuilder.addPath(REFERENCE_BAM, false, REFERENCE_BAM_DESC);
        configBuilder.addConfigItem(SAMPLE_IDS, true, "List of sample IDs, separated by ','");
        configBuilder.addConfigItem(SAMPLE_BAMS, false, "List of sample BAMs, separated by ','");
        configBuilder.addConfigItem(OUTPUT_VCF, true, "Output VCF filename");

        configBuilder.addPaths(
                JUNCTION_FILES, false, "List of SvPrep junction files, separated by ',', default is to match by sample name");

        addRefGenomeConfig(configBuilder, true);
        configBuilder.addPath(REF_GENOME_IMAGE, false, REFERENCE_BAM_DESC);
        configBuilder.addConfigItem(HTML_SUMMARY_DIR, false, "Directory for HTML summaries, default 'html'");
        configBuilder.addFlag(PLOT_DIAGRAMS, "Create HTML files containing SVGs");

        String writeTypes = Arrays.stream(WriteType.values()).map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM));
        configBuilder.addConfigItem(WRITE_TYPES, false, "Write types from list: " + writeTypes);

        configBuilder.addFlag(PERF_DEBUG, PERF_DEBUG_DESC);
        configBuilder.addDecimal(PERF_LOG_TIME, "Log performance data for routine exceeding specified time (0 = disabled)", 0);
        configBuilder.addFlag(OTHER_DEBUG, "Various other debugging");

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }

}
