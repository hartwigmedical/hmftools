package com.hartwig.hmftools.bamtools.markdups;

import static com.hartwig.hmftools.bamtools.BmConfig.BAM_FILE;
import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.BmConfig.ITEM_DELIM;
import static com.hartwig.hmftools.bamtools.BmConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.bamtools.BmConfig.PARTITION_SIZE;
import static com.hartwig.hmftools.bamtools.BmConfig.PERF_DEBUG;
import static com.hartwig.hmftools.bamtools.BmConfig.SAMPLE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificChromsomesOrRegions;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class MarkDupsConfig
{
    public final String SampleId;
    public final String BamFile;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final int PartitionSize;
    public final int BufferSize;

    public final String OutputDir;
    public final String OutputId;
    public final boolean WriteBam;
    public final boolean UseInterimFiles;
    public final int Threads;

    // debug
    public final List<String> SpecificChromosomes;
    public final List<String> LogReadIds;
    public final List<ChrBaseRegion> SpecificRegions;
    public final FilterReadsType SpecificRegionsFilterType;
    public final ReadOutput LogReadType;
    public final boolean PerfDebug;
    public final boolean RunChecks;

    private boolean mIsValid;

    // config strings
    private static final String BUFFER_SIZE = "buffer_size";
    private static final String READ_OUTPUTS = "read_output";
    private static final String WRITE_BAM = "write_bam";
    private static final String RUN_CHECKS = "run_checks";
    private static final String USE_INTERIM_FILES = "use_interim_files";
    private static final String SPECIFIC_REGION_FILTER_TYPE = "specific_region_filter";

    private static final int DEFAULT_PARTITION_SIZE = 1000000;
    private static final int DEFAULT_POS_BUFFER_SIZE = 1000;

    public MarkDupsConfig(final CommandLine cmd)
    {
        mIsValid = true;
        SampleId = cmd.getOptionValue(SAMPLE);
        BamFile = cmd.getOptionValue(BAM_FILE);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME);
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);

        if(SampleId == null || BamFile == null || OutputDir == null || RefGenomeFile == null)
        {
            BM_LOGGER.error("missing config: sample({}) bam({}) refGenome({}) outputDir({})",
                    SampleId != null, BamFile != null, RefGenomeFile != null, OutputDir != null);
            mIsValid = false;
        }

        RefGenVersion = cmd.hasOption(REF_GENOME_VERSION) ? RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION)) : V37;

        BM_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        BM_LOGGER.info("output({})", OutputDir);

        PartitionSize = Integer.parseInt(cmd.getOptionValue(PARTITION_SIZE, String.valueOf(DEFAULT_PARTITION_SIZE)));
        BufferSize = Integer.parseInt(cmd.getOptionValue(BUFFER_SIZE, String.valueOf(DEFAULT_POS_BUFFER_SIZE)));

        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();

        try
        {
            loadSpecificChromsomesOrRegions(cmd, SpecificChromosomes, SpecificRegions, BM_LOGGER);
        }
        catch(ParseException e)
        {
            mIsValid = false;
        }

        SpecificRegionsFilterType = !SpecificChromosomes.isEmpty() || !SpecificRegions.isEmpty() ?
                FilterReadsType.valueOf(cmd.getOptionValue(SPECIFIC_REGION_FILTER_TYPE, FilterReadsType.READ.toString())) :
                FilterReadsType.NONE;

        WriteBam = cmd.hasOption(WRITE_BAM) || !cmd.hasOption(READ_OUTPUTS);
        LogReadType = ReadOutput.valueOf(cmd.getOptionValue(READ_OUTPUTS, ReadOutput.NONE.toString()));

        LogReadIds = cmd.hasOption(LOG_READ_IDS) ?
                Arrays.stream(cmd.getOptionValue(LOG_READ_IDS).split(ITEM_DELIM, -1)).collect(Collectors.toList()) : Lists.newArrayList();

        Threads = parseThreads(cmd);

        PerfDebug = cmd.hasOption(PERF_DEBUG);
        RunChecks = cmd.hasOption(RUN_CHECKS);
        UseInterimFiles = cmd.hasOption(USE_INTERIM_FILES);
    }

    public boolean isValid()
    {
        if(!mIsValid)
            return false;

        if(!Files.exists(Paths.get(BamFile)))
        {
            BM_LOGGER.error("invalid bam file path: {}", BamFile);
            return false;
        }

        if(!Files.exists(Paths.get(RefGenomeFile)))
        {
            BM_LOGGER.error("invalid ref genome file: {}", RefGenomeFile);
            return false;
        }

        return true;
    }

    public String formFilename(final String fileType)
    {
        String filename = OutputDir + SampleId;

        filename += "." + fileType;

        if(OutputId != null)
            filename += "." + OutputId;

        filename += ".csv";

        return filename;
    }

    public boolean runReadChecks() { return RunChecks && (!SpecificChromosomes.isEmpty() || !SpecificRegions.isEmpty()); }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        addOutputOptions(options);
        addLoggingOptions(options);

        options.addOption(SAMPLE, true, "Tumor sample ID");
        options.addOption(BAM_FILE, true, "BAM file location");
        addRefGenomeConfig(options);;
        options.addOption(PARTITION_SIZE, true, "Partition size, default: " + DEFAULT_PARTITION_SIZE);
        options.addOption(BUFFER_SIZE, true, "Read buffer size, default: " + DEFAULT_POS_BUFFER_SIZE);
        options.addOption(READ_OUTPUTS, true, "Write reads: NONE (default), 'MISMATCHES', 'DUPLICATES', 'ALL'");
        options.addOption(WRITE_BAM, false, "Write BAM, default true if not write read output");
        addThreadOptions(options);

        addSpecificChromosomesRegionsConfig(options);
        options.addOption(LOG_READ_IDS, true, "Log specific read IDs, separated by ';'");
        options.addOption(PERF_DEBUG, false, "Detailed performance tracking and logging");
        options.addOption(RUN_CHECKS, false, "Run duplicate mismatch checks");
        options.addOption(USE_INTERIM_FILES, false, "Write candidate duplicate reads to file");
        options.addOption(SPECIFIC_REGION_FILTER_TYPE, true, "Used with specific regions, to filter mates or supps");

        return options;
    }

    public MarkDupsConfig()
    {
        this(DEFAULT_PARTITION_SIZE, DEFAULT_POS_BUFFER_SIZE);
    }

    public MarkDupsConfig(int partitionSize, int bufferSize)
    {
        mIsValid = true;
        SampleId = "";
        BamFile = null;
        RefGenomeFile = null;
        OutputDir = null;
        OutputId = "";
        RefGenVersion = V37;

        PartitionSize = partitionSize;
        BufferSize = bufferSize;

        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();
        SpecificRegionsFilterType = FilterReadsType.MATE_AND_SUPP;

        WriteBam = false;
        LogReadType = ReadOutput.NONE;

        LogReadIds = Lists.newArrayList();
        Threads = 0;
        PerfDebug = false;
        RunChecks = false;
        UseInterimFiles = false;
    }
}
