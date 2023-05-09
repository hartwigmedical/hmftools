package com.hartwig.hmftools.markdups;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.samtools.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_CHROMOSOMES;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificChromsomesOrRegions;
import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_PARTITION_SIZE;
import static com.hartwig.hmftools.markdups.common.Constants.DEFAULT_POS_BUFFER_SIZE;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.samtools.BamUtils;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.markdups.common.FilterReadsType;
import com.hartwig.hmftools.markdups.umi.UmiConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class MarkDupsConfig
{
    public final String SampleId;
    public final String BamFile;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;
    public final RefGenomeInterface RefGenome;

    public final int PartitionSize;
    public final int BufferSize;
    public final ValidationStringency BamStringency;

    // UMI group config
    public final UmiConfig UMIs;

    public final String OutputDir;
    public final String OutputId;
    public final boolean WriteBam;
    public final boolean NoMateCigar;
    public final int Threads;

    // debug
    public final List<String> SpecificChromosomes;
    public final List<String> LogReadIds;
    public final List<ChrBaseRegion> SpecificRegions;
    public final FilterReadsType SpecificRegionsFilterType;
    public final ReadOutput LogReadType;
    public final boolean PerfDebug;
    public final boolean RunChecks;
    public final boolean WriteStats;

    private boolean mIsValid;

    public static final Logger MD_LOGGER = LogManager.getLogger(MarkDupsConfig.class);

    // config strings
    private static final String SAMPLE = "sample";
    private  static final String BAM_FILE = "bam_file";
    private static final String PARTITION_SIZE = "partition_size";
    private static final String BUFFER_SIZE = "buffer_size";
    private static final String READ_OUTPUTS = "read_output";
    private static final String NO_MATE_CIGAR = "no_mate_cigar";
    private static final String WRITE_BAM = "write_bam";

    private static final String LOG_READ_IDS = "log_read_ids";
    private static final String PERF_DEBUG = "perf_debug";
    private static final String RUN_CHECKS = "run_checks";
    private static final String WRITE_STATS = "write_stats";
    private static final String SPECIFIC_REGION_FILTER_TYPE = "specific_region_filter";

    private static final String ITEM_DELIM = ";";

    public MarkDupsConfig(final CommandLine cmd)
    {
        mIsValid = true;
        SampleId = cmd.getOptionValue(SAMPLE);
        BamFile = cmd.getOptionValue(BAM_FILE);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME);
        RefGenome = loadRefGenome(RefGenomeFile);
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);

        if(SampleId == null || BamFile == null || OutputDir == null || RefGenomeFile == null)
        {
            MD_LOGGER.error("missing config: sample({}) bam({}) refGenome({}) outputDir({})",
                    SampleId != null, BamFile != null, RefGenomeFile != null, OutputDir != null);
            mIsValid = false;
        }

        RefGenVersion = RefGenomeVersion.from(cmd);

        MD_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        MD_LOGGER.info("output({})", OutputDir);

        PartitionSize = Integer.parseInt(cmd.getOptionValue(PARTITION_SIZE, String.valueOf(DEFAULT_PARTITION_SIZE)));
        BufferSize = Integer.parseInt(cmd.getOptionValue(BUFFER_SIZE, String.valueOf(DEFAULT_POS_BUFFER_SIZE)));
        NoMateCigar = cmd.hasOption(NO_MATE_CIGAR);
        BamStringency = BamUtils.validationStringency(cmd);

        UMIs = UmiConfig.from(cmd);

        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();

        try
        {
            loadSpecificChromsomesOrRegions(cmd, SpecificChromosomes, SpecificRegions, MD_LOGGER);
            Collections.sort(SpecificRegions);
        }
        catch(ParseException e)
        {
            MD_LOGGER.error("invalid specific regions({}) chromosomes({}) config",
                    cmd.getOptionValue(SPECIFIC_REGIONS, ""), cmd.getOptionValue(SPECIFIC_CHROMOSOMES, ""));
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

        WriteStats = cmd.hasOption(WRITE_STATS);
        PerfDebug = cmd.hasOption(PERF_DEBUG);
        RunChecks = cmd.hasOption(RUN_CHECKS);

        if(RunChecks || UMIs.Debug || UMIs.HighlightConsensus)
        {
            MD_LOGGER.info("running debug options: read-checks({}) umi-validation({}) consensus-highlight({})",
                    RunChecks, UMIs.Debug, UMIs.HighlightConsensus);
        }
    }

    public boolean isValid()
    {
        if(!mIsValid)
            return false;

        if(!Files.exists(Paths.get(BamFile)))
        {
            MD_LOGGER.error("invalid bam file path: {}", BamFile);
            return false;
        }

        if(!Files.exists(Paths.get(RefGenomeFile)))
        {
            MD_LOGGER.error("invalid ref genome file: {}", RefGenomeFile);
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
        options.addOption(NO_MATE_CIGAR, false, "Mate CIGAR not set by aligner, make no attempt to use it");
        addValidationStringencyOption(options);
        UmiConfig.addCommandLineOptions(options);
        addThreadOptions(options);

        addSpecificChromosomesRegionsConfig(options);
        options.addOption(LOG_READ_IDS, true, "Log specific read IDs, separated by ';'");
        options.addOption(PERF_DEBUG, false, "Detailed performance tracking and logging");
        options.addOption(RUN_CHECKS, false, "Run duplicate mismatch checks");
        options.addOption(WRITE_STATS, false, "Write duplicate and UMI-group stats");
        options.addOption(SPECIFIC_REGION_FILTER_TYPE, true, "Used with specific regions, to filter mates or supps");

        return options;
    }

    public MarkDupsConfig(int partitionSize, int bufferSize, final RefGenomeInterface refGenome, boolean umiEnabled)
    {
        mIsValid = true;
        SampleId = "";
        BamFile = null;
        RefGenomeFile = null;
        OutputDir = null;
        OutputId = "";
        RefGenVersion = V37;
        RefGenome = refGenome;

        PartitionSize = partitionSize;
        BufferSize = bufferSize;
        UMIs = new UmiConfig(umiEnabled);
        NoMateCigar = false;
        BamStringency = ValidationStringency.STRICT;

        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();
        SpecificRegionsFilterType = FilterReadsType.MATE_AND_SUPP;

        WriteBam = false;
        LogReadType = ReadOutput.NONE;

        LogReadIds = Lists.newArrayList();
        Threads = 0;
        PerfDebug = false;
        RunChecks = true;
        WriteStats = false;
    }
}
