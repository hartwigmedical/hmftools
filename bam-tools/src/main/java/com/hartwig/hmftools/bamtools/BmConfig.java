package com.hartwig.hmftools.bamtools;

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

import java.io.IOException;
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
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BmConfig
{
    public final String SampleId;
    public final String BamFile;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final int MapQualityThreshold;
    public final int BaseQualityThreshold;
    public final int MaxCoverage;

    public final int PartitionSize;

    // metrics capture config
    public final boolean ExcludeZeroCoverage;
    public final boolean WriteOldStyle;

    public final String OutputDir;
    public final String OutputId;

    public final int Threads;

    // debug
    public final List<String> SpecificChromosomes;
    public final List<String> LogReadIds;
    public final List<ChrBaseRegion> SpecificRegions;
    public final boolean PerfDebug;

    private final List<BamFunction> mFunctions;
    private boolean mIsValid;

    public static final Logger BM_LOGGER = LogManager.getLogger(BmConfig.class);

    // config strings
    public static final String SAMPLE = "sample";
    private static final String BAM_FILE = "bam_file";
    private static final String FUNCTIONS = "functions";
    private static final String PARTITION_SIZE = "partition_size";
    private static final String MAP_QUAL_THRESHOLD = "map_qual_threshold";
    private static final String BASE_QUAL_THRESHOLD = "base_qual_threshold";
    private static final String MAX_COVERAGE = "max_coverage";
    private static final String EXCLUDE_ZERO_COVERAGE = "exclude_zero_coverage";
    private static final String WRITE_OLD_STYLE = "write_old_style";
    private static final String SLICE_BED = "slice_bed";

    private static final String LOG_READ_IDS = "log_read_ids";
    private static final String PERF_DEBUG = "perf_debug";

    public static final String ITEM_DELIM = ";";

    private static final int DEFAULT_CHR_PARTITION_SIZE = 1000000;
    private static final int DEFAULT_MAP_QUAL_THRESHOLD = 20;
    private static final int DEFAULT_BASE_QUAL_THRESHOLD = 10;
    private static final int DEFAULT_MAX_COVERAGE = 250;

    public BmConfig(final CommandLine cmd)
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

        mFunctions = Lists.newArrayList();

        if(cmd.hasOption(FUNCTIONS))
        {
            Arrays.stream(cmd.getOptionValue(FUNCTIONS).split(ITEM_DELIM)).forEach(x -> mFunctions.add(BamFunction.valueOf(x)));
        }
        else
        {
            mFunctions.add(BamFunction.METRICS);
        }

        BM_LOGGER.info("running functions: {}", mFunctions);

        PartitionSize = Integer.parseInt(cmd.getOptionValue(PARTITION_SIZE, String.valueOf(DEFAULT_CHR_PARTITION_SIZE)));
        MapQualityThreshold = Integer.parseInt(cmd.getOptionValue(MAP_QUAL_THRESHOLD, String.valueOf(DEFAULT_MAP_QUAL_THRESHOLD)));
        BaseQualityThreshold = Integer.parseInt(cmd.getOptionValue(BASE_QUAL_THRESHOLD, String.valueOf(DEFAULT_BASE_QUAL_THRESHOLD)));
        MaxCoverage = Integer.parseInt(cmd.getOptionValue(MAX_COVERAGE, String.valueOf(DEFAULT_MAX_COVERAGE)));
        ExcludeZeroCoverage = cmd.hasOption(EXCLUDE_ZERO_COVERAGE);
        WriteOldStyle = cmd.hasOption(WRITE_OLD_STYLE);

        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();

        if(runSlicing() && cmd.hasOption(SLICE_BED))
        {
            loadSliceRegions(cmd.getOptionValue(SLICE_BED));
        }
        else
        {
            try
            {
                loadSpecificChromsomesOrRegions(cmd, SpecificChromosomes, SpecificRegions, BM_LOGGER);

                if(runSlicing())
                    SpecificChromosomes.clear();
            }
            catch(ParseException e)
            {
                mIsValid = false;
            }
        }

        if(runSlicing() && SpecificRegions.isEmpty())
        {
            BM_LOGGER.error("missing specific regions or slice BED file for slicing");
            mIsValid = false;
        }

        LogReadIds = cmd.hasOption(LOG_READ_IDS) ?
                Arrays.stream(cmd.getOptionValue(LOG_READ_IDS).split(ITEM_DELIM, -1)).collect(Collectors.toList()) : Lists.newArrayList();

        Threads = parseThreads(cmd);

        PerfDebug = cmd.hasOption(PERF_DEBUG);
    }

    public boolean runFunction(final BamFunction function) { return mFunctions.contains(function); }
    public boolean runSlicing() { return runFunction(BamFunction.BAM_SLICE) || runFunction(BamFunction.BAM_READS); }
    public boolean runMetrics() { return runFunction(BamFunction.METRICS); }

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

        filename += ".bam_" + fileType;

        if(OutputId != null)
            filename += "." + OutputId;

        filename += ".csv";

        return filename;
    }

    private void loadSliceRegions(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));

            for(String line : lines)
            {
                if(line.contains("Chromosome"))
                    continue;

                final String[] values = line.split("\t", -1);

                if(values.length < 3)
                {
                    BM_LOGGER.error("invalid slice BED entry: {}", line);
                    mIsValid = false;
                    return;
                }

                String chromosome = values[0];
                int posStart = Integer.parseInt(values[1]) + 1;
                int posEnd = Integer.parseInt(values[2]);
                SpecificRegions.add(new ChrBaseRegion(chromosome, posStart, posEnd));
            }

            BM_LOGGER.info("loaded {} slice regions from file", SpecificRegions.size(), filename);
        }
        catch(IOException e)
        {
            BM_LOGGER.error("failed to load hotspot data file({}): {}", filename, e.toString());
            mIsValid = false;
        }
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        addOutputOptions(options);
        addLoggingOptions(options);

        options.addOption(SAMPLE, true, "Tumor sample ID");
        options.addOption(BAM_FILE, true, "BAM file location");
        addRefGenomeConfig(options);;
        options.addOption(FUNCTIONS, true, "Functions to run: {}" + BamFunction.values());
        options.addOption(SLICE_BED, true, "BED file defining region to slice");
        options.addOption(PARTITION_SIZE, true, "Partition size, default: " + DEFAULT_CHR_PARTITION_SIZE);
        options.addOption(MAP_QUAL_THRESHOLD, true, "Map quality threshold, default: " + DEFAULT_MAP_QUAL_THRESHOLD);
        options.addOption(BASE_QUAL_THRESHOLD, true, "Base quality threshold, default: " + DEFAULT_BASE_QUAL_THRESHOLD);
        options.addOption(MAX_COVERAGE, true, "Max coverage, default: " + DEFAULT_MAX_COVERAGE);
        options.addOption(EXCLUDE_ZERO_COVERAGE, false, "Exclude bases with zero coverage");
        options.addOption(WRITE_OLD_STYLE, false, "Write data in same format as Picard CollectWgsMetrics");
        addThreadOptions(options);

        addSpecificChromosomesRegionsConfig(options);
        options.addOption(LOG_READ_IDS, true, "Log specific read IDs, separated by ';'");
        options.addOption(PERF_DEBUG, false, "Detailed performance tracking and logging");

        return options;
    }
}
