package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.SAMPLE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.addCommonCommandOptions;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.checkFileExists;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.loadSpecificRegionsConfig;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.metrics.MetricsConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.bamtools.metrics.MetricsConfig.PERF_DEBUG;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.DEFAULT_CHR_PARTITION_SIZE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.common.CommonUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class SliceConfig
{
    public final String SampleId;
    public final String BamFile;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final int PartitionSize;

    public final String OutputDir;
    public final String OutputId;

    public final boolean WriteBam;
    public final boolean WriteReads;
    public final int Threads;

    // debug
    public final List<String> SpecificChromosomes;
    public final List<ChrBaseRegion> SpecificRegions;
    public final boolean PerfDebug;

    private boolean mIsValid;

    private static final String WRITE_BAM = "write_bam";
    private static final String WRITE_READS = "write_reads";

    public SliceConfig(final CommandLine cmd)
    {
        mIsValid = true;

        SampleId = cmd.getOptionValue(SAMPLE);
        BamFile = cmd.getOptionValue(BAM_FILE);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME);
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);
        WriteReads = cmd.hasOption(WRITE_READS);
        WriteBam = cmd.hasOption(WRITE_BAM) || !WriteReads;

        if(BamFile == null || OutputDir == null || RefGenomeFile == null)
        {
            BT_LOGGER.error("missing config: bam({}) refGenome({}) outputDir({})",
                    BamFile != null, RefGenomeFile != null, OutputDir != null);
            mIsValid = false;
        }

        RefGenVersion = RefGenomeVersion.from(cmd);

        BT_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        BT_LOGGER.info("output({})", OutputDir);

        PartitionSize = Integer.parseInt(cmd.getOptionValue(PARTITION_SIZE, String.valueOf(DEFAULT_CHR_PARTITION_SIZE)));

        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();

        mIsValid &= loadSpecificRegionsConfig(cmd, SpecificChromosomes, SpecificRegions);

        if(SpecificRegions.isEmpty())
        {
            BT_LOGGER.error("missing specific regions or slice BED file for slicing");
            mIsValid = false;
        }

        Threads = parseThreads(cmd);

        PerfDebug = cmd.hasOption(PERF_DEBUG);
    }

    public boolean isValid()
    {
        if(!mIsValid)
            return false;

        mIsValid = checkFileExists(BamFile) && checkFileExists(RefGenomeFile);
        return mIsValid;
    }

    public String formFilename(final String fileType) { return CommonUtils.formFilename(SampleId, BamFile, OutputDir, OutputId, fileType); }


    public static Options createCmdLineOptions()
    {
        final Options options = new Options();

        addCommonCommandOptions(options);

        options.addOption(PARTITION_SIZE, true, "Partition size, default: " + DEFAULT_CHR_PARTITION_SIZE);
        options.addOption(WRITE_BAM, false, "Write BAM file for sliced region");
        options.addOption(WRITE_READS, false, "Write CSV reads file for sliced region");

        options.addOption(LOG_READ_IDS, true, "Log specific read IDs, separated by ';'");

        return options;
    }
}
