package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.addCommonCommandOptions;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.checkFileExists;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.loadSpecificRegionsConfig;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.DEFAULT_CHR_PARTITION_SIZE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PERF_DEBUG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.common.CommonUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

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

    public SliceConfig(final ConfigBuilder configBuilder)
    {
        mIsValid = true;

        SampleId = configBuilder.getValue(SAMPLE);
        BamFile = configBuilder.getValue(BAM_FILE);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        WriteReads = configBuilder.hasFlag(WRITE_READS);
        WriteBam = configBuilder.hasFlag(WRITE_BAM) || !WriteReads;

        if(BamFile == null || OutputDir == null || RefGenomeFile == null)
        {
            BT_LOGGER.error("missing config: bam({}) refGenome({}) outputDir({})",
                    BamFile != null, RefGenomeFile != null, OutputDir != null);
            mIsValid = false;
        }

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        BT_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        BT_LOGGER.info("output({})", OutputDir);

        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);

        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();

        mIsValid &= loadSpecificRegionsConfig(configBuilder, SpecificChromosomes, SpecificRegions);

        if(SpecificRegions.isEmpty())
        {
            BT_LOGGER.error("missing specific regions or slice BED file for slicing");
            mIsValid = false;
        }

        Threads = parseThreads(configBuilder);

        PerfDebug = configBuilder.hasFlag(PERF_DEBUG);
    }

    public boolean isValid()
    {
        if(!mIsValid)
            return false;

        mIsValid = checkFileExists(BamFile) && checkFileExists(RefGenomeFile);
        return mIsValid;
    }

    public String formFilename(final String fileType) { return CommonUtils.formFilename(SampleId, BamFile, OutputDir, OutputId, fileType); }


    public static void addConfig(final ConfigBuilder configBuilder)
    {
        addCommonCommandOptions(configBuilder);

        configBuilder.addInteger(PARTITION_SIZE, "Partition size", DEFAULT_CHR_PARTITION_SIZE);
        configBuilder.addFlag(WRITE_BAM, "Write BAM file for sliced region");
        configBuilder.addFlag(WRITE_READS, "Write CSV reads file for sliced region");
        configBuilder.addFlag(PERF_DEBUG, "Detailed performance tracking and logging");
    }
}
