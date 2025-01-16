package com.hartwig.hmftools.bamtools.remapper;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.bwa.BwaUtils.LIBBWA_PATH;
import static com.hartwig.hmftools.common.bwa.BwaUtils.loadAlignerLibrary;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.utils.TaskExecutor.THREADS;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

public class AltContigRemapperConfig
{
    public final String mOutputFile;
    public final String mOrigBamFile;
    public final String mRefGenomeFile;
    public final RefGenomeVersion mRefGenVersion;
    public final String mBamToolPath;
    public final int mThreads;

    private static final String OUTPUT_FILE = "output_file";
    private static final String ORIG_BAM_FILE = "orig_bam_file";

    public AltContigRemapperConfig(final ConfigBuilder configBuilder)
    {
        mOutputFile = configBuilder.getValue(OUTPUT_FILE);
        mOrigBamFile = configBuilder.getValue(ORIG_BAM_FILE);
        mRefGenomeFile = configBuilder.getValue(REF_GENOME);
        mBamToolPath = configBuilder.getValue(BAMTOOL_PATH);
        mThreads = configBuilder.getInteger(THREADS);

        if(mOrigBamFile == null || mOutputFile == null)
        {
            BT_LOGGER.error("missing config: bam(orig={} new={})",
                    mOrigBamFile != null, mOutputFile != null);
            System.exit(1);
        }
        loadAlignerLibrary(null);

        mRefGenVersion = BamUtils.deriveRefGenomeVersion(mOrigBamFile);

        BT_LOGGER.info("origBam({}) outputFile({})", mOrigBamFile, mOutputFile);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output comparison file");
        configBuilder.addPath(ORIG_BAM_FILE, true, "Original BAM file");
        configBuilder.addPath(LIBBWA_PATH, false, "Path to BWA library");
        configBuilder.addPath(BAMTOOL_PATH, true, "Path to BWA library");

        addRefGenomeFile(configBuilder, true);
        addThreadOptions(configBuilder, 1);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }

    public PairAligner aligner()
    {
        return new BwaPairAligner(mRefGenomeFile);
    }

    @VisibleForTesting
    public AltContigRemapperConfig(String origBamFile, String outputFile, String bamToolPath)
    {
        mOrigBamFile = origBamFile;
        mOutputFile = outputFile;
        mBamToolPath = bamToolPath;
        mRefGenomeFile = null;
        mRefGenVersion = RefGenomeVersion.V38;
        mThreads = 1;
    }
}
