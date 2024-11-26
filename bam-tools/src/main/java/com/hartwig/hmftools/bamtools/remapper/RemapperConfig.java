package com.hartwig.hmftools.bamtools.remapper;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bwa.BwaUtils.LIBBWA_PATH;
import static com.hartwig.hmftools.common.bwa.BwaUtils.loadAlignerLibrary;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

public class RemapperConfig
{
    public final String OutputFile;
    public final String OrigBamFile;
    public final String RefGenomeFile;
//    public final RefGenomeVersion RefGenVersion;

    private static final String OUTPUT_FILE = "output_file";
    private static final String ORIG_BAM_FILE = "orig_bam_file";

    public RemapperConfig(final ConfigBuilder configBuilder)
    {
        OutputFile =  configBuilder.getValue(OUTPUT_FILE);
        OrigBamFile =  configBuilder.getValue(ORIG_BAM_FILE);
        RefGenomeFile =  configBuilder.getValue(REF_GENOME);

        if(OrigBamFile == null || OutputFile == null)
        {
            BT_LOGGER.error("missing config: bam(orig={} new={})",
                    OrigBamFile != null, OutputFile != null);
            System.exit(1);
        }
        loadAlignerLibrary(null);


//        RefGenVersion = BamUtils.deriveRefGenomeVersion(OrigBamFile);

        BT_LOGGER.info("origBam({}) outputFile({})", OrigBamFile, OutputFile);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output comparison file");
        configBuilder.addPath(ORIG_BAM_FILE, true, "Original BAM file");
        configBuilder.addPath(LIBBWA_PATH, false, "Path to BWA library");

        addRefGenomeFile(configBuilder, true);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }

    @VisibleForTesting
    public RemapperConfig()
    {
        OutputFile = null;
        OrigBamFile = null;
        RefGenomeFile = null;
//        RefGenVersion = null;
    }
}
