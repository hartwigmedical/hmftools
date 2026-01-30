package com.hartwig.hmftools.bamtools.remapper;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.perf.TaskExecutor.THREADS;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.umccr.java.hellbender.utils.bwa.BwaMemIndex;

import htsjdk.samtools.SAMFileHeader;

public class AltContigRemapperConfig
{
    public final String OutputFile;
    public final String OrigBamFile;
    public final String RefGenomeFile;
    public final boolean SliceHlaRegionsOnly;
    public final RefGenomeVersion RefGenVersion;
    public final String BamToolPath;
    public final int Threads;

    private static final String OUTPUT_FILE = "output_file";
    private static final String ORIG_BAM_FILE = "orig_bam_file";
    private static final String HLA_REGIONS_ONLY = "hla_regions_only";
    private static final String REF_GENOME_IMAGE_EXTENSION = ".img";

    public AltContigRemapperConfig(final ConfigBuilder configBuilder)
    {
        OutputFile = configBuilder.getValue(OUTPUT_FILE);
        OrigBamFile = configBuilder.getValue(ORIG_BAM_FILE);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        SliceHlaRegionsOnly = configBuilder.hasFlag(HLA_REGIONS_ONLY);
        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);
        Threads = configBuilder.getInteger(THREADS);

        if(OrigBamFile == null || OutputFile == null)
        {
            BT_LOGGER.error("missing config: bam(orig={} new={})",
                    OrigBamFile != null, OutputFile != null);
            System.exit(1);
        }

        RefGenVersion = BamUtils.deriveRefGenomeVersion(OrigBamFile);

        BT_LOGGER.info("origBam({}) outputFile({})", OrigBamFile);
        BT_LOGGER.info("outputBam({})", OutputFile);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output comparison file");
        configBuilder.addPath(ORIG_BAM_FILE, true, "Original BAM file");
        configBuilder.addFlag(HLA_REGIONS_ONLY, "Slice HLA regions only");
        configBuilder.addPath(BAMTOOL_PATH, false, "Path to BWA library");

        addRefGenomeFile(configBuilder, true);
        addThreadOptions(configBuilder, 1);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }

    public PairAligner pairAligner()
    {
        return new BwaPairAligner(bwaMemIndex());
    }

    public SingleRecordAligner singleRecordAligner(SAMFileHeader header)
    {
        return new BwaSingleRecordAligner(bwaMemIndex(), header);
    }

    private BwaMemIndex bwaMemIndex()
    {
        try
        {
            String imgFileName = RefGenomeFile + REF_GENOME_IMAGE_EXTENSION;
            return new BwaMemIndex(imgFileName);
        }
        catch(Exception e)
        {
            BT_LOGGER.error("failed to initialise BWA aligner: {}", e.toString());
            return null;
        }
    }

    @VisibleForTesting
    public AltContigRemapperConfig(String origBamFile, String outputFile, String bamToolPath, boolean sliceHlaRegionsOnly)
    {
        OrigBamFile = origBamFile;
        OutputFile = outputFile;
        BamToolPath = bamToolPath;
        SliceHlaRegionsOnly = sliceHlaRegionsOnly;
        RefGenomeFile = null;
        RefGenVersion = RefGenomeVersion.V38;
        Threads = 1;
    }
}
