package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class TeloConfig
{
    public final int ThreadCount;

    public final String SampleId;

    public final String SampleBamFile;
    public final String RefGenomeFile;
    public final String OutputDir;
    public final boolean WriteReads;

    public static String SAMPLE_ID = "sample";

    private static String THREADS = "threads";
    private static String SAMPLE_BAM_FILE = "sample_bam";
    private static String REF_GENOME = "ref_genome";
    private static String WRITE_READS = "write_reads";

    public static final Logger TE_LOGGER = LogManager.getLogger(com.hartwig.hmftools.telo.TeloConfig.class);

    public TeloConfig(final CommandLine cmd)
    {
        ThreadCount = getConfigValue(cmd, THREADS, 1);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME, "");

        SampleId = cmd.getOptionValue(SAMPLE_ID);
        SampleBamFile = cmd.getOptionValue(SAMPLE_BAM_FILE);
        OutputDir = parseOutputDir(cmd);

        WriteReads = cmd.hasOption(WRITE_READS);
    }

    public boolean isValid()
    {
        if(SampleBamFile == null || SampleId == null)
        {
            TE_LOGGER.error( "missing sampleId or BAM file");
            return false;
        }

        if(!checkCreateOutputDir(OutputDir))
        {
            TE_LOGGER.error("failed to create output directory({})", OutputDir);
            return false;
        }

        return true;
    }

    public static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(THREADS, true, "Number of threads");
        options.addOption(SAMPLE_ID, true, "Name of tumor sample");
        options.addOption(SAMPLE_BAM_FILE, true, "Path to tumor bam file");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(REF_GENOME, true, "Path to reference genome fasta file if using CRAM files");
        options.addOption(WRITE_READS, false, "Write BAM read data to CSV file");

        return options;
    }
}
