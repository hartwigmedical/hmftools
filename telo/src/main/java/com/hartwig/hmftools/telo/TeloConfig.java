package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;
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

    public final List<String> SpecificChromosomes;
    public final int Threads;

    public static String SAMPLE_ID = "sample";

    private static String THREADS = "threads";
    private static String SAMPLE_BAM_FILE = "sample_bam";
    private static String REF_GENOME = "ref_genome";
    private static String WRITE_READS = "write_reads";
    private static final String SPECIFIC_CHR = "specific_chr";

    public static final Logger TE_LOGGER = LogManager.getLogger(com.hartwig.hmftools.telo.TeloConfig.class);

    public TeloConfig(final CommandLine cmd)
    {
        ThreadCount = getConfigValue(cmd, THREADS, 1);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME, "");

        SampleId = cmd.getOptionValue(SAMPLE_ID);
        SampleBamFile = cmd.getOptionValue(SAMPLE_BAM_FILE);
        OutputDir = parseOutputDir(cmd);

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));

        SpecificChromosomes = Lists.newArrayList();

        if(cmd.hasOption(SPECIFIC_CHR))
        {
            final String chromosomes = cmd.getOptionValue(SPECIFIC_CHR);
            TE_LOGGER.info("filtering for specific chromosomes: {}", chromosomes);
            SpecificChromosomes.addAll(Arrays.stream(chromosomes.split(";")).collect(Collectors.toList()));
        }

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
        options.addOption(SPECIFIC_CHR, true, "Optional: list of chromosomes separated by ;");
        options.addOption(WRITE_READS, false, "Write BAM read data to CSV file");
        options.addOption(LOG_DEBUG, false, "Log verbose");

        return options;
    }
}
