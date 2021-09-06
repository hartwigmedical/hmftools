package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.ConfigUtils;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class TeloConfig
{
    public final int ThreadCount;

    public final String BamFile;
    public final String RefGenomeFile;
    public final String OutputDir;

    public final List<String> SpecificChromosomes;
    public final int Threads;

    private static String THREADS = "threads";
    private static String BAM_FILE = "bam_file";
    private static String REF_GENOME = "ref_genome";
    private static final String SPECIFIC_CHR = "specific_chr";

    public static final Logger TE_LOGGER = LogManager.getLogger(com.hartwig.hmftools.telo.TeloConfig.class);

    public TeloConfig(final CommandLine cmd)
    {
        ThreadCount = getConfigValue(cmd, THREADS, 1);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME, "");

        BamFile = cmd.getOptionValue(BAM_FILE);
        OutputDir = parseOutputDir(cmd);

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));

        SpecificChromosomes = Lists.newArrayList();

        if(cmd.hasOption(SPECIFIC_CHR))
        {
            final String chromosomes = cmd.getOptionValue(SPECIFIC_CHR);
            TE_LOGGER.info("filtering for specific chromosomes: {}", chromosomes);
            SpecificChromosomes.addAll(Arrays.stream(chromosomes.split(";")).collect(Collectors.toList()));
        }
    }

    public boolean isValid()
    {
        if(BamFile == null)
        {
            TE_LOGGER.error( "missing BAM / CRAM file");
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
        options.addOption(BAM_FILE, true, "Path to bam file");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(REF_GENOME, true, "Path to reference genome fasta file if using CRAM files");
        options.addOption(SPECIFIC_CHR, true, "Optional: list of chromosomes separated by ;");
        ConfigUtils.addLoggingOptions(options);

        return options;
    }
}
