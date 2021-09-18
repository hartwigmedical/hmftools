package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.ConfigUtils;
import org.apache.commons.cli.ParseException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class TeloConfig
{
    public final int ThreadCount;

    public final String SampleId;

    public final String BamFile;
    public final String RefGenomeFile;
    public final String OutputDir;

    public final String SampleType;
    public final List<String> SpecificChromosomes;

    private static String SAMPLE_ID = "sample_id";
    private static String BAM_FILE = "bam_file";
    private static String SAMPLE_TYPE = "sample_type";
    private static String REF_GENOME = "ref_genome";
    private static String THREADS = "threads";
    private static final String SPECIFIC_CHR = "specific_chr";

    public static final Logger TE_LOGGER = LogManager.getLogger(com.hartwig.hmftools.telo.TeloConfig.class);

    public TeloConfig(final CommandLine cmd) throws ParseException
    {
        RefGenomeFile = cmd.getOptionValue(REF_GENOME);

        SampleId = cmd.getOptionValue(SAMPLE_ID);
        BamFile = cmd.getOptionValue(BAM_FILE);
        OutputDir = parseOutputDir(cmd);
        SampleType = cmd.getOptionValue(SAMPLE_TYPE);

        try
        {
            ThreadCount = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));
        }
        catch (NumberFormatException e)
        {
            throw new ParseException(String.format("unable to parse %s arg: %s", THREADS, cmd.getOptionValue(THREADS)));
        }

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
        options.addOption(Option.builder(SAMPLE_ID).hasArg().required().desc("ID of tumor sample").build());
        options.addOption(Option.builder(BAM_FILE).hasArg().required().desc("Path to bam/cram file").build());
        options.addOption(Option.builder(OUTPUT_DIR).hasArg().required().desc("Output directory").build());
        options.addOption(Option.builder(SAMPLE_TYPE).hasArg().required().desc("Type of sample (germline / somatic)").build());
        options.addOption(Option.builder(THREADS).hasArg().type(Number.class).desc("Number of bam reader threads (default = 1)").build());
        options.addOption(Option.builder(REF_GENOME).hasArg().desc("Path to reference genome fasta file if using CRAM files").build());
        options.addOption(Option.builder(SPECIFIC_CHR).hasArg().desc("Optional: list of chromosomes separated by ;").build());
        ConfigUtils.addLoggingOptions(options);

        return options;
    }
}
