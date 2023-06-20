package com.hartwig.hmftools.svprep.append;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.samtools.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificRegions;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConfig.BAM_FILE;
import static com.hartwig.hmftools.svprep.SvConfig.LOG_READ_IDS;
import static com.hartwig.hmftools.svprep.SvConfig.SAMPLE;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.reads.ReadFilterConfig;
import com.hartwig.hmftools.svprep.reads.ReadFilters;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class AppendConfig
{
    public final String SampleId;
    public final String BamFile;
    public final String InputVcf;
    public final String OutputVcf;

    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final ReadFilters ReadFiltering;

    public final int Threads;

    // debug
    public final List<String> LogReadIds;
    public final List<ChrBaseRegion> SpecificRegions;

    public static final String INPUT_VCF = "input_vcf";
    public static final String OUTPUT_VCF = "output_vcf";

    public AppendConfig(final CommandLine cmd)
    {
        SampleId = cmd.getOptionValue(SAMPLE);
        BamFile = cmd.getOptionValue(BAM_FILE);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME);
        InputVcf = cmd.getOptionValue(INPUT_VCF);
        OutputVcf = cmd.getOptionValue(OUTPUT_VCF);

        RefGenVersion = RefGenomeVersion.from(cmd);
        SV_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);

        ReadFiltering = new ReadFilters(ReadFilterConfig.from(cmd));

        SpecificRegions = Lists.newArrayList();

        try
        {
            SpecificRegions.addAll(loadSpecificRegions(cmd));
        }
        catch(ParseException e)
        {
            System.exit(1);
        }

        LogReadIds = cmd.hasOption(LOG_READ_IDS) ?
                Arrays.stream(cmd.getOptionValue(LOG_READ_IDS).split(ITEM_DELIM, -1)).collect(Collectors.toList()) : Lists.newArrayList();

        Threads = parseThreads(cmd);
    }

    public boolean isValid()
    {
        return checkFileExists(BamFile) && checkFileExists(InputVcf) && checkFileExists(RefGenomeFile);
    }

    private static boolean checkFileExists(final String filename)
    {
        if(!Files.exists(Paths.get(filename)))
        {
            SV_LOGGER.error("missing required file({})", filename);
            return false;
        }

        return true;
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        addOutputOptions(options);
        addLoggingOptions(options);

        options.addOption(SAMPLE, true, "Sample ID to append to existing VCF");
        options.addOption(BAM_FILE, true, "BAM file location");
        options.addOption(INPUT_VCF, true, "Input VCF");
        options.addOption(OUTPUT_VCF, true, "Output VCF");
        addRefGenomeConfig(options);
        addValidationStringencyOption(options);
        ReadFilterConfig.addCmdLineArgs(options);

        addSpecificChromosomesRegionsConfig(options);
        options.addOption(LOG_READ_IDS, true, "Log specific read IDs, separated by ';'");
        addThreadOptions(options);

        return options;
    }

}
