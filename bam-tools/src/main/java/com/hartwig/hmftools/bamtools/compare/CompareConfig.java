package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.ITEM_DELIM;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.LOG_READ_IDS;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.PARTITION_SIZE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.SAMPLE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.addCommonCommandOptions;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.loadSpecificRegionsConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class CompareConfig
{
    public final String OutputFile;
    public final String RefBamFile;
    public final String NewBamFile;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final int PartitionSize;

    public final int Threads;
    public final List<String> LogReadIds;

    // debug
    public final List<String> SpecificChromosomes;
    public final List<ChrBaseRegion> SpecificRegions;

    private static final String OUTPUT_FILE = "output_file";
    private static final String REF_BAM_FILE = "ref_bam";
    private static final String NEW_BAM_FILE = "new_bam";

    private static final int DEFAULT_CHR_PARTITION_SIZE = 100000;

    public CompareConfig(final CommandLine cmd)
    {
        OutputFile = cmd.getOptionValue(OUTPUT_FILE);
        RefBamFile = cmd.getOptionValue(REF_BAM_FILE);
        NewBamFile = cmd.getOptionValue(NEW_BAM_FILE);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME);

        if(RefBamFile == null || NewBamFile == null || OutputFile == null || RefGenomeFile == null)
        {
            BT_LOGGER.error("missing config: bam(ref={} new={}) refGenome({}) outputDir({})",
                    RefBamFile != null, NewBamFile != null, RefGenomeFile != null, OutputFile != null);
            System.exit(1);
        }

        RefGenVersion = RefGenomeVersion.from(cmd);

        BT_LOGGER.info("refGenome({}), bam({})", RefGenVersion, RefBamFile);
        BT_LOGGER.info("output file({})", OutputFile);

        PartitionSize = Integer.parseInt(cmd.getOptionValue(PARTITION_SIZE, String.valueOf(DEFAULT_CHR_PARTITION_SIZE)));

        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();

        if(!loadSpecificRegionsConfig(cmd, SpecificChromosomes, SpecificRegions))
            System.exit(1);

        Threads = parseThreads(cmd);

        LogReadIds = cmd.hasOption(LOG_READ_IDS) ?
                Arrays.stream(cmd.getOptionValue(LOG_READ_IDS).split(ITEM_DELIM, -1)).collect(Collectors.toList()) : Lists.newArrayList();
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();

        addCommonCommandOptions(options);

        options.addOption(OUTPUT_FILE, true, "Output comparison file");
        options.addOption(PARTITION_SIZE, true, "Partition size, default: " + DEFAULT_CHR_PARTITION_SIZE);
        options.addOption(REF_BAM_FILE, true, "Ref BAM file");
        options.addOption(NEW_BAM_FILE, true, "New BAM file");
        options.addOption(LOG_READ_IDS, true, "Log specific read IDs, separated by ';'");

        return options;
    }
}
