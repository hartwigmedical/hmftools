package com.hartwig.hmftools.svtools.sv_prep;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS_DESC;
import static com.hartwig.hmftools.svtools.sv_prep.SvCommon.ITEM_DELIM;
import static com.hartwig.hmftools.svtools.sv_prep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.DEFAULT_BUCKET_SIZE;
import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.DEFAULT_CHR_PARTITION_SIZE;

import java.util.Arrays;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.util.Strings;

public class SvConfig
{
    public final String SampleId;
    public final String BamFile;
    public final String RefGenomeFile;
    // public final RefGenomeSource RefGenome;
    public final RefGenomeVersion RefGenVersion;

    public final ReadFilterConfig ReadFilters;

    public final Set<String> SpecificChromosomes;
    public final List<ChrBaseRegion> SpecificRegions;

    public final String OutputDir;
    public final String OutputId;

    public final int PartitionSize;
    public final int BucketSize;
    public final int Threads;

    public final Set<WriteType> WriteTypes;

    private boolean mIsValid;

    // config strings
    public static final String SAMPLE = "sample";
    private static final String BAM_FILE = "bam_file";

    private static final String SPECIFIC_CHROMOSOMES = "specific_chr";
    private static final String WRITE_TYPES = "write_types";

    private static final String THREADS = "threads";

    public SvConfig(final CommandLine cmd)
    {
        mIsValid = true;

        SampleId = cmd.getOptionValue(SAMPLE);
        BamFile = cmd.getOptionValue(BAM_FILE);
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);

        RefGenomeFile = cmd.getOptionValue(REF_GENOME);
        // RefGenome = loadRefGenome(RefGenomeFile);

        RefGenVersion = cmd.hasOption(REF_GENOME_VERSION) ? RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION)) : V37;

        PartitionSize = DEFAULT_CHR_PARTITION_SIZE;
        BucketSize = DEFAULT_BUCKET_SIZE;

        ReadFilters = new ReadFilterConfig();

        WriteTypes = Sets.newHashSet();

        String[] writeTypes = cmd.getOptionValue(WRITE_TYPES).split(ITEM_DELIM, -1);
        Arrays.stream(writeTypes).forEach(x -> WriteTypes.add(WriteType.valueOf(x)));

        SpecificChromosomes = Sets.newHashSet();
        SpecificRegions = Lists.newArrayList();

        if(cmd.hasOption(SPECIFIC_REGIONS))
        {
            try
            {
                SpecificRegions.addAll(ChrBaseRegion.loadSpecificRegions(cmd));

                for(ChrBaseRegion region : SpecificRegions)
                {
                    SV_LOGGER.info("filtering for specific region: {}", region);
                    SpecificChromosomes.add(region.Chromosome);
                }
            }
            catch(Exception e)
            {
                SV_LOGGER.error("invalid specific regions: {}", cmd.getOptionValue(SPECIFIC_REGIONS));
                mIsValid = false;
            }
        }
        else if(cmd.hasOption(SPECIFIC_CHROMOSOMES))
        {
            final String chromosomeList = cmd.getOptionValue(SPECIFIC_CHROMOSOMES, Strings.EMPTY);
            if(!chromosomeList.isEmpty())
            {
                SpecificChromosomes.addAll(Lists.newArrayList(chromosomeList.split(ITEM_DELIM)));
            }
        }

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        addOutputDir(options);
        addLoggingOptions(options);

        options.addOption(SAMPLE, true, "Tumor sample ID");
        options.addOption(BAM_FILE, true, "RNA BAM file location");

        // reference data
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version - accepts 37 (default) or 38");

        options.addOption(WRITE_TYPES, true, "Write types: " + WriteType.values().toString());
        options.addOption(SPECIFIC_CHROMOSOMES, true, "Specific chromosomes separated by ';'");
        options.addOption(SPECIFIC_REGIONS, true, SPECIFIC_REGIONS_DESC);
        options.addOption(THREADS, true, "Thread count");

        return options;
    }
}
