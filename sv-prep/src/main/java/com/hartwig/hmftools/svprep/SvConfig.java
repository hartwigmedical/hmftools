package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS_DESC;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SUB_ITEM_DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.ITEM_DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConstants.DEFAULT_BUCKET_SIZE;
import static com.hartwig.hmftools.svprep.SvConstants.DEFAULT_CHR_PARTITION_SIZE;
import static com.hartwig.hmftools.svprep.SvConstants.DEFAULT_READ_LENGTH;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

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

    public final int PartitionSize;
    public final int BucketSize;
    public final int ReadLength;
    public final boolean CalcFragmentLength;

    public final String OutputDir;
    public final String OutputId;

    public final int Threads;
    public final Set<String> SpecificChromosomes;
    public final List<String> LogReadIds;
    public final List<ChrBaseRegion> SpecificRegions;

    public final Set<WriteType> WriteTypes;

    private boolean mIsValid;

    // config strings
    public static final String SAMPLE = "sample";
    private static final String BAM_FILE = "bam_file";

    public static final String SPECIFIC_CHROMOSOMES = "specific_chr";
    private static final String LOG_READ_IDS = "log_read_ids";
    private static final String WRITE_TYPES = "write_types";

    private static final String THREADS = "threads";
    private static final String READ_LENGTH = "read_length";
    private static final String FRAG_LENGTH_RANGE = "fragment_length_range";
    private static final String CALC_FRAG_LENGTH = "calc_fragment_length";

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

        ReadLength = Integer.parseInt(cmd.getOptionValue(READ_LENGTH, String.valueOf(DEFAULT_READ_LENGTH)));

        ReadFilters = new ReadFilterConfig();

        CalcFragmentLength = cmd.hasOption(CALC_FRAG_LENGTH);

        if(!CalcFragmentLength && cmd.hasOption(FRAG_LENGTH_RANGE))
        {
            String[] fragRange = cmd.getOptionValue(FRAG_LENGTH_RANGE).split(SUB_ITEM_DELIM, 2);
            ReadFilters.setFragmentLengthMin(Integer.parseInt(fragRange[0]), Integer.parseInt(fragRange[1]));
        }

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

        LogReadIds = cmd.hasOption(LOG_READ_IDS) ?
                Arrays.stream(cmd.getOptionValue(LOG_READ_IDS).split(ITEM_DELIM, -1)).collect(Collectors.toList()) : Lists.newArrayList();

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));
    }

    public String formFilename(final WriteType writeType)
    {
        String filename = OutputDir + SampleId;

        filename += ".sv_prep.";

        if(OutputId != null)
            filename += OutputId + ".";

        switch(writeType)
        {
            case SV_BED: return filename + "bed";
            case READS: return filename + "reads.csv";
            case BUCKET_STATS: return filename + "buckets.csv";
            case BAM: return filename + "bam";
            case FRAGMENT_LENGTH_DIST: return filename + "fragment_lengths";
        }

        return null;
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        addOutputOptions(options);
        addLoggingOptions(options);

        options.addOption(SAMPLE, true, "Tumor sample ID");
        options.addOption(BAM_FILE, true, "RNA BAM file location");
        options.addOption(READ_LENGTH, true, "Read length");
        options.addOption(CALC_FRAG_LENGTH, false, "Calculate distribution for fragment length");
        options.addOption(FRAG_LENGTH_RANGE, true, "Empirical fragment length range: Min:Max");

        // reference data
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version - accepts 37 (default) or 38");

        options.addOption(WRITE_TYPES, true, "Write types: " + WriteType.values().toString());
        options.addOption(SPECIFIC_CHROMOSOMES, true, "Specific chromosomes separated by ';'");
        options.addOption(SPECIFIC_REGIONS, true, SPECIFIC_REGIONS_DESC);
        options.addOption(LOG_READ_IDS, true, "Log specific read IDs, separated by ';'");
        options.addOption(THREADS, true, "Thread count");

        return options;
    }
}
