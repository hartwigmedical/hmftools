package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SUB_ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificChromsomesOrRegions;
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
import com.hartwig.hmftools.svprep.reads.ReadFilterConfig;
import com.hartwig.hmftools.svprep.reads.ReadFilters;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class SvConfig
{
    public final String SampleId;
    public final String BamFile;
    public final String RefGenomeFile;
    // public final RefGenomeSource RefGenome;
    public final RefGenomeVersion RefGenVersion;

    public final ReadFilters ReadFiltering;
    public final HotspotCache Hotspots;

    public final int PartitionSize;
    public final int BucketSize;
    public final int ReadLength;
    public final boolean CalcFragmentLength;

    public final String OutputDir;
    public final String OutputId;
    public final Set<WriteType> WriteTypes;

    public final int Threads;

    // debug
    public final List<String> SpecificChromosomes;
    public final List<String> LogReadIds;
    public final int MaxPartitionReads;
    public final List<ChrBaseRegion> SpecificRegions;

    private boolean mIsValid;

    // config strings
    public static final String SAMPLE = "sample";
    private static final String BAM_FILE = "bam_file";
    private static final String KNOWN_FUSION_BED = "known_fusion_bed";

    private static final String WRITE_TYPES = "write_types";

    private static final String THREADS = "threads";
    private static final String READ_LENGTH = "read_length";
    private static final String FRAG_LENGTH_RANGE = "fragment_length_range";
    private static final String CALC_FRAG_LENGTH = "calc_fragment_length";

    private static final String LOG_READ_IDS = "log_read_ids";
    private static final String MAX_PARTITION_READS = "max_partition_reads";

    public SvConfig(final CommandLine cmd)
    {
        mIsValid = true;

        SampleId = cmd.getOptionValue(SAMPLE);
        BamFile = cmd.getOptionValue(BAM_FILE);
        RefGenomeFile = cmd.getOptionValue(REF_GENOME);
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);

        if(SampleId == null || BamFile == null || OutputDir == null || RefGenomeFile == null)
        {
            SV_LOGGER.error("missing config: sample({}) bam({}) refGenome({}) outputDir({})",
                    SampleId != null, BamFile != null, RefGenomeFile != null, OutputDir != null);
            mIsValid = false;
        }

        RefGenVersion = cmd.hasOption(REF_GENOME_VERSION) ? RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION)) : V37;

        Hotspots = new HotspotCache(cmd.getOptionValue(KNOWN_FUSION_BED));

        PartitionSize = DEFAULT_CHR_PARTITION_SIZE;
        BucketSize = DEFAULT_BUCKET_SIZE;

        ReadLength = Integer.parseInt(cmd.getOptionValue(READ_LENGTH, String.valueOf(DEFAULT_READ_LENGTH)));

        ReadFiltering = new ReadFilters(ReadFilterConfig.from(cmd));

        CalcFragmentLength = cmd.hasOption(CALC_FRAG_LENGTH);

        if(!CalcFragmentLength && cmd.hasOption(FRAG_LENGTH_RANGE))
        {
            String[] fragRange = cmd.getOptionValue(FRAG_LENGTH_RANGE).split(SUB_ITEM_DELIM, 2);
            ReadFiltering.config().setFragmentLengthMin(Integer.parseInt(fragRange[0]), Integer.parseInt(fragRange[1]));
        }

        WriteTypes = Sets.newHashSet();

        if(cmd.hasOption(WRITE_TYPES))
        {
            String[] writeTypes = cmd.getOptionValue(WRITE_TYPES).split(ITEM_DELIM, -1);
            Arrays.stream(writeTypes).forEach(x -> WriteTypes.add(WriteType.valueOf(x)));
        }
        else
        {
            WriteTypes.add(WriteType.JUNCTIONS);
        }

        SpecificChromosomes = Lists.newArrayList();
        SpecificRegions = Lists.newArrayList();

        try
        {
            loadSpecificChromsomesOrRegions(cmd, SpecificChromosomes, SpecificRegions, SV_LOGGER);
        }
        catch(ParseException e)
        {
            mIsValid = false;
        }

        LogReadIds = cmd.hasOption(LOG_READ_IDS) ?
                Arrays.stream(cmd.getOptionValue(LOG_READ_IDS).split(ITEM_DELIM, -1)).collect(Collectors.toList()) : Lists.newArrayList();

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));
        MaxPartitionReads = Integer.parseInt(cmd.getOptionValue(MAX_PARTITION_READS, "0"));
    }

    public boolean isValid() { return mIsValid; }

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
            case BAM: return filename + "bam";
            case JUNCTIONS: return filename + "junctions.csv";
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
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version - accepts 37 (default) or 38");
        options.addOption(KNOWN_FUSION_BED, true, "Known fusion hotspot BED file");
        options.addOption(READ_LENGTH, true, "Read length");
        options.addOption(CALC_FRAG_LENGTH, false, "Calculate distribution for fragment length");
        options.addOption(FRAG_LENGTH_RANGE, true, "Empirical fragment length range: Min:Max");
        options.addOption(WRITE_TYPES, true, "Write types: " + WriteType.values().toString());
        addSpecificChromosomesRegionsConfig(options);
        options.addOption(LOG_READ_IDS, true, "Log specific read IDs, separated by ';'");
        options.addOption(MAX_PARTITION_READS, true, "Limit to stop processing reads in partition, for debug");
        options.addOption(THREADS, true, "Thread count");
        ReadFilterConfig.addCmdLineArgs(options);

        return options;
    }
}
