package com.hartwig.hmftools.isofox;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.HG19;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.ConfigUtils.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadDelimitedIdFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_FRAG_LENGTH_MIN_COUNT;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_MAX_FRAGMENT_SIZE;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_SINGLE_MAP_QUALITY;
import static com.hartwig.hmftools.isofox.IsofoxFunction.EXPECTED_GC_COUNTS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.EXPECTED_TRANS_COUNTS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.FUSIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.ALT_SPLICE_JUNCTIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.RETAINED_INTRONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.STATISTICS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.expression.ExpectedRatesGenerator.FL_LENGTH;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;
import static com.hartwig.hmftools.isofox.unmapped.UmrCohortFrequency.UMR_COHORT_FREQUENCY_FILE;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.common.GeneRegionFilters;
import com.hartwig.hmftools.isofox.expression.ExpectedRatesGenerator;
import com.hartwig.hmftools.isofox.fusion.FusionConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class IsofoxConfig
{
    // config items
    public static final String SAMPLE = "sample";
    public static final String FUNCTIONS = "functions";

    private static final String CANONICAL_ONLY = "canonical_only";

    private static final String WRITE_EXON_DATA = "write_exon_data";
    private static final String WRITE_READ_DATA = "write_read_data";
    private static final String WRITE_SPLICE_SITE_DATA = "write_splice_sites";

    private static final String BAM_FILE = "bam_file";
    private static final String LONG_FRAGMENT_LIMIT = "long_frag_limit";
    private static final String DROP_DUPLICATES = "drop_dups";
    private static final String MARK_DUPLICATES = "mark_dups";
    private static final String SINGLE_MAP_QUAL = "single_map_qual";
    public static final String GENE_ID_FILE = "gene_id_file";

    private static final String WRITE_FRAG_LENGTHS = "write_frag_lengths";
    private static final String FRAG_LENGTH_MIN_COUNT = "frag_length_min_count";
    private static final String FRAG_LENGTHS_BY_GENE = "frag_length_by_gene";

    private static final String WRITE_GC_DATA = "write_gc_data";

    // expected expression config
    private static final String EXP_COUNTS_FILE = "exp_counts_file";
    private static final String EXP_GC_RATIOS_FILE = "exp_gc_ratios_file";
    private static final String APPLY_EXP_RATES = "apply_exp_rates";
    private static final String READ_LENGTH = "read_length";
    private static final String ER_FRAGMENT_LENGTHS = "exp_rate_frag_lengths";
    private static final String APPLY_GC_BIAS_ADJUSTMENT = "apply_gc_bias_adjust";
    private static final String WRITE_EXPECTED_RATES = "write_exp_rates";
    private static final String WRITE_TRANS_COMBO_DATA = "write_trans_combo_data";
    private static final String APPLY_MQ_ADJUST = "apply_map_qual_adjust";

    // neo-epitopes
    private static final String NEO_EPITOPE_FILE = "neoepitope_file";

    // debug and performance
    private static final String GENE_READ_LIMIT = "gene_read_limit";
    private static final String RUN_VALIDATIONS = "validate";
    private static final String PERF_CHECKS = "run_perf_checks";
    private static final String THREADS = "threads";
    private static final String FILTER_READS_FILE = "filter_reads_file";

    public final String SampleId;

    public final String OutputDir;
    public final String OutputIdentifier; // optionally include extra identifier in output files

    public final List<IsofoxFunction> Functions;

    public final boolean CanonicalTranscriptOnly;
    public final String BamFile;
    public final File RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;
    public final RefGenomeInterface RefGenome;
    public final GeneRegionFilters Filters;
    public final int GeneReadLimit;
    public int MaxFragmentLength;
    public final boolean DropDuplicates;
    public final boolean MarkDuplicates;

    public final boolean WriteExonData;
    public final boolean WriteReadData;
    public final boolean WriteSpliceSiteData;

    public final String ExpCountsFile;
    public final String NeoEpitopeFile;
    public final String UnmappedCohortFreqFile;
    public final String ExpGcRatiosFile;
    public final boolean ApplyFragmentLengthAdjust;
    public int ReadLength;
    public final List<FragmentSize> FragmentSizeData;
    public final boolean WriteExpectedRates;
    public final boolean WriteTransComboData;

    public final boolean WriteFragmentLengths;
    public final int FragmentLengthSamplingCount;
    public final boolean WriteFragmentLengthsByGene;

    public final boolean WriteGcData;

    public final FusionConfig Fusions;

    // debugging and performance options
    public final boolean RunValidations;
    public final boolean RunPerfChecks;
    public final int Threads;
    public final List<String> FilteredReadIds;

    public static final Logger ISF_LOGGER = LogManager.getLogger(IsofoxConfig.class);

    public IsofoxConfig(final CommandLine cmd) throws Exception
    {
        SampleId = cmd.getOptionValue(SAMPLE);

        Functions = Lists.newArrayList();

        if(cmd.hasOption(FUNCTIONS))
        {
            final String[] functionsStr = cmd.getOptionValue(FUNCTIONS).split(ITEM_DELIM);

            for(final String functionStr : functionsStr)
            {
                Functions.add(IsofoxFunction.valueOf(functionStr));
                ISF_LOGGER.info("running function(s): {}", functionStr);
            }
        }
        else
        {
            Functions.add(TRANSCRIPT_COUNTS);
            Functions.add(ALT_SPLICE_JUNCTIONS);
            Functions.add(FUSIONS);
        }

        CanonicalTranscriptOnly = cmd.hasOption(CANONICAL_ONLY);

        OutputDir = parseOutputDir(cmd);
        OutputIdentifier = cmd.getOptionValue(OUTPUT_ID);

        BamFile = cmd.getOptionValue(BAM_FILE);

        final String refGenomeFilename = cmd.getOptionValue(REF_GENOME);
        RefGenomeFile = refGenomeFilename != null ? new File(refGenomeFilename) : null;
        RefGenome = loadRefGenome(refGenomeFilename);

        if(cmd.hasOption(REF_GENOME_VERSION))
        {
            RefGenVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION));
        }
        else
        {
            RefGenomeVersion refGenVersionOverride = checkRefGenomeVersion();
            RefGenVersion = refGenVersionOverride != null ? refGenVersionOverride : V37;
        }

        Filters = new GeneRegionFilters(RefGenVersion);
        Filters.loadConfig(cmd);

        GeneReadLimit = Integer.parseInt(cmd.getOptionValue(GENE_READ_LIMIT, "0"));

        MaxFragmentLength = Integer.parseInt(cmd.getOptionValue(LONG_FRAGMENT_LIMIT, String.valueOf(DEFAULT_MAX_FRAGMENT_SIZE)));
        IsofoxConstants.SINGLE_MAP_QUALITY = Short.parseShort(cmd.getOptionValue(SINGLE_MAP_QUAL, String.valueOf(DEFAULT_SINGLE_MAP_QUALITY)));
        DropDuplicates = cmd.hasOption(DROP_DUPLICATES);
        MarkDuplicates = cmd.hasOption(MARK_DUPLICATES);

        WriteExonData = cmd.hasOption(WRITE_EXON_DATA);
        WriteFragmentLengths = cmd.hasOption(WRITE_FRAG_LENGTHS);
        WriteFragmentLengthsByGene = cmd.hasOption(FRAG_LENGTHS_BY_GENE);
        WriteReadData = cmd.hasOption(WRITE_READ_DATA);
        WriteSpliceSiteData = cmd.hasOption(WRITE_SPLICE_SITE_DATA);
        WriteTransComboData = cmd.hasOption(WRITE_TRANS_COMBO_DATA);
        WriteGcData = cmd.hasOption(WRITE_GC_DATA);

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));

        if(Functions.contains(TRANSCRIPT_COUNTS))
        {
            ExpCountsFile = cmd.getOptionValue(EXP_COUNTS_FILE);
            ExpGcRatiosFile = cmd.getOptionValue(EXP_GC_RATIOS_FILE);
        }
        else
        {
            ExpCountsFile = null;
            ExpGcRatiosFile = null;
        }

        NeoEpitopeFile = cmd.getOptionValue(NEO_EPITOPE_FILE);
        UnmappedCohortFreqFile = cmd.getOptionValue(UMR_COHORT_FREQUENCY_FILE);

        WriteExpectedRates = cmd.hasOption(WRITE_EXPECTED_RATES);

        ApplyFragmentLengthAdjust = ExpCountsFile != null;
        int defaultFragLengthSamplingCount = ApplyFragmentLengthAdjust ? DEFAULT_FRAG_LENGTH_MIN_COUNT : 0;
        FragmentLengthSamplingCount = Integer.parseInt(cmd.getOptionValue(FRAG_LENGTH_MIN_COUNT, String.valueOf(defaultFragLengthSamplingCount)));

        ReadLength = Integer.parseInt(cmd.getOptionValue(READ_LENGTH, "0"));
        FragmentSizeData = Lists.newArrayList();

        if(cmd.hasOption(ER_FRAGMENT_LENGTHS))
        {
            String[] fragLengths = cmd.getOptionValue(ER_FRAGMENT_LENGTHS).split(ITEM_DELIM);
            for(int i = 0; i < fragLengths.length; ++i)
            {
                String[] flItem = fragLengths[i].split("-");
                int fragLength = Integer.parseInt(flItem[FL_LENGTH]);
                int fragFrequency = max(Integer.parseInt(flItem[ExpectedRatesGenerator.FL_FREQUENCY]), 1);
                FragmentSizeData.add(new FragmentSize(fragLength, fragFrequency));
            }
        }

        Fusions = Functions.contains(FUSIONS) ? new FusionConfig(cmd) : new FusionConfig();

        RunValidations = cmd.hasOption(RUN_VALIDATIONS);
        RunPerfChecks = cmd.hasOption(PERF_CHECKS);

        FilteredReadIds = cmd.hasOption(FILTER_READS_FILE) ?
                loadDelimitedIdFile(cmd.getOptionValue(FILTER_READS_FILE), "FilteredReadIds", CSV_DELIM) : null;
    }

    public boolean isValid()
    {
        if(OutputDir == null)
        {
            ISF_LOGGER.error("no output directory specified");
            return false;
        }

        if(!FragmentSizeData.isEmpty())
        {
            if(FragmentSizeData.stream().anyMatch(x -> x.Length <= 0 || x.Frequency < 0))
            {
                ISF_LOGGER.error("invalid fragment lengths");
                return false;
            }
        }

        if(runFunction(EXPECTED_TRANS_COUNTS))
        {
            if(ReadLength == 0 || FragmentSizeData.isEmpty())
            {
                ISF_LOGGER.error("invalid read or fragment lengths for generating expected trans counts");
                return false;
            }

            return true;
        }

        if(runFunction(EXPECTED_GC_COUNTS))
        {
            if(ReadLength == 0 || RefGenomeFile == null)
            {
                ISF_LOGGER.error("invalid read length or ref genome for generating expected GC ratio counts");
                return false;
            }

            return true;
        }

        if(BamFile == null || BamFile.isEmpty() || !Files.exists(Paths.get(BamFile)))
        {
            if(!runFusionsOnly() || Fusions.ChimericReadsFile == null)
            {
                ISF_LOGGER.error("BAM file({}) missing or not found", BamFile);
                return false;
            }
        }

        if(SampleId == null || SampleId.isEmpty())
        {
            ISF_LOGGER.error("sampleId missing");
            return false;
        }

        if(RefGenomeFile == null)
        {
            ISF_LOGGER.error("ref genome missing");
            return false;
        }

        if(Filters.RestrictedGeneIds.isEmpty() && (WriteExonData || WriteReadData))
        {
            ISF_LOGGER.warn("writing exon and/or read data for all transcripts may be slow and generate large output files");
        }

        return true;
    }

    public static boolean validConfigPaths(final CommandLine cmd)
    {
        return configPathValid(cmd, OUTPUT_DIR) && configPathValid(cmd, REF_GENOME)  && configPathValid(cmd, ENSEMBL_DATA_DIR)
                && configPathValid(cmd, GENE_ID_FILE) && configPathValid(cmd, BAM_FILE) && configPathValid(cmd, EXP_COUNTS_FILE)
                && configPathValid(cmd, EXP_GC_RATIOS_FILE) && configPathValid(cmd, NEO_EPITOPE_FILE);
    }

    public static boolean configPathValid(final CommandLine cmd, final String configItem)
    {
        if(!cmd.hasOption(configItem))
            return true;

        final String filePath = cmd.getOptionValue(configItem);
        if(!Files.exists(Paths.get(filePath)))
        {
            ISF_LOGGER.error("invalid config path: {} = {}", configItem, filePath);
            return false;
        }

        return true;
    }

    public boolean runFunction(IsofoxFunction function) { return Functions.contains(function); }
    public boolean runFusionsOnly() { return Functions.contains(FUSIONS) && Functions.size() == 1; }
    public boolean runStatisticsOnly() { return Functions.contains(STATISTICS) && Functions.size() == 1; }

    public boolean generateExpectedDataOnly()
    {
        return runFunction(EXPECTED_GC_COUNTS) || runFunction(EXPECTED_TRANS_COUNTS);
    }

    public boolean requireFragmentLengthCalcs()
    {
        return WriteFragmentLengths || FragmentLengthSamplingCount > 0;
    }

    public boolean requireGcRatioCalcs() { return WriteGcData || ExpGcRatiosFile != null; }
    public boolean applyGcBiasAdjust() { return ExpGcRatiosFile != null; }
    public boolean applyExpectedCounts() { return ExpCountsFile != null; }

    public String formOutputFile(final String fileId)
    {
        if(OutputIdentifier != null)
            return OutputDir + SampleId + ISF_FILE_ID + OutputIdentifier + "." + fileId;
        else
            return OutputDir + SampleId + ISF_FILE_ID + fileId;
    }

    private RefGenomeVersion checkRefGenomeVersion()
    {
        if(BamFile == null || !Files.exists(Paths.get(BamFile)) || RefGenomeFile == null)
            return null;

        final SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(RefGenomeFile).open(new File(BamFile));

        if(samReader == null)
            return null;

        if(RefGenomeFunctions.samReaderUsesChrInContigs(samReader))
            return HG19;

        return V37;
    }

    public boolean skipFilteredRead(final String readId) { return FilteredReadIds != null && !FilteredReadIds.contains(readId); }

    public IsofoxConfig()
    {
        SampleId = "TEST";

        Functions = Lists.newArrayList();
        Functions.add(TRANSCRIPT_COUNTS);
        Functions.add(ALT_SPLICE_JUNCTIONS);
        Functions.add(RETAINED_INTRONS);

        Filters = new GeneRegionFilters(V37);
        OutputDir = null;
        BamFile = null;
        RefGenomeFile = null;
        RefGenVersion = V37;
        RefGenome = new MockRefGenome();
        CanonicalTranscriptOnly = false;
        GeneReadLimit = 0;
        MaxFragmentLength = DEFAULT_MAX_FRAGMENT_SIZE;
        DropDuplicates = false;
        MarkDuplicates = false;

        ReadLength = 0;
        FragmentSizeData = Lists.newArrayList();
        ExpCountsFile = null;
        ExpGcRatiosFile = null;
        NeoEpitopeFile = null;
        UnmappedCohortFreqFile = null;

        WriteExonData = false;
        WriteReadData = false;
        WriteSpliceSiteData = false;
        WriteFragmentLengths = false;
        WriteTransComboData = false;
        WriteGcData = false;
        Fusions = new FusionConfig();

        WriteExpectedRates = false;
        ApplyFragmentLengthAdjust = false;
        OutputIdentifier = null;
        WriteFragmentLengthsByGene = false;
        FragmentLengthSamplingCount = 0;

        RunValidations = true;
        RunPerfChecks = false;
        Threads = 0;
        FilteredReadIds = null;
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample ID");
        addEnsemblDir(options);
        addOutputDir(options);
        addLoggingOptions(options);

        options.addOption(FUNCTIONS, true, "Optional: list of functional routines to run (see documentation)");
        options.addOption(CANONICAL_ONLY, false, "Check all transcripts, not just canonical");
        options.addOption(GENE_READ_LIMIT, true, "Per-gene limit on max reads processed (default=0, not applied)");
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version - accepts 37 (default) or 38");
        options.addOption(LONG_FRAGMENT_LIMIT, true, "Max RNA fragment size");
        options.addOption(DROP_DUPLICATES, false, "Include duplicate fragments in expression calculations");
        options.addOption(MARK_DUPLICATES, false, "Manually identify duplicate fragments");
        options.addOption(FRAG_LENGTH_MIN_COUNT, true, "Fragment length measurement - min read fragments required");
        options.addOption(FRAG_LENGTHS_BY_GENE, false, "Write fragment lengths by gene");
        options.addOption(BAM_FILE, true, "RNA BAM file location");
        options.addOption(GENE_ID_FILE, true, "Optional CSV file of genes to analyse");
        options.addOption(WRITE_EXON_DATA, false, "Exon region data");
        options.addOption(WRITE_READ_DATA, false, "BAM read data");
        options.addOption(WRITE_SPLICE_SITE_DATA, false, "Write support info for each splice site");
        options.addOption(WRITE_TRANS_COMBO_DATA, false, "Write transcript group data for EM algo");
        options.addOption(WRITE_FRAG_LENGTHS, false, "Write intronic fragment lengths to log");

        options.addOption(WRITE_GC_DATA, false, "Write GC ratio counts from all genic reads");

        options.addOption(APPLY_EXP_RATES, false, "Generate expected expression rates for transcripts");
        options.addOption(EXP_COUNTS_FILE, true, "File with generated expected expression rates per transcript");
        options.addOption(EXP_GC_RATIOS_FILE, true, "File with generated expected GC ratios per transcript");
        options.addOption(NEO_EPITOPE_FILE, true, "File with neo-epitopes to measure fragment support");
        options.addOption(UMR_COHORT_FREQUENCY_FILE, true, "Unmapped reads cohort frequency file");
        options.addOption(READ_LENGTH, true, "Sample sequencing read length (eg 76 or 151 bases");
        options.addOption(SINGLE_MAP_QUAL, true, "Optional - map quality for reads mapped to a single location (default=255)");
        options.addOption(APPLY_MQ_ADJUST, false, "Use read map quality to adjust expected counts");
        options.addOption(APPLY_GC_BIAS_ADJUSTMENT, false, "Use GC Bias adjustments in expected rate calcs");

        options.addOption(ER_FRAGMENT_LENGTHS, true,
                "Fragment sizes and weights for expected transcript calcs (format: length1-freq1;length3-freq2 eg 100-10;150-20) in integer terms");

        options.addOption(WRITE_EXPECTED_RATES, false, "Write sample expected expression rates to file");

        options.addOption(OUTPUT_ID, true, "Optionally add identifier to output files");
        options.addOption(THREADS, true, "Number of threads to use (default=0, single-threaded)");
        options.addOption(RUN_VALIDATIONS, false, "Run auto-validations");
        options.addOption(PERF_CHECKS, false, "Run performance logging routines");
        options.addOption(FILTER_READS_FILE, true, "Only process reads in this file");

        GeneRegionFilters.addCommandLineOptions(options);
        FusionConfig.addCommandLineOptions(options);

        return options;
    }
}
