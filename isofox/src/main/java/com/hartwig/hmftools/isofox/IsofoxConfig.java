package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.NEO_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.NEO_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadDelimitedIdFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_EXPECTED_RATE_LENGTHS;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_FRAG_LENGTH_MIN_COUNT;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_MAX_FRAGMENT_SIZE;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_SINGLE_MAP_QUALITY;
import static com.hartwig.hmftools.isofox.IsofoxFunction.FUSIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.ALT_SPLICE_JUNCTIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.RETAINED_INTRONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.STATISTICS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.expression.ExpectedRatesCommon.ER_FRAGMENT_LENGTHS;
import static com.hartwig.hmftools.isofox.expression.ExpectedRatesCommon.ER_FRAGMENT_LENGTHS_DESC;
import static com.hartwig.hmftools.isofox.expression.ExpectedRatesCommon.loadFragmentSizeConfig;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.common.GeneRegionFilters;
import com.hartwig.hmftools.isofox.fusion.FusionConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class IsofoxConfig
{
    // config items
    public static final String FUNCTIONS = "functions";

    private static final String CANONICAL_ONLY = "canonical_only";

    private static final String BAM_FILE = "bam_file";
    public static final String LONG_FRAGMENT_LIMIT = "long_frag_limit";
    public static final String READ_LENGTH = "read_length";

    private static final String WRITE_EXON_DATA = "write_exon_data";
    private static final String WRITE_READ_DATA = "write_read_data";
    private static final String WRITE_SPLICE_SITE_DATA = "write_splice_sites";
    private static final String WRITE_SPLICE_JUNC_DATA = "write_splice_junctions";
    private static final String WRITE_FRAG_LENGTHS = "write_frag_lengths";
    private static final String FRAG_LENGTH_MIN_COUNT = "frag_length_min_count";
    private static final String FRAG_LENGTHS_BY_GENE = "frag_length_by_gene";
    private static final String WRITE_GC_DATA = "write_gc_data";
    private static final String WRITE_TRANS_COMBO_DATA = "write_trans_combo_data";

    // expected expression config
    private static final String EXP_COUNTS_FILE = "exp_counts_file";
    private static final String EXP_GC_RATIOS_FILE = "exp_gc_ratios_file";
    private static final String PANEL_TPM_NORM_FILE = "panel_tpm_norm_file";

    private static final String DROP_DUPLICATES = "drop_dups";
    private static final String SINGLE_MAP_QUAL = "single_map_qual";

    // debug and performance
    private static final String GENE_READ_LIMIT = "gene_read_limit";
    private static final String RUN_VALIDATIONS = "validate";
    private static final String PERF_CHECKS = "run_perf_checks";
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
    public int ReadLength;
    public int MaxFragmentLength;
    public final boolean DropDuplicates;

    // expression
    public final String ExpCountsFile;
    public final String ExpGcRatiosFile;
    public final String PanelTpmNormFile;
    public final String NeoDir;
    public final boolean ApplyFragmentLengthAdjust;
    public final List<FragmentSize> FragmentSizeData;

    public final boolean WriteExonData;
    public final boolean WriteSpliceJunctions;
    public final boolean WriteReadData;
    public final boolean WriteSpliceSiteData;
    public final boolean WriteTransComboData;
    public final boolean WriteFragmentLengths;
    public final int FragmentLengthSamplingCount;
    public final boolean WriteFragmentLengthsByGene;
    public final boolean WriteGcData;

    public final FusionConfig Fusions;

    // debugging and performance options
    public final int GeneReadLimit;
    public final boolean RunValidations;
    public final boolean RunPerfChecks;
    public final int Threads;
    public final List<String> FilteredReadIds;

    public static final Logger ISF_LOGGER = LogManager.getLogger(IsofoxConfig.class);

    public IsofoxConfig(final ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(SAMPLE);

        Functions = Lists.newArrayList();

        if(configBuilder.hasValue(FUNCTIONS))
        {
            final String[] functionsStr = configBuilder.getValue(FUNCTIONS).split(ITEM_DELIM);

            for(final String functionStr : functionsStr)
            {
                Functions.add(IsofoxFunction.valueOf(functionStr));
            }
        }
        else
        {
            Functions.add(TRANSCRIPT_COUNTS);
            Functions.add(ALT_SPLICE_JUNCTIONS);
            Functions.add(FUSIONS);
        }

        ISF_LOGGER.info("running function(s): {}",
                Functions.stream().map(x -> x.toString()).collect(Collectors.joining(",")));

        CanonicalTranscriptOnly = configBuilder.hasValue(CANONICAL_ONLY);

        OutputDir = parseOutputDir(configBuilder);
        OutputIdentifier = configBuilder.getValue(OUTPUT_ID);

        BamFile = configBuilder.getValue(BAM_FILE);

        final String refGenomeFilename = configBuilder.getValue(REF_GENOME);
        RefGenomeFile = refGenomeFilename != null ? new File(refGenomeFilename) : null;
        RefGenome = loadRefGenome(refGenomeFilename);

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        Filters = new GeneRegionFilters(RefGenVersion);
        Filters.loadConfig(configBuilder);

        GeneReadLimit = Integer.parseInt(configBuilder.getValue(GENE_READ_LIMIT, "0"));

        MaxFragmentLength = configBuilder.getInteger(LONG_FRAGMENT_LIMIT);
        IsofoxConstants.SINGLE_MAP_QUALITY = (short)configBuilder.getInteger(SINGLE_MAP_QUAL);
        DropDuplicates = configBuilder.hasValue(DROP_DUPLICATES);

        WriteExonData = configBuilder.hasFlag(WRITE_EXON_DATA);
        WriteSpliceJunctions = configBuilder.hasFlag(WRITE_SPLICE_JUNC_DATA);
        WriteFragmentLengths = configBuilder.hasFlag(WRITE_FRAG_LENGTHS);
        WriteFragmentLengthsByGene = configBuilder.hasFlag(FRAG_LENGTHS_BY_GENE);
        WriteReadData = configBuilder.hasFlag(WRITE_READ_DATA);
        WriteSpliceSiteData = configBuilder.hasFlag(WRITE_SPLICE_SITE_DATA);
        WriteTransComboData = configBuilder.hasFlag(WRITE_TRANS_COMBO_DATA);
        WriteGcData = configBuilder.hasFlag(WRITE_GC_DATA);

        Threads = parseThreads(configBuilder);

        if(Functions.contains(TRANSCRIPT_COUNTS))
        {
            ExpCountsFile = configBuilder.getValue(EXP_COUNTS_FILE);
            ExpGcRatiosFile = configBuilder.getValue(EXP_GC_RATIOS_FILE);
        }
        else
        {
            ExpCountsFile = null;
            ExpGcRatiosFile = null;
        }

        NeoDir = configBuilder.getValue(NEO_DIR_CFG);
        PanelTpmNormFile = configBuilder.getValue(PANEL_TPM_NORM_FILE);

        ApplyFragmentLengthAdjust = ExpCountsFile != null;

        int defaultFragLengthSamplingCount = ApplyFragmentLengthAdjust ? DEFAULT_FRAG_LENGTH_MIN_COUNT : 0;
        FragmentLengthSamplingCount = configBuilder.hasValue(FRAG_LENGTH_MIN_COUNT) ?
                configBuilder.getInteger(FRAG_LENGTH_MIN_COUNT) : defaultFragLengthSamplingCount;

        ReadLength = configBuilder.getInteger(READ_LENGTH);
        FragmentSizeData = loadFragmentSizeConfig(configBuilder);

        Fusions = Functions.contains(FUSIONS) ? new FusionConfig(configBuilder) : new FusionConfig();

        RunValidations = configBuilder.hasValue(RUN_VALIDATIONS);
        RunPerfChecks = configBuilder.hasValue(PERF_CHECKS);

        FilteredReadIds = configBuilder.hasValue(FILTER_READS_FILE) ?
                loadDelimitedIdFile(configBuilder.getValue(FILTER_READS_FILE), "FilteredReadIds", CSV_DELIM) : null;
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

        if(BamFile == null || BamFile.isEmpty() || !Files.exists(Paths.get(BamFile)))
        {
            if(!runFusionsOnly())
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

    public boolean runFunction(IsofoxFunction function) { return Functions.contains(function); }
    public boolean runFusionsOnly() { return Functions.contains(FUSIONS) && Functions.size() == 1; }
    public boolean runStatisticsOnly() { return Functions.contains(STATISTICS) && Functions.size() == 1; }

    public boolean requireFragmentLengthCalcs()
    {
        return WriteFragmentLengths || FragmentLengthSamplingCount > 0;
    }

    public boolean requireGcRatioCalcs() { return WriteGcData || ExpGcRatiosFile != null; }
    public boolean applyGcBiasAdjust() { return ExpGcRatiosFile != null; }

    public String formOutputFile(final String fileId)
    {
        if(OutputIdentifier != null)
            return OutputDir + SampleId + ISF_FILE_ID + OutputIdentifier + "." + fileId;
        else
            return OutputDir + SampleId + ISF_FILE_ID + fileId;
    }

    public boolean skipFilteredRead(final String readId) { return FilteredReadIds != null && !FilteredReadIds.contains(readId); }

    public IsofoxConfig(final RefGenomeInterface refGenome)
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
        RefGenome = refGenome;
        CanonicalTranscriptOnly = false;
        GeneReadLimit = 0;
        MaxFragmentLength = DEFAULT_MAX_FRAGMENT_SIZE;
        DropDuplicates = false;

        ReadLength = 0;
        FragmentSizeData = Lists.newArrayList();
        ExpCountsFile = null;
        ExpGcRatiosFile = null;
        NeoDir = null;
        PanelTpmNormFile = null;

        WriteExonData = false;
        WriteSpliceJunctions = false;
        WriteReadData = false;
        WriteSpliceSiteData = false;
        WriteFragmentLengths = false;
        WriteTransComboData = false;
        WriteGcData = false;
        Fusions = new FusionConfig();

        ApplyFragmentLengthAdjust = false;
        OutputIdentifier = null;
        WriteFragmentLengthsByGene = false;
        FragmentLengthSamplingCount = 0;

        RunValidations = true;
        RunPerfChecks = false;
        Threads = 0;
        FilteredReadIds = null;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        addEnsemblDir(configBuilder);
        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.addConfigItem(FUNCTIONS, false, "List of functional routines to run (see documentation)");
        configBuilder.addFlag(CANONICAL_ONLY, "Check all transcripts, not just canonical");
        configBuilder.addInteger(GENE_READ_LIMIT, "Per-gene limit on max reads processed (0 = not applied)", 0);
        addRefGenomeConfig(configBuilder, true);
        configBuilder.addInteger(LONG_FRAGMENT_LIMIT, "Max RNA fragment size", DEFAULT_MAX_FRAGMENT_SIZE);
        configBuilder.addFlag(DROP_DUPLICATES, "Include duplicate fragments in expression calculations");

        configBuilder.addInteger(
                FRAG_LENGTH_MIN_COUNT, "Fragment length measurement - min read fragments required", DEFAULT_FRAG_LENGTH_MIN_COUNT);

        configBuilder.addFlag(FRAG_LENGTHS_BY_GENE, "Write fragment lengths by gene");
        configBuilder.addPath(BAM_FILE, true, "RNA BAM file location");
        configBuilder.addFlag(WRITE_EXON_DATA, "Exon region data");
        configBuilder.addFlag(WRITE_SPLICE_JUNC_DATA, "Write canonical splice junction counts");
        configBuilder.addFlag(WRITE_READ_DATA, "BAM read data");
        configBuilder.addFlag(WRITE_SPLICE_SITE_DATA, "Write support info for each splice site");
        configBuilder.addFlag(WRITE_TRANS_COMBO_DATA, "Write transcript group data for EM algo");
        configBuilder.addFlag(WRITE_FRAG_LENGTHS, "Write intronic fragment lengths to log");
        configBuilder.addFlag(WRITE_GC_DATA, "Write GC ratio counts from all genic reads");

        configBuilder.addPath(EXP_COUNTS_FILE, false, "File with generated expected expression rates per transcript");
        configBuilder.addPath(EXP_GC_RATIOS_FILE, false, "File with generated expected GC ratios per transcript");
        configBuilder.addPath(NEO_DIR_CFG, false, NEO_DIR_DESC);
        configBuilder.addPath(PANEL_TPM_NORM_FILE, false, "Panel TPM normalisation file");
        configBuilder.addInteger(READ_LENGTH, "Sample sequencing read length, if 0 then is inferred from reads", 0);
        configBuilder.addInteger(SINGLE_MAP_QUAL, "Map quality for reads mapped to a single location", DEFAULT_SINGLE_MAP_QUALITY);

        configBuilder.addConfigItem(ER_FRAGMENT_LENGTHS, false, ER_FRAGMENT_LENGTHS_DESC, DEFAULT_EXPECTED_RATE_LENGTHS);

        configBuilder.addFlag(RUN_VALIDATIONS, "Run auto-validations");
        configBuilder.addFlag(PERF_CHECKS, "Run performance logging routines");
        configBuilder.addPath(FILTER_READS_FILE, false, "Only process reads in this file");

        GeneRegionFilters.registerConfig(configBuilder);
        FusionConfig.registerConfig(configBuilder);
        addThreadOptions(configBuilder);
    }
}
