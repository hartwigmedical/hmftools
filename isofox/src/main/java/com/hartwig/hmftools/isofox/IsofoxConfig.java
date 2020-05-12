package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_FRAG_LENGTH_MIN_COUNT;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_GC_RATIO_BUCKET;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_MAX_FRAGMENT_SIZE;
import static com.hartwig.hmftools.isofox.IsofoxFunction.EXPECTED_GC_COUNTS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.EXPECTED_TRANS_COUNTS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.NOVEL_LOCATIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.FL_FREQUENCY;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.FL_LENGTH;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator;
import com.hartwig.hmftools.isofox.fusion.FusionConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class IsofoxConfig
{
    // config items
    public static final String GENE_TRANSCRIPTS_DIR = "gene_transcripts_dir";
    public static final String SAMPLE = "sample";
    public static final String DATA_OUTPUT_DIR = "output_dir";
    public static final String OUTPUT_ID = "output_id";
    public static final String FUNCTIONS = "functions";

    public static final String GENE_ID_FILE = "gene_id_file";
    public static final String EXCLUDED_GENE_ID_FILE = "excluded_gene_id_file";
    private static final String ENRICHED_GENE_IDS = "enriched_gene_ids";
    private static final String CANONICAL_ONLY = "canonical_only";
    private static final String WRITE_EXON_DATA = "write_exon_data";
    private static final String WRITE_READ_DATA = "write_read_data";
    private static final String WRITE_TRANS_COMBO_DATA = "write_trans_combo_data";

    public static final String REF_GENOME = "ref_genome";
    private static final String BAM_FILE = "bam_file";
    private static final String LONG_FRAGMENT_LIMIT = "long_frag_limit";
    private static final String DROP_DUPLICATES = "drop_dups";
    private static final String MARK_DUPLICATES = "mark_dups";

    private static final String WRITE_FRAG_LENGTHS = "write_frag_lengths";
    private static final String FRAG_LENGTH_MIN_COUNT = "frag_length_min_count";
    private static final String FRAG_LENGTHS_BY_GENE = "frag_length_by_gene";

    private static final String WRITE_GC_DATA = "write_gc_data";
    private static final String GC_RATIO_BUCKET_SIZE = "gc_ratio_bucket";
    private static final String GC_ADJUSTS_FILE = "gc_adjusts_file";

    // expected expression config
    private static final String EXP_COUNTS_FILE = "exp_counts_file";
    private static final String EXP_GC_RATIOS_FILE = "exp_gc_ratios_file";
    private static final String APPLY_EXP_RATES = "apply_exp_rates";
    private static final String READ_LENGTH = "read_length";
    private static final String ER_FRAGMENT_LENGTHS = "exp_rate_frag_lengths";
    private static final String APPLY_FRAG_LENGTH_ADJUSTMENT = "apply_calc_frag_lengths";
    private static final String APPLY_GC_BIAS_ADJUSTMENT = "apply_gc_bias_adjust";
    private static final String UNSPLICED_WEIGHT = "unspliced_weight";
    private static final String WRITE_EXPECTED_RATES = "write_exp_rates";

    private static final String SPECIFIC_TRANS_IDS = "specific_trans";
    private static final String SPECIFIC_CHR = "specific_chr";
    private static final String GENE_READ_LIMIT = "gene_read_limit";
    private static final String RUN_VALIDATIONS = "validate";
    private static final String PERF_CHECKS = "perf_checks";
    private static final String THREADS = "threads";
    public static final String LOG_DEBUG = "log_debug";
    public static final String LOG_LEVEL = "log_level";

    public final String SampleId;
    public final List<String> RestrictedGeneIds; // specific set of genes to process
    public final List<String> ExcludedGeneIds; // genes to ignore
    public final String OutputDir;
    public final String OutputIdentifier; // optionally include extra identifier in output files

    public final List<IsofoxFunction> Functions;

    public final boolean CanonicalTranscriptOnly;
    public final String BamFile;
    public final File RefGenomeFile;
    public IndexedFastaSequenceFile RefFastaSeqFile;
    public final int GeneReadLimit;
    public int MaxFragmentLength;
    public final boolean DropDuplicates;
    public final boolean MarkDuplicates;

    public final boolean WriteExonData;
    public final boolean WriteReadData;
    public final boolean WriteTransComboData;

    public final String ExpCountsFile;
    public final String ExpGcRatiosFile;
    public final boolean ApplyExpectedRates;
    public final boolean ApplyFragmentLengthAdjust;
    public final boolean ApplyGcBiasAdjust;
    public int ReadLength;
    public final List<int[]> FragmentLengthData;
    public final double UnsplicedWeight;
    public final boolean WriteExpectedRates;

    public final boolean WriteFragmentLengths;
    public final int FragmentLengthMinCount;
    public final boolean WriteFragmentLengthsByGene;

    public final boolean WriteGcData;
    public final String GcAdjustmentsFile;
    public static double GC_RATIO_BUCKET = DEFAULT_GC_RATIO_BUCKET;

    public final List<String> EnrichedGeneIds;

    public final FusionConfig Fusions;

    public final List<String> SpecificTransIds;
    public final List<String> SpecificChromosomes;
    public final boolean RunValidations;
    public final boolean RunPerfChecks;
    public final int Threads;

    public static final Logger ISF_LOGGER = LogManager.getLogger(IsofoxConfig.class);

    public IsofoxConfig(final CommandLine cmd)
    {
        SampleId = cmd.getOptionValue(SAMPLE);

        RestrictedGeneIds = Lists.newArrayList();
        ExcludedGeneIds = Lists.newArrayList();

        EnrichedGeneIds = cmd.hasOption(ENRICHED_GENE_IDS) ?
                Arrays.stream(cmd.getOptionValue(ENRICHED_GENE_IDS).split(";")).collect(Collectors.toList()) : Lists.newArrayList();

        if(cmd.hasOption(GENE_ID_FILE))
        {
            final String inputFile = cmd.getOptionValue(GENE_ID_FILE);
            loadGeneIdsFile(inputFile, RestrictedGeneIds);
            ISF_LOGGER.info("file({}) loaded {} restricted genes", inputFile, RestrictedGeneIds.size());
        }

        if(cmd.hasOption(EXCLUDED_GENE_ID_FILE))
        {
            final String inputFile = cmd.getOptionValue(EXCLUDED_GENE_ID_FILE);
            loadGeneIdsFile(inputFile, ExcludedGeneIds);
            ISF_LOGGER.info("file({}) loaded {} excluded genes", inputFile, ExcludedGeneIds.size());
        }

        Functions = Lists.newArrayList();

        if(cmd.hasOption(FUNCTIONS))
        {
            final String[] functionsStr = cmd.getOptionValue(FUNCTIONS).split(";");

            for(final String functionStr : functionsStr)
            {
                Functions.add(IsofoxFunction.valueOf(functionStr));
                ISF_LOGGER.info("running function(s): {}", functionStr);
            }
        }
        else
        {
            Functions.add(TRANSCRIPT_COUNTS);
            Functions.add(NOVEL_LOCATIONS);
        }

        CanonicalTranscriptOnly = cmd.hasOption(CANONICAL_ONLY);

        String outputdir = cmd.getOptionValue(DATA_OUTPUT_DIR);
        if(!outputdir.endsWith(File.separator))
            outputdir += File.separator;
        OutputDir = outputdir;
        OutputIdentifier = cmd.getOptionValue(OUTPUT_ID);

        BamFile = cmd.getOptionValue(BAM_FILE);
        GcAdjustmentsFile = cmd.getOptionValue(GC_ADJUSTS_FILE);

        final String refGenomeFilename = cmd.getOptionValue(REF_GENOME);
        RefGenomeFile = refGenomeFilename != null ? new File(refGenomeFilename) : null;

        RefFastaSeqFile = null;

        if(RefGenomeFile != null)
        {
            try
            {
                ISF_LOGGER.debug("loading indexed fasta reference file");
                RefFastaSeqFile = new IndexedFastaSequenceFile(new File(refGenomeFilename));
            }
            catch (IOException e)
            {
                ISF_LOGGER.error("Reference file loading failed: {}", e.toString());
            }
        }

        GeneReadLimit = Integer.parseInt(cmd.getOptionValue(GENE_READ_LIMIT, "0"));
        MaxFragmentLength = Integer.parseInt(cmd.getOptionValue(LONG_FRAGMENT_LIMIT, String.valueOf(DEFAULT_MAX_FRAGMENT_SIZE)));
        DropDuplicates = cmd.hasOption(DROP_DUPLICATES);
        MarkDuplicates = cmd.hasOption(MARK_DUPLICATES);
        FragmentLengthMinCount = Integer.parseInt(cmd.getOptionValue(FRAG_LENGTH_MIN_COUNT, String.valueOf(DEFAULT_FRAG_LENGTH_MIN_COUNT)));

        WriteExonData = cmd.hasOption(WRITE_EXON_DATA);
        WriteFragmentLengths = cmd.hasOption(WRITE_FRAG_LENGTHS);
        WriteFragmentLengthsByGene = cmd.hasOption(FRAG_LENGTHS_BY_GENE);
        WriteReadData = cmd.hasOption(WRITE_READ_DATA);
        WriteTransComboData = cmd.hasOption(WRITE_TRANS_COMBO_DATA);
        WriteGcData = cmd.hasOption(WRITE_GC_DATA);

        GC_RATIO_BUCKET = cmd.hasOption(GC_RATIO_BUCKET_SIZE) ?
                Double.parseDouble(cmd.getOptionValue(GC_RATIO_BUCKET_SIZE)) : DEFAULT_GC_RATIO_BUCKET;

        SpecificTransIds = cmd.hasOption(SPECIFIC_TRANS_IDS) ?
                Arrays.stream(cmd.getOptionValue(SPECIFIC_TRANS_IDS).split(";")).collect(Collectors.toList())
                : Lists.newArrayList();

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));
        RunValidations = cmd.hasOption(RUN_VALIDATIONS);
        RunPerfChecks = cmd.hasOption(PERF_CHECKS);
        SpecificChromosomes = cmd.hasOption(SPECIFIC_CHR) ? Arrays.stream(cmd.getOptionValue(SPECIFIC_CHR).split(";")).collect(Collectors.toList())
                : Lists.newArrayList();

        ApplyExpectedRates = cmd.hasOption(APPLY_EXP_RATES);
        ExpCountsFile = cmd.getOptionValue(EXP_COUNTS_FILE);
        ExpGcRatiosFile = cmd.getOptionValue(EXP_GC_RATIOS_FILE);

        WriteExpectedRates = cmd.hasOption(WRITE_EXPECTED_RATES);
        ApplyFragmentLengthAdjust = cmd.hasOption(APPLY_FRAG_LENGTH_ADJUSTMENT);
        ApplyGcBiasAdjust = cmd.hasOption(APPLY_GC_BIAS_ADJUSTMENT);
        ReadLength = Integer.parseInt(cmd.getOptionValue(READ_LENGTH, "0"));
        FragmentLengthData = Lists.newArrayList();
        UnsplicedWeight = 1; // Double.parseDouble(cmd.getOptionValue(UNSPLICED_WEIGHT, "1.0"));

        if(cmd.hasOption(ER_FRAGMENT_LENGTHS))
        {
            String[] fragLengths = cmd.getOptionValue(ER_FRAGMENT_LENGTHS).split(";");
            for(int i = 0; i < fragLengths.length; ++i)
            {
                String[] flItem = fragLengths[i].split("-");
                int fragLength = Integer.parseInt(flItem[FL_LENGTH]);
                int fragFrequency = Integer.parseInt(flItem[ExpectedRatesGenerator.FL_FREQUENCY]);
                FragmentLengthData.add(new int[]{ fragLength, fragFrequency });
            }
        }

        Fusions = new FusionConfig(cmd);
    }

    public boolean isValid()
    {
        if(OutputDir == null)
        {
            ISF_LOGGER.error("not output directory specified");
            return false;
        }

        if(!FragmentLengthData.isEmpty())
        {
            if(FragmentLengthData.stream().anyMatch(x -> x[FL_LENGTH] <= 0 || x[FL_FREQUENCY] < 0))
            {
                ISF_LOGGER.error("invalid fragment lengths");
                return false;
            }
        }

        if(runFunction(EXPECTED_TRANS_COUNTS))
        {
            if(ReadLength == 0 || FragmentLengthData.isEmpty())
            {
                ISF_LOGGER.error("invalid read or fragment lengths for generating expected trans counts");
                return false;
            }

            return true;
        }

        if(runFunction(EXPECTED_GC_COUNTS))
        {
            if(ReadLength == 0 || RefFastaSeqFile == null)
            {
                ISF_LOGGER.error("invalid read length or ref genome for generating expected GC ratio counts");
                return false;
            }

            return true;
        }

        if(ApplyExpectedRates && ExpCountsFile == null)
        {
            if(ApplyFragmentLengthAdjust)
            {
                if(FragmentLengthMinCount == 0)
                {
                    ISF_LOGGER.error("invalid min fragment length count for fragment length sampling");
                    return false;
                }
            }
            else
            {
                if(ReadLength == 0 || FragmentLengthData.isEmpty())
                {
                    ISF_LOGGER.error("invalid read or fragment lengths for generating expected trans rates");
                    return false;
                }
            }

            return true;
        }

        if(BamFile == null || BamFile.isEmpty() || !Files.exists(Paths.get(BamFile)))
        {
            ISF_LOGGER.error("BAM file({}) missing or not found", BamFile);
            return false;
        }

        if(SampleId == null || SampleId.isEmpty())
        {
            ISF_LOGGER.error("sampleId missing");
            return false;
        }

        if(RefFastaSeqFile == null)
        {
            ISF_LOGGER.error("ref genome missing");
            return false;
        }

        if(RestrictedGeneIds.isEmpty() && (WriteExonData || WriteReadData))
        {
            ISF_LOGGER.warn("writing exon and/or read data for all transcripts may be slow and generate large output files");
        }

        return true;
    }

    public static boolean validConfigPaths(final CommandLine cmd)
    {
        return configPathValid(cmd, DATA_OUTPUT_DIR) && configPathValid(cmd, REF_GENOME)  && configPathValid(cmd, GENE_TRANSCRIPTS_DIR)
                && configPathValid(cmd, GENE_ID_FILE) && configPathValid(cmd, EXCLUDED_GENE_ID_FILE) && configPathValid(cmd, GC_ADJUSTS_FILE)
                && configPathValid(cmd, BAM_FILE) && configPathValid(cmd, EXP_COUNTS_FILE) && configPathValid(cmd, EXP_GC_RATIOS_FILE);
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

    public boolean skipChromosome(final String chromosome)
    {
        return !SpecificChromosomes.isEmpty() && !SpecificChromosomes.contains(chromosome);
    }

    public boolean generateExpectedDataOnly()
    {
        return runFunction(EXPECTED_GC_COUNTS) || runFunction(EXPECTED_TRANS_COUNTS);
    }

    public boolean requireFragmentLengthCalcs()
    {
        return WriteFragmentLengths || ApplyFragmentLengthAdjust;
    }
    public boolean requireGcRatioCalcs() { return WriteGcData || ApplyGcBiasAdjust; }

    public String formOutputFile(final String fileId)
    {
        if(OutputIdentifier != null)
            return OutputDir + SampleId + ".isf." + OutputIdentifier + "." + fileId;
        else
            return OutputDir + SampleId + ".isf." + fileId;
    }

    public IsofoxConfig()
    {
        SampleId = "TEST";

        Functions = Lists.newArrayList();
        Functions.add(TRANSCRIPT_COUNTS);
        Functions.add(NOVEL_LOCATIONS);

        RestrictedGeneIds = Lists.newArrayList();
        ExcludedGeneIds = Lists.newArrayList();
        EnrichedGeneIds = Lists.newArrayList();
        OutputDir = null;
        BamFile = null;
        RefGenomeFile = null;
        RefFastaSeqFile = null;
        CanonicalTranscriptOnly = false;
        GeneReadLimit = 0;
        GcAdjustmentsFile = "";
        MaxFragmentLength = DEFAULT_MAX_FRAGMENT_SIZE;
        DropDuplicates = false;
        MarkDuplicates = false;

        ApplyExpectedRates = false;
        ReadLength = 0;
        FragmentLengthData = Lists.newArrayList();
        UnsplicedWeight = 1;
        ExpCountsFile = null;
        ExpGcRatiosFile = null;

        WriteExonData = false;
        WriteReadData = false;
        WriteFragmentLengths = false;
        WriteTransComboData = false;
        WriteGcData = false;
        Fusions = new FusionConfig();

        WriteExpectedRates = false;
        ApplyFragmentLengthAdjust = false;
        ApplyGcBiasAdjust = false;
        OutputIdentifier = null;
        WriteFragmentLengthsByGene = false;
        FragmentLengthMinCount = 0;

        SpecificTransIds = Lists.newArrayList();
        SpecificChromosomes = Lists.newArrayList();
        RunValidations = true;
        RunPerfChecks = false;
        Threads = 0;
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample ID");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Path to Ensembl data cache");

        options.addOption(FUNCTIONS, true, "Optional: list of functional routines to run (see documentation)");
        options.addOption(CANONICAL_ONLY, false, "Check all transcripts, not just canonical");
        options.addOption(GENE_ID_FILE, true, "Optional CSV file of genes to analyse");
        options.addOption(EXCLUDED_GENE_ID_FILE, true, "Optional CSV file of genes to ignore");
        options.addOption(ENRICHED_GENE_IDS, true, "Optional list of geneIds to treat as enriched");
        options.addOption(DATA_OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        options.addOption(LOG_LEVEL, true, "Logging: INFO(default), DEBUG or TRACE (verbose)");
        options.addOption(GENE_READ_LIMIT, true, "Per-gene limit on max reads processed (default=0, not applied)");
        options.addOption(REF_GENOME, true, "Ref genome file location");
        options.addOption(LONG_FRAGMENT_LIMIT, true, "Max RNA fragment size");
        options.addOption(DROP_DUPLICATES, false, "Include duplicate fragments in expression calculations");
        options.addOption(MARK_DUPLICATES, false, "Manually identify duplicate fragments");
        options.addOption(FRAG_LENGTH_MIN_COUNT, true, "Fragment length measurement - min read fragments required");
        options.addOption(FRAG_LENGTHS_BY_GENE, false, "Write fragment lengths by gene");
        options.addOption(BAM_FILE, true, "RNA BAM file location");
        options.addOption(WRITE_EXON_DATA, false, "Exon region data");
        options.addOption(WRITE_READ_DATA, false, "BAM read data");
        options.addOption(WRITE_TRANS_COMBO_DATA, false, "Write transcript group data for EM algo");
        options.addOption(WRITE_FRAG_LENGTHS, false, "Write intronic fragment lengths to log");

        options.addOption(WRITE_GC_DATA, false, "Write GC ratio counts from all genic reads");
        options.addOption(GC_ADJUSTS_FILE, true, "GC-bias file, generate if not found");
        options.addOption(GC_RATIO_BUCKET_SIZE, true, "Rounding size for GC-calcs (default=0.01");

        options.addOption(APPLY_EXP_RATES, false, "Generate expected expression rates for transcripts");
        options.addOption(EXP_COUNTS_FILE, true, "File with generated expected expression rates per transcript");
        options.addOption(EXP_GC_RATIOS_FILE, true, "File with generated expected GC ratios per transcript");
        options.addOption(READ_LENGTH, true, "Sample sequencing read length (eg 76 or 151 bases");
        options.addOption(UNSPLICED_WEIGHT, true, "Weighting for unspliced expected fragments");
        options.addOption(APPLY_FRAG_LENGTH_ADJUSTMENT, false, "Use sample fragment length distribution in expected rate calcs");
        options.addOption(APPLY_GC_BIAS_ADJUSTMENT, false, "Use GC Bias adjustments in expected rate calcs");

        options.addOption(ER_FRAGMENT_LENGTHS, true,
                "Fragment sizes and weights for expected transcript calcs (format: length1-freq1;length3-freq2 eg 100-10;150-20) in integer terms");

        options.addOption(WRITE_EXPECTED_RATES, false, "Write sample expected expression rates to file");

        options.addOption(OUTPUT_ID, true, "Optionally add identifier to output files");
        options.addOption(SPECIFIC_TRANS_IDS, true, "List of transcripts separated by ';'");
        options.addOption(SPECIFIC_CHR, true, "Specify a single chromosome to analyse");
        options.addOption(THREADS, true, "Number of threads to use (default=0, single-threaded)");
        options.addOption(RUN_VALIDATIONS, false, "Run auto-validations");
        options.addOption(PERF_CHECKS, false, "Run performance logging routines");

        FusionConfig.addCommandLineOptions(options);

        return options;
    }

    private static final int COL_GENE_ID = 0;

    public static void loadGeneIdsFile(final String filename, final List<String> geneIdList)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            ISF_LOGGER.warn("invalid gene ID file({})", filename);
            return;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return;

            if (fileContents.get(0).contains("GeneId"))
            {
                // check for header row
                fileContents.remove(0);
            }

            geneIdList.addAll(fileContents.stream()
                    .filter(x -> !x.contains("GeneId"))
                    .filter(x -> !x.startsWith("#"))
                    .map(x -> x.split(",")[COL_GENE_ID]).collect(Collectors.toList()));
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load gene ID file({}): {}", filename, e.toString());
        }
    }

}
