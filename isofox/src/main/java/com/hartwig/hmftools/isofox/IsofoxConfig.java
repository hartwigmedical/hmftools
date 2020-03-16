package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.FL_FREQUENCY;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.FL_LENGTH;
import static com.hartwig.hmftools.linx.LinxConfig.formOutputPath;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator;
import com.hartwig.hmftools.isofox.gc.GcRatioCounts;

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

    private static final String GENE_ID_FILE = "gene_id_file";
    private static final String EXCLUDED_GENE_ID_FILE = "excluded_gene_id_file";
    private static final String CANONICAL_ONLY = "canonical_only";
    private static final String WRITE_TRANS_DATA = "write_trans_data";
    private static final String WRITE_EXON_DATA = "write_exon_data";
    private static final String WRITE_READ_DATA = "write_read_data";
    private static final String WRITE_TRANS_COMBO_DATA = "write_trans_combo_data";
    private static final String OUTPUT_ID = "output_id";

    public static final String REF_GENOME = "ref_genome";
    private static final String BAM_FILE = "bam_file";
    private static final String LONG_FRAGMENT_LIMIT = "long_frag_limit";
    private static final String KEEP_DUPLICATES = "keep_dups";
    private static final String MARK_DUPLICATES = "mark_dups";

    private static final String WRITE_FRAG_LENGTHS = "write_frag_lengths";
    private static final String WRITE_FRAG_LENGTHS_ONLY = "write_frag_lengths_only";
    private static final String FRAG_LENGTH_MIN_COUNT = "frag_length_min_count";
    private static final String FRAG_LENGTHS_BY_GENE = "frag_length_by_gene";
    private static final String WRITE_FRAG_READS = "write_frag_length_reads";

    private static final String GC_BIAS_FILE = "gcbias_file";
    private static final String WRITE_READ_GC_RATIOS = "write_read_gc_ratios";
    private static final String GC_RATIO_BUCKET_SIZE = "gc_ratio_bucket";

    // expected expression config
    private static final String EXP_COUNTS_FILE = "exp_counts_file";
    private static final String APPLY_EXP_RATES = "apply_exp_rates";
    private static final String READ_LENGTH = "read_length";
    private static final String ER_FRAGMENT_LENGTHS = "exp_rate_frag_lengths";
    private static final String ER_CALC_FRAG_LENGTHS = "use_calc_frag_lengths";
    private static final String UNSPLICED_WEIGHT = "unspliced_weight";
    private static final String WRITE_EXPECTED_RATES = "write_exp_rates";
    private static final String WRITE_EXPECTED_COUNTS = "write_exp_counts";

    private static final String SPECIFIC_TRANS_IDS = "specific_trans";
    private static final String SPECIFIC_CHR = "specific_chr";
    private static final String READ_COUNT_LIMIT = "read_count_limit";
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
    public final boolean CanonicalTranscriptOnly;
    public final String BamFile;
    public final File RefGenomeFile;
    public IndexedFastaSequenceFile RefFastaSeqFile;
    public final int ReadCountLimit;
    public int MaxFragmentLength;
    public final boolean KeepDuplicates;
    public final boolean MarkDuplicates;

    public final boolean WriteTransData;
    public final boolean WriteExonData;
    public final boolean WriteReadData;
    public final boolean WriteTransComboData;

    public final String ExpCountsFile;
    public final boolean ApplyExpectedRates;
    public final boolean UseCalculatedFragmentLengths;
    public int ReadLength;
    public final List<int[]> FragmentLengthData;
    public final double UnsplicedWeight;
    public final boolean WriteExpectedRates;
    public final boolean WriteExpectedCounts;

    public final boolean WriteFragmentLengths;
    public final int FragmentLengthMinCount;
    public final boolean FragmentLengthsByGene;
    public final boolean WriteFragmentLengthsOnly;
    public final boolean WriteFragmentReads;

    public final boolean WriteReadGcRatios;
    public final String GcBiasFile;
    public static double GC_RATIO_BUCKET = GcRatioCounts.DEFAULT_GC_RATIO_BUCKET;

    public final List<String> SpecificTransIds;
    public final List<String> SpecificChromosomes;
    public final boolean RunValidations;
    public final boolean RunPerfChecks;
    public final int Threads;

    public static final int DEFAULT_MAX_READ_COUNT = 100000;
    public static final int DEFAULT_MAX_FRAGMENT_SIZE = 550;

    public static final int GENE_FRAGMENT_BUFFER = 1000; // width around a gene within which to search for reads

    public static final Logger ISF_LOGGER = LogManager.getLogger(IsofoxConfig.class);

    public IsofoxConfig(final CommandLine cmd)
    {
        SampleId = cmd.getOptionValue(SAMPLE);

        RestrictedGeneIds = Lists.newArrayList();
        ExcludedGeneIds = Lists.newArrayList();

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

        CanonicalTranscriptOnly = cmd.hasOption(CANONICAL_ONLY);

        OutputDir = formOutputPath(cmd.getOptionValue(DATA_OUTPUT_DIR));
        OutputIdentifier = cmd.getOptionValue(OUTPUT_ID);

        BamFile = cmd.getOptionValue(BAM_FILE);
        GcBiasFile = cmd.getOptionValue(GC_BIAS_FILE, "");

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

        ReadCountLimit = Integer.parseInt(cmd.getOptionValue(READ_COUNT_LIMIT, "0"));
        MaxFragmentLength = Integer.parseInt(cmd.getOptionValue(LONG_FRAGMENT_LIMIT, String.valueOf(DEFAULT_MAX_FRAGMENT_SIZE)));
        KeepDuplicates = cmd.hasOption(KEEP_DUPLICATES);
        MarkDuplicates = cmd.hasOption(MARK_DUPLICATES);
        FragmentLengthMinCount = Integer.parseInt(cmd.getOptionValue(FRAG_LENGTH_MIN_COUNT, "0"));
        FragmentLengthsByGene = cmd.hasOption(FRAG_LENGTHS_BY_GENE);

        WriteExonData = cmd.hasOption(WRITE_EXON_DATA);
        WriteFragmentLengths = cmd.hasOption(WRITE_FRAG_LENGTHS);
        WriteFragmentReads = cmd.hasOption(WRITE_FRAG_READS);
        WriteFragmentLengthsOnly = cmd.hasOption(WRITE_FRAG_LENGTHS_ONLY);
        WriteReadData = cmd.hasOption(WRITE_READ_DATA);
        WriteTransComboData = cmd.hasOption(WRITE_TRANS_COMBO_DATA);
        WriteReadGcRatios = cmd.hasOption(WRITE_READ_GC_RATIOS);
        WriteTransData = Boolean.parseBoolean(cmd.getOptionValue(WRITE_TRANS_DATA, "true"));

        GC_RATIO_BUCKET = cmd.hasOption(GC_RATIO_BUCKET_SIZE) ?
                Double.parseDouble(cmd.getOptionValue(GC_RATIO_BUCKET_SIZE)) : GcRatioCounts.DEFAULT_GC_RATIO_BUCKET;

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

        WriteExpectedCounts = cmd.hasOption(WRITE_EXPECTED_COUNTS);
        WriteExpectedRates = !WriteExpectedCounts && cmd.hasOption(WRITE_EXPECTED_RATES);
        UseCalculatedFragmentLengths = cmd.hasOption(ER_CALC_FRAG_LENGTHS);
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

        if(WriteExpectedCounts)
        {
            if(ReadLength == 0 || FragmentLengthData.isEmpty())
            {
                ISF_LOGGER.error("invalid read or fragment lengths for generating expected trans rates");
                return false;
            }

            return true;
        }

        if(ApplyExpectedRates && ExpCountsFile == null)
        {
            if(!UseCalculatedFragmentLengths && (ReadLength == 0 || FragmentLengthData.isEmpty()))
            {
                ISF_LOGGER.error("invalid read or fragment lengths for generating expected trans rates");
                return false;
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

        if(WriteFragmentLengthsOnly && FragmentLengthMinCount == 0)
        {
            ISF_LOGGER.error("min frag count missing for frag length distribution logging");
            return false;
        }

        if(RestrictedGeneIds.isEmpty() && (WriteExonData || WriteReadData))
        {
            ISF_LOGGER.warn("writing exon and/or read data for all transcripts may be slow and generate large output files");
        }

        return true;
    }

    public boolean skipChromosome(final String chromosome)
    {
        return !SpecificChromosomes.isEmpty() && !SpecificChromosomes.contains(chromosome);
    }

    public boolean generateExpRatesOnly()
    {
        return WriteExpectedCounts && !UseCalculatedFragmentLengths && !ApplyExpectedRates;
    }

    public boolean writeExpectedRateData()
    {
        return WriteExpectedCounts || WriteExpectedRates;
    }

    public boolean requireFragmentLengthCalcs()
    {
        return WriteFragmentLengths || UseCalculatedFragmentLengths;
    }

    public String formOutputFile(final String fileId)
    {
        if(OutputIdentifier != null)
            return OutputDir + SampleId + "." + OutputIdentifier + "." + fileId;
        else
            return OutputDir + SampleId + "." + fileId;
    }

    public IsofoxConfig()
    {
        SampleId = "TEST";

        RestrictedGeneIds = Lists.newArrayList();
        ExcludedGeneIds = Lists.newArrayList();
        OutputDir = "";
        BamFile = null;
        RefGenomeFile = null;
        RefFastaSeqFile = null;
        CanonicalTranscriptOnly = false;
        ReadCountLimit = DEFAULT_MAX_READ_COUNT;
        GcBiasFile = "";
        MaxFragmentLength = DEFAULT_MAX_FRAGMENT_SIZE;
        KeepDuplicates = false;
        MarkDuplicates = false;

        ApplyExpectedRates = false;
        ReadLength = 0;
        FragmentLengthData = Lists.newArrayList();
        UnsplicedWeight = 1;
        ExpCountsFile = null;

        WriteTransData = true;
        WriteExonData = false;
        WriteReadData = false;
        WriteFragmentLengths = false;
        WriteFragmentLengthsOnly = false;
        WriteTransComboData = false;
        WriteFragmentReads = false;
        WriteReadGcRatios = false;

        WriteExpectedRates = false;
        WriteExpectedCounts = false;
        UseCalculatedFragmentLengths = false;
        OutputIdentifier = null;
        FragmentLengthsByGene = false;
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

        options.addOption(CANONICAL_ONLY, false, "Check all transcripts, not just canonical");
        options.addOption(GENE_ID_FILE, true, "Optional CSV file of genes to analyse");
        options.addOption(EXCLUDED_GENE_ID_FILE, true, "Optional CSV file of genes to ignore");
        options.addOption(DATA_OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        options.addOption(LOG_LEVEL, true, "Logging: INFO(default), DEBUG or TRACE (verbose)");
        options.addOption(READ_COUNT_LIMIT, true, "Cap read-processing for genes with depth greater than this");
        options.addOption(REF_GENOME, true, "Ref genome file location");
        options.addOption(LONG_FRAGMENT_LIMIT, true, "Max RNA fragment size");
        options.addOption(KEEP_DUPLICATES, false, "Process duplicate reads (if marked as such eg by picard)");
        options.addOption(MARK_DUPLICATES, false, "Manually identify duplicate reads");
        options.addOption(FRAG_LENGTH_MIN_COUNT, true, "Fragment length measurement - min read fragments required");
        options.addOption(FRAG_LENGTHS_BY_GENE, false, "Write fragment lengths by gene");
        options.addOption(WRITE_FRAG_READS, false, "Write fragment read data from length determination");
        options.addOption(BAM_FILE, true, "RNA BAM file location");
        options.addOption(WRITE_TRANS_DATA, true, "Produce transcript-level counts data (default=true)");
        options.addOption(WRITE_EXON_DATA, false, "Exon region data");
        options.addOption(WRITE_READ_DATA, false, "BAM read data");
        options.addOption(WRITE_TRANS_COMBO_DATA, false, "Write transcript group data for EM algo");
        options.addOption(WRITE_FRAG_LENGTHS, false, "Write intronic fragment lengths to log");
        options.addOption(WRITE_FRAG_LENGTHS_ONLY, false, "Only write intronic fragment lengths then exit");

        options.addOption(WRITE_READ_GC_RATIOS, false, "Write GC Ratio counts from all genic reads");
        options.addOption(GC_BIAS_FILE, true, "GC-bias file, generate if not found");
        options.addOption(GC_RATIO_BUCKET_SIZE, true, "Rounding size for GC-calcs (default=0.01");

        options.addOption(APPLY_EXP_RATES, false, "Generate expected expression rates for transcripts");
        options.addOption(EXP_COUNTS_FILE, true, "File with generated expected expression rates for transcripts");
        options.addOption(READ_LENGTH, true, "Sample sequencing read length (eg 76 or 151 bases");
        options.addOption(UNSPLICED_WEIGHT, true, "Weighting for unspliced expected fragments");
        options.addOption(ER_CALC_FRAG_LENGTHS, false, "Use sample fragment length distribution in expected rate calcs");

        options.addOption(ER_FRAGMENT_LENGTHS, true,
                "Fragment sizes and weights for expected transcript calcs (format: length1-freq1;length3-freq2 eg 100-10;150-20) in integer terms");

        options.addOption(WRITE_EXPECTED_RATES, false, "Write sample expected expression rates to file");
        options.addOption(WRITE_EXPECTED_COUNTS, false, "Write expected expression counts from common frag lengths to file");

        options.addOption(OUTPUT_ID, true, "Optionally add identifier to output files");
        options.addOption(SPECIFIC_TRANS_IDS, true, "List of transcripts separated by ';'");
        options.addOption(SPECIFIC_CHR, true, "Specify a single chromosome to analyse");
        options.addOption(THREADS, true, "Number of threads to use (default=0, single-threaded)");
        options.addOption(RUN_VALIDATIONS, false, "Run auto-validations");
        options.addOption(PERF_CHECKS, false, "Run performance logging routines");

        return options;
    }

    private static final int COL_GENE_ID = 0;

    private static void loadGeneIdsFile(final String filename, final List<String> geneIdList)
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

            geneIdList.addAll(fileContents.stream().map(x -> x.split(",")[COL_GENE_ID]).collect(Collectors.toList()));
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load gene ID file({}): {}", filename, e.toString());
        }
    }

}
