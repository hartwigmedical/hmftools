package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.linx.LinxConfig.formOutputPath;
import static com.hartwig.hmftools.svtools.common.ConfigUtils.DATA_OUTPUT_DIR;
import static com.hartwig.hmftools.svtools.common.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedRatesGenerator.FL_FREQUENCY;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedRatesGenerator.FL_LENGTH;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RnaExpConfig
{
    // config items
    public static final String GENE_TRANSCRIPTS_DIR = "gene_transcripts_dir";
    public static final String SAMPLE = "sample";

    private static final String GENE_ID_FILE = "gene_id_file";
    private static final String EXCLUDED_GENE_ID_FILE = "excluded_gene_id_file";
    private static final String CANONICAL_ONLY = "canonical_only";
    private static final String WRITE_EXON_DATA = "write_exon_data";
    private static final String WRITE_READ_DATA = "write_read_data";
    private static final String WRITE_TRANS_COMBO_DATA = "write_trans_combo_data";
    private static final String GENE_STATS_ONLY = "gene_stats_only";
    private static final String OUTPUT_ID = "output_id";

    public static final String REF_GENOME = "ref_genome";
    public static final String BAM_FILE = "bam_file";
    public static final String GC_BIAS_FILE = "gcbias_file";
    public static final String LONG_FRAGMENT_LIMIT = "long_frag_limit";
    public static final String KEEP_DUPLICATES = "keep_dups";
    public static final String MARK_DUPLICATES = "mark_dups";

    private static final String WRITE_FRAGMENT_LENGTHS = "write_frag_lengths";
    public static final String FRAG_LENGTH_MIN_COUNT = "frag_length_min_count";
    public static final String FRAG_LENGTHS_BY_GENE = "frag_length_by_gene";

    // expected expression config
    public static final String EXP_RATES_FILE = "exp_rates_file";
    public static final String APPLY_EXP_RATES = "apply_exp_rates";
    public static final String READ_LENGTH = "read_length";
    public static final String ER_FRAGMENT_LENGTHS = "exp_rate_frag_lengths";
    public static final String ER_CALC_FRAG_LENGTHS = "use_calc_frag_lengths";
    public static final String UNSPLICED_WEIGHT = "unspliced_weight";
    public static final String WRITE_EXPECTED_RATES = "write_exp_rates";

    public static final String SPECIFIC_TRANS_IDS = "specific_trans";
    public static final String SPECIFIC_CHR = "specific_chr";
    public static final String READ_COUNT_LIMIT = "read_count_limit";
    public static final String RUN_VALIDATIONS = "validate";
    public static final String THREADS = "threads";

    public final String SampleId;
    public final List<String> RestrictedGeneIds; // specific set of genes to process
    public final List<String> ExcludedGeneIds; // genes to ignore
    public final String OutputDir;
    public final String OutputIdentifier; // optionally include extra identifier in output files
    public final String GcBiasFile;
    public final boolean CanonicalTranscriptOnly;
    public final String BamFile;
    public final File RefGenomeFile;
    public IndexedFastaSequenceFile RefFastaSeqFile;
    public final int ReadCountLimit;
    public int MaxFragmentLength;
    public final boolean KeepDuplicates;
    public final boolean MarkDuplicates;

    public final boolean WriteExonData;
    public final boolean WriteReadData;
    public final boolean WriteFragmentLengths;
    public final boolean WriteTransComboData;

    public final String ExpRatesFile;
    public final boolean ApplyExpectedRates;
    public final boolean UseCalculatedFragmentLengths;
    public int ReadLength;
    public final List<int[]> ExpRateFragmentLengths;
    public final double UnsplicedWeight;
    public final boolean WriteExpectedRates;

    public final boolean GeneStatsOnly;
    public final int FragmentLengthMinCount;
    public final boolean FragmentLengthsByGene;

    public final List<String> SpecificTransIds;
    public final List<String> SpecificChromosomes;
    public final boolean RunValidations;
    public final int Threads;

    public static final int DEFAULT_MAX_READ_COUNT = 100000;
    public static final int DEFAULT_MAX_FRAGMENT_SIZE = 550;

    public static final int GENE_FRAGMENT_BUFFER = 1000; // width around a gene within which to search for reads

    public static final Logger RE_LOGGER = LogManager.getLogger(RnaExpConfig.class);

    public RnaExpConfig(final CommandLine cmd)
    {
        SampleId = cmd.getOptionValue(SAMPLE);

        RestrictedGeneIds = Lists.newArrayList();
        ExcludedGeneIds = Lists.newArrayList();

        if(cmd.hasOption(GENE_ID_FILE))
        {
            final String inputFile = cmd.getOptionValue(GENE_ID_FILE);
            loadGeneIdsFile(inputFile, RestrictedGeneIds);
            RE_LOGGER.info("file({}) load {} restricted genes", inputFile, RestrictedGeneIds.size());
        }

        if(cmd.hasOption(EXCLUDED_GENE_ID_FILE))
        {
            final String inputFile = cmd.getOptionValue(EXCLUDED_GENE_ID_FILE);
            loadGeneIdsFile(inputFile, ExcludedGeneIds);
            RE_LOGGER.info("file({}) load {} excluded genes", inputFile, ExcludedGeneIds.size());
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
                RE_LOGGER.debug("loading indexed fasta reference file");
                RefFastaSeqFile = new IndexedFastaSequenceFile(new File(refGenomeFilename));
            }
            catch (IOException e)
            {
                RE_LOGGER.error("Reference file loading failed: {}", e.toString());
            }
        }

        ReadCountLimit = Integer.parseInt(cmd.getOptionValue(READ_COUNT_LIMIT, "0"));
        MaxFragmentLength = Integer.parseInt(cmd.getOptionValue(LONG_FRAGMENT_LIMIT, String.valueOf(DEFAULT_MAX_FRAGMENT_SIZE)));
        KeepDuplicates = cmd.hasOption(KEEP_DUPLICATES);
        MarkDuplicates = cmd.hasOption(MARK_DUPLICATES);
        FragmentLengthMinCount = Integer.parseInt(cmd.getOptionValue(FRAG_LENGTH_MIN_COUNT, "0"));
        FragmentLengthsByGene = cmd.hasOption(FRAG_LENGTHS_BY_GENE);

        WriteExonData = cmd.hasOption(WRITE_EXON_DATA);
        WriteFragmentLengths = cmd.hasOption(WRITE_FRAGMENT_LENGTHS);
        WriteReadData = cmd.hasOption(WRITE_READ_DATA);
        WriteTransComboData = cmd.hasOption(WRITE_TRANS_COMBO_DATA);
        GeneStatsOnly = cmd.hasOption(GENE_STATS_ONLY);

        SpecificTransIds = cmd.hasOption(SPECIFIC_TRANS_IDS) ?
                Arrays.stream(cmd.getOptionValue(SPECIFIC_TRANS_IDS).split(";")).collect(Collectors.toList())
                : Lists.newArrayList();

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));
        RunValidations = cmd.hasOption(RUN_VALIDATIONS);
        SpecificChromosomes = cmd.hasOption(SPECIFIC_CHR) ? Arrays.stream(cmd.getOptionValue(SPECIFIC_CHR).split(";")).collect(Collectors.toList())
                : Lists.newArrayList();

        ApplyExpectedRates = cmd.hasOption(APPLY_EXP_RATES);
        ExpRatesFile = cmd.getOptionValue(EXP_RATES_FILE);

        WriteExpectedRates = cmd.hasOption(WRITE_EXPECTED_RATES);
        UseCalculatedFragmentLengths = cmd.hasOption(ER_CALC_FRAG_LENGTHS);
        ReadLength = Integer.parseInt(cmd.getOptionValue(READ_LENGTH, "0"));
        ExpRateFragmentLengths = Lists.newArrayList();
        UnsplicedWeight = 1; // Double.parseDouble(cmd.getOptionValue(UNSPLICED_WEIGHT, "1.0"));

        if(cmd.hasOption(ER_FRAGMENT_LENGTHS))
        {
            String[] fragLengths = cmd.getOptionValue(ER_FRAGMENT_LENGTHS).split(";");
            for(int i = 0; i < fragLengths.length; ++i)
            {
                String[] flItem = fragLengths[i].split("-");
                int fragLength = Integer.parseInt(flItem[FL_LENGTH]);
                int fragFrequency = Integer.parseInt(flItem[FL_FREQUENCY]);
                ExpRateFragmentLengths.add(new int[]{ fragLength, fragFrequency });
            }
        }
    }

    public boolean isValid()
    {
        if(OutputDir == null)
        {
            RE_LOGGER.error("not output directory specified");
            return false;
        }

        if(WriteExpectedRates)
        {
            if(ReadLength == 0 || ExpRateFragmentLengths.isEmpty())
            {
                RE_LOGGER.error("invalid read or fragment lengths for generating expected trans rates");
                return false;
            }

            return true;
        }

        if(BamFile == null || BamFile.isEmpty() || !Files.exists(Paths.get(BamFile)))
        {
            RE_LOGGER.error("BAM file({}) missing or not found", BamFile);
            return false;
        }

        if(SampleId == null || SampleId.isEmpty())
        {
            RE_LOGGER.error("sampleId missing");
            return false;
        }

        if(RefFastaSeqFile == null)
        {
            RE_LOGGER.error("ref genome missing");
            return false;
        }

        return true;
    }

    public boolean skipChromosome(final String chromosome)
    {
        return !SpecificChromosomes.isEmpty() && !SpecificChromosomes.contains(chromosome);
    }

    public String formOutputFile(final String fileId)
    {
        if(OutputIdentifier != null)
            return OutputDir + SampleId + "." + OutputIdentifier + "." + fileId;
        else
            return OutputDir + SampleId + "." + fileId;
    }

    public RnaExpConfig()
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
        ExpRateFragmentLengths = Lists.newArrayList();
        UnsplicedWeight = 1;
        ExpRatesFile = null;

        WriteExonData = false;
        WriteReadData = false;
        WriteFragmentLengths = false;
        WriteTransComboData = false;

        WriteExpectedRates = false;
        UseCalculatedFragmentLengths = false;
        OutputIdentifier = null;
        GeneStatsOnly = false;
        FragmentLengthsByGene = false;
        FragmentLengthMinCount = 0;
        SpecificTransIds = Lists.newArrayList();
        SpecificChromosomes = Lists.newArrayList();
        RunValidations = true;
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
        options.addOption(READ_COUNT_LIMIT, true, "Cap read-processing for genes with depth greater than this");
        options.addOption(GC_BIAS_FILE, true, "GC-bias file, generate if not found");
        options.addOption(REF_GENOME, true, "Ref genome file location");
        options.addOption(LONG_FRAGMENT_LIMIT, true, "Max RNA fragment size");
        options.addOption(KEEP_DUPLICATES, false, "Process duplicate reads (if marked as such eg by picard)");
        options.addOption(MARK_DUPLICATES, false, "Manually identify duplicate reads");
        options.addOption(FRAG_LENGTH_MIN_COUNT, true, "Fragment length measurement - min read fragments required");
        options.addOption(FRAG_LENGTHS_BY_GENE, false, "Write fragment lengths by gene");
        options.addOption(BAM_FILE, true, "RNA BAM file location");
        options.addOption(WRITE_EXON_DATA, false, "Exon region data");
        options.addOption(WRITE_READ_DATA, false, "BAM read data");
        options.addOption(WRITE_TRANS_COMBO_DATA, false, "Write transcript group data for EM algo");
        options.addOption(GENE_STATS_ONLY, false, "Skip all processing except gene summary data");
        options.addOption(WRITE_FRAGMENT_LENGTHS, false, "Write intronic fragment lengths to log");

        options.addOption(APPLY_EXP_RATES, false, "Generate expected expression rates for transcripts");
        options.addOption(EXP_RATES_FILE, true, "File with generated expected expression rates for transcripts");
        options.addOption(READ_LENGTH, true, "Sample sequencing read length (eg 76 or 151 bases");
        options.addOption(UNSPLICED_WEIGHT, true, "Weighting for unspliced expected fragments");
        options.addOption(ER_CALC_FRAG_LENGTHS, false, "Use sample fragment length distribution in expected rate calcs");

        options.addOption(ER_FRAGMENT_LENGTHS, true,
                "Fragment sizes and weights for expected transcript calcs (format: length1-freq1;length3-freq2 eg 100-10;150-20) in integer terms");

        options.addOption(WRITE_EXPECTED_RATES, false, "Write expected transcript rates to file");

        options.addOption(OUTPUT_ID, true, "Optionally add identifier to output files");
        options.addOption(SPECIFIC_TRANS_IDS, true, "List of transcripts separated by ';'");
        options.addOption(SPECIFIC_CHR, true, "Specify a single chromosome to analyse");
        options.addOption(THREADS, true, "Number of threads to use (default=0, single-threaded)");
        options.addOption(RUN_VALIDATIONS, false, "Run auto-validations");

        return options;
    }

    private static final int COL_GENE_ID = 0;

    private static void loadGeneIdsFile(final String filename, final List<String> geneIdList)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            RE_LOGGER.warn("invalid gene ID file({})", filename);
            return;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return;

            if (fileContents.get(COL_GENE_ID).contains("GeneId"))
            {
                // check for header row
                fileContents.remove(0);
            }

            geneIdList.addAll(fileContents.stream().map(x -> x.split(",")[COL_GENE_ID]).collect(Collectors.toList()));
        }
        catch (IOException e)
        {
            RE_LOGGER.warn("failed to load gene ID file({}): {}", filename, e.toString());
        }
    }

}
