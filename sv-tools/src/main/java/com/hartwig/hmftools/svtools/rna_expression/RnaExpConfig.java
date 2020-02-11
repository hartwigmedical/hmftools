package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.svtools.common.ConfigUtils.DATA_OUTPUT_DIR;
import static com.hartwig.hmftools.svtools.common.ConfigUtils.LOG_DEBUG;

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
    private static final String ALL_TRANSCRIPTS = "all_transcripts";
    private static final String WRITE_EXON_DATA = "write_exon_data";
    private static final String WRITE_FRAGMENT_LENGTHS = "write_frag_lengths";
    private static final String WRITE_READ_DATA = "write_read_data";
    private static final String GENE_STATS_ONLY = "gene_stats_only";

    public static final String REF_GENOME = "ref_genome";
    public static final String BAM_FILE = "bam_file";
    public static final String GC_BIAS_FILE = "gcbias_file";
    public static final String READ_COUNT_LIMIT = "read_count_limit";
    public static final String LONG_FRAGMENT_LIMIT = "long_frag_limit";
    public static final String KEEP_DUPLICATES = "keep_dups";
    public static final String MARK_DUPLICATES = "mark_dups";
    public static final String FRAG_LENGTH_MIN_COUNT = "frag_length_min_count";
    public static final String FRAG_LENGTHS_BY_GENE = "frag_length_by_gene";

    public static final String SPECIFIC_TRANS_IDS = "specific_trans";

    public final List<String> RestrictedGeneIds;
    public final String OutputDir;
    public final String GcBiasFile;
    public final boolean AllTranscripts;
    public final String BamFile;
    public final File RefGenomeFile;
    public IndexedFastaSequenceFile RefFastaSeqFile;
    public final int ReadCountLimit;
    public final int LongFragmentLimit;
    public final boolean KeepDuplicates;
    public final boolean MarkDuplicates;

    public final boolean WriteExonData;
    public final boolean WriteReadData;
    public final boolean WriteFragmentLengths;
    public final boolean GeneStatsOnly;
    public final int FragmentLengthMinCount;
    public final boolean FragmentLengthsByGene;

    public final List<String> SpecificTransIds;

    public static final int DEFAULT_MAX_READ_COUNT = 100000;
    public static final int DEFAULT_MAX_FRAGMENT_SIZE = 1000;

    public static final int GENE_FRAGMENT_BUFFER = 1000; // width around a gene within which to search for reads

    private static final Logger LOGGER = LogManager.getLogger(RnaExpConfig.class);

    public RnaExpConfig(final CommandLine cmd)
    {
        RestrictedGeneIds = Lists.newArrayList();

        if(cmd.hasOption(GENE_ID_FILE))
        {
            loadGeneIdsFile(cmd.getOptionValue(GENE_ID_FILE));
        }

        AllTranscripts = cmd.hasOption(ALL_TRANSCRIPTS);
        OutputDir = cmd.getOptionValue(DATA_OUTPUT_DIR);

        BamFile = cmd.getOptionValue(BAM_FILE);
        GcBiasFile = cmd.getOptionValue(GC_BIAS_FILE, "");

        final String refGenomeFilename = cmd.getOptionValue(REF_GENOME);
        RefGenomeFile = new File(refGenomeFilename);

        RefFastaSeqFile = null;

        try
        {
            LOGGER.debug("loading indexed fasta reference file");
            RefFastaSeqFile = new IndexedFastaSequenceFile(new File(refGenomeFilename));
        }
        catch (IOException e)
        {
            LOGGER.error("Reference file loading failed: {}", e.toString());
        }

        ReadCountLimit = Integer.parseInt(cmd.getOptionValue(READ_COUNT_LIMIT, "0"));
        LongFragmentLimit = Integer.parseInt(cmd.getOptionValue(LONG_FRAGMENT_LIMIT, String.valueOf(DEFAULT_MAX_FRAGMENT_SIZE)));
        KeepDuplicates = cmd.hasOption(KEEP_DUPLICATES);
        MarkDuplicates = cmd.hasOption(MARK_DUPLICATES);
        FragmentLengthMinCount = Integer.parseInt(cmd.getOptionValue(FRAG_LENGTH_MIN_COUNT, "0"));
        FragmentLengthsByGene = cmd.hasOption(FRAG_LENGTHS_BY_GENE);

        WriteExonData = cmd.hasOption(WRITE_EXON_DATA);
        WriteFragmentLengths = cmd.hasOption(WRITE_FRAGMENT_LENGTHS);
        WriteReadData = cmd.hasOption(WRITE_READ_DATA);
        GeneStatsOnly = cmd.hasOption(GENE_STATS_ONLY);

        SpecificTransIds = cmd.hasOption(SPECIFIC_TRANS_IDS) ?
                Arrays.stream(cmd.getOptionValue(SPECIFIC_TRANS_IDS).split(";")).collect(Collectors.toList())
                : Lists.newArrayList();
    }

    public RnaExpConfig()
    {
        RestrictedGeneIds = Lists.newArrayList();
        OutputDir = "";
        BamFile = "";
        RefGenomeFile = null;
        RefFastaSeqFile = null;
        AllTranscripts = true;
        ReadCountLimit = DEFAULT_MAX_READ_COUNT;
        GcBiasFile = "";
        LongFragmentLimit = DEFAULT_MAX_FRAGMENT_SIZE;
        KeepDuplicates = false;
        MarkDuplicates = false;

        WriteExonData = false;
        WriteReadData = false;
        WriteFragmentLengths = false;
        GeneStatsOnly = false;
        FragmentLengthsByGene = false;
        FragmentLengthMinCount = 0;
        SpecificTransIds = Lists.newArrayList();
    }

    public static boolean checkValid(final CommandLine cmd)
    {
        if(!cmd.hasOption(SAMPLE) || !cmd.hasOption(GENE_TRANSCRIPTS_DIR))
            return false;

        return cmd.hasOption(REF_GENOME) && cmd.hasOption(BAM_FILE);
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample ID");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Path to Ensembl data cache");

        options.addOption(ALL_TRANSCRIPTS, false, "Check all transcripts, not just canonical");
        options.addOption(GENE_ID_FILE, true, "Optional CSV file of genes to analyse");
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
        options.addOption(GENE_STATS_ONLY, false, "Skip all processing except gene summary data");
        options.addOption(WRITE_FRAGMENT_LENGTHS, false, "Write intronic fragment lengths to log");

        options.addOption(SPECIFIC_TRANS_IDS, true, "List of transcripts separated by ';'");

        return options;
    }

    private void loadGeneIdsFile(final String filename)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            LOGGER.warn("invalid gene ID file({})", filename);
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

            RestrictedGeneIds.addAll(fileContents.stream().map(x -> x.split(",")[0]).collect(Collectors.toList()));

            LOGGER.info("file({}) analysing {} specific genes", filename, RestrictedGeneIds.size());
        }
        catch (IOException e)
        {
            LOGGER.warn("failed to load gene ID file({}): {}", filename, e.toString());
        }
    }

}
