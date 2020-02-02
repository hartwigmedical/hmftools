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

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RnaExpConfig
{
    // config items
    public static final String GENE_TRANSCRIPTS_DIR = "gene_transcripts_dir";
    public static final String SAMPLE = "sample";

    private static final String GENE_ID_FILE = "gene_id_file";
    private static final String ALL_TRANSCRIPTS = "all_transcripts";
    private static final String WRITE_MATCHED_READS = "write_matched_reads";

    public static final String REF_GENOME = "ref_genome";
    public static final String BAM_FILE = "bam_file";
    public static final String GC_BIAS_FILE = "gcbias_file";
    public static final String READ_COUNT_LIMIT = "read_count_limit";

    public static final String SPECIFIC_TRANS_IDS = "specific_trans";

    public final List<String> RestrictedGeneIds;
    public final String OutputDir;
    public final String GcBiasFile;
    public final boolean AllTranscripts;
    public final String BamFile;
    public final File RefGenomeFile;
    public IndexedFastaSequenceFile RefFastaSeqFile;
    public final int ReadCountLimit;

    public final List<String> SpecificTransIds;

    public static final int MAX_READ_COUNT = 100000;

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

        ReadCountLimit = Integer.parseInt(cmd.getOptionValue(READ_COUNT_LIMIT, String.valueOf(MAX_READ_COUNT)));

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
        ReadCountLimit = MAX_READ_COUNT;
        GcBiasFile = "";
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
        options.addOption(BAM_FILE, true, "RNA BAM file location");

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
