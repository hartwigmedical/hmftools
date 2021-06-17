package com.hartwig.hmftools.neo;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.NeoCommon.LOG_DEBUG;
import static com.hartwig.hmftools.neo.NeoCommon.loadGeneIdsFile;
import static com.hartwig.hmftools.neo.NeoCommon.loadSampleDataFile;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class NeoConfig
{
    public final List<SampleData> Samples;
    public final int[] PeptideLengths;
    public final int PeptideFlanks;

    public final RefGenomeInterface RefGenome;

    public final List<String> RestrictedGeneIds;
    public final List<String> CommonHlaTypes;
    public final String MutationsFile;

    public final int RequiredAminoAcids;
    public final boolean WriteTransData;
    public final boolean WriteCohortFiles;
    public final String SvFusionsDir;
    public final String OutputDir;
    public final int Threads;

    public static final String SAMPLE = "sample";
    public static final String CANCER_TYPE = "cancer_type";
    public static final String PEPTIDE_LENGTHS = "peptide_lengths";
    public static final String PEPTIDE_FLANKS = "peptide_flanks";

    public static final String HLA_TYPES_FILE = "hla_types_file";
    public static final String MUTATIONS_FILE = "mutations_file";

    public static final String GENE_TRANSCRIPTS_DIR = "gene_transcripts_dir";
    public static final String SV_FUSION_DATA_DIR = "sv_fusion_data_dir";
    public static final String GENE_ID_FILE = "gene_id_file";
    public static final String REF_GENOME = "ref_genome";
    public static final String CANCER_TPM_FILE = "cancer_tpm_file";
    public static final String REQ_AMINO_ACIDS = "req_amino_acids";
    public static final String WRITE_TRANS_DATA = "write_trans_data";
    public static final String WRITE_COHORT_FILES = "write_cohort_files";
    public static final String THREADS = "threads";

    public static final String NEO_EPITOPE_FILE_ID = ".imu.neo_epitopes.csv";
    public static final String HLA_PEPTIDE_FILE_ID = ".imu.hla_peptides.csv";

    public static final int DEFAULT_AMINO_ACID_REF_COUNT = 18;

    public NeoConfig(final CommandLine cmd)
    {
        Samples = Lists.newArrayList();

        final String sampleIdConfig = cmd.getOptionValue(SAMPLE);

        if(sampleIdConfig.contains(".csv"))
        {
            loadSampleDataFile(sampleIdConfig, Samples);
        }
        else
        {
            if(sampleIdConfig.contains(DELIMITER))
            {
                SampleData sample = SampleData.fromCsv(sampleIdConfig);

                if(sample == null)
                {
                    NE_LOGGER.error("invalid sample data: {}", sampleIdConfig);
                }
                else
                {
                    Samples.add(sample);
                }
            }
            else
            {
                Samples.add(new SampleData(sampleIdConfig, "", Lists.newArrayList()));
            }
        }

        PeptideLengths = new int[SE_PAIR];

        if(cmd.hasOption(PEPTIDE_LENGTHS))
        {
            String[] lengths = cmd.getOptionValue(PEPTIDE_LENGTHS).split("-");

            if(lengths.length == 1)
            {
                PeptideLengths[SE_START] = PeptideLengths[SE_END] = Integer.parseInt(lengths[0]);
            }
            else if(lengths.length == 2)
            {
                PeptideLengths[SE_START] = Integer.parseInt(lengths[SE_START]);
                PeptideLengths[SE_END] = Integer.parseInt(lengths[SE_END]);
            }
            else
            {
                NE_LOGGER.error("invalid peptide lengths: {}", cmd.getOptionValue(PEPTIDE_LENGTHS));
            }
        }

        PeptideFlanks = Integer.parseInt(cmd.getOptionValue(PEPTIDE_FLANKS, "0"));

        final String refGenomeFilename = cmd.getOptionValue(REF_GENOME);
        RefGenome = loadRefGenome(refGenomeFilename);

        SvFusionsDir = cmd.getOptionValue(SV_FUSION_DATA_DIR);
        OutputDir = parseOutputDir(cmd);

        RequiredAminoAcids = Integer.parseInt(cmd.getOptionValue(REQ_AMINO_ACIDS, String.valueOf(DEFAULT_AMINO_ACID_REF_COUNT)));

        CommonHlaTypes = Lists.newArrayList();

        if(cmd.hasOption(HLA_TYPES_FILE))
            loadCommonHlaTypes(cmd.getOptionValue(HLA_TYPES_FILE));

        MutationsFile = cmd.getOptionValue(MUTATIONS_FILE);

        WriteTransData = cmd.hasOption(WRITE_TRANS_DATA);
        WriteCohortFiles = cmd.hasOption(WRITE_COHORT_FILES);

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));

        RestrictedGeneIds = Lists.newArrayList();
        if(cmd.hasOption(GENE_ID_FILE))
            loadGeneIdsFile(cmd.getOptionValue(GENE_ID_FILE), RestrictedGeneIds);
    }

    public boolean isMultiSample() { return Samples.size() > 1; }

    public NeoConfig(
            final List<SampleData> samples, final RefGenomeInterface refGenome, final List<String> restrictedGeneIds,
            final int requiredAminoAcids)
    {
        Samples = samples;
        PeptideLengths = new int[SE_PAIR];
        PeptideFlanks = 0;
        RefGenome = refGenome;
        RestrictedGeneIds = restrictedGeneIds;
        RequiredAminoAcids = requiredAminoAcids;
        CommonHlaTypes = Lists.newArrayList();
        MutationsFile = null;
        SvFusionsDir = "";
        OutputDir = "";
        WriteTransData = false;
        WriteCohortFiles = false;
        Threads = 0;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE, true, "Sample - Id(s) separated by ';' or CSV file");
        options.addOption(CANCER_TYPE, true, "Tumor cancer type (optional) - to retrieve cancer median TPM");
        options.addOption(HLA_TYPES_FILE, true, "File with a list of HLA types");
        options.addOption(MUTATIONS_FILE, true, "File with a list of point mutations");
        options.addOption(PEPTIDE_LENGTHS, true, "Peptide length min-max, separated by '-', eg 8-12");
        options.addOption(PEPTIDE_FLANKS, true, "Peptide flanking amino acids");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Ensembl data cache directory");
        options.addOption(GENE_ID_FILE, true, "Restrict to specific genes");
        options.addOption(REF_GENOME, true, "Ref genome");
        options.addOption(SV_FUSION_DATA_DIR, true, "SV fusion file (single sample or cohort)");
        options.addOption(CANCER_TPM_FILE, true, "TPM per cancer type and pan-cancer");
        options.addOption(WRITE_TRANS_DATA, false, "Write transcript data for each neo-epitope");
        options.addOption(WRITE_COHORT_FILES, false, "Write cohort files for multiple samples");
        options.addOption(REQ_AMINO_ACIDS, true, "Number of amino acids in neo-epitopes (default: 18)");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(THREADS, true, "Thread count");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        DatabaseAccess.addDatabaseCmdLineArgs(options);
    }

    public void loadCommonHlaTypes(final String filename)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.warn("invalid HLA types file({})", filename);
            return;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return;

            if(fileContents.get(0).contains("HlaType"))
                fileContents.remove(0);

            CommonHlaTypes.addAll(fileContents);
        }
        catch (IOException e)
        {
            NE_LOGGER.warn("failed to load common HLA types file({}): {}", filename, e.toString());
        }
    }

}
