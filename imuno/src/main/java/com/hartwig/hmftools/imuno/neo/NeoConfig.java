package com.hartwig.hmftools.imuno.neo;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.LOG_DEBUG;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.loadGeneIdsFile;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.loadSampleIdsFile;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class NeoConfig
{
    public final List<String> SampleIds;

    public final RefGenomeInterface RefGenome;

    public final List<String> RestrictedGeneIds;

    public final int RequiredAminoAcids;
    public final boolean WriteTransData;
    public final String OutputDir;

    public static final String SAMPLE = "sample";

    public static final String GENE_TRANSCRIPTS_DIR = "gene_transcripts_dir";
    public static final String SV_FUSION_FILE = "sv_fusion_file";
    public static final String GENE_ID_FILE = "gene_id_file";
    public static final String REF_GENOME = "ref_genome";
    public static final String REQ_AMINO_ACIDS = "req_amino_acids";
    public static final String WRITE_TRANS_DATA = "write_trans_data";

    public static final int DEFAULT_AMINO_ACID_REF_COUNT = 18;

    public NeoConfig(final CommandLine cmd)
    {
        SampleIds = Lists.newArrayList();

        final String sampleIdConfig = cmd.getOptionValue(SAMPLE);

        if(sampleIdConfig.contains(".csv"))
        {
            loadSampleIdsFile(sampleIdConfig, SampleIds);
        }
        else if(sampleIdConfig.contains(";"))
        {
            Arrays.stream(sampleIdConfig.split(";")).forEach(x -> SampleIds.add(x));
        }
        else
        {
            SampleIds.add(sampleIdConfig);
        }

        final String refGenomeFilename = cmd.getOptionValue(REF_GENOME);
        RefGenome = loadRefGenome(refGenomeFilename);

        RestrictedGeneIds = Lists.newArrayList();
        OutputDir = parseOutputDir(cmd);

        RequiredAminoAcids = Integer.parseInt(cmd.getOptionValue(REQ_AMINO_ACIDS, String.valueOf(DEFAULT_AMINO_ACID_REF_COUNT)));
        WriteTransData = cmd.hasOption(WRITE_TRANS_DATA);

        if(cmd.hasOption(GENE_ID_FILE))
        {
            loadGeneIdsFile(cmd.getOptionValue(GENE_ID_FILE), RestrictedGeneIds);
        }
    }

    public boolean isMultiSample() { return SampleIds.size() > 1; }

    public NeoConfig(
            final List<String> sampleIds, final RefGenomeInterface refGenome, final List<String> restrictedGeneIds,
            final int requiredAminoAcids)
    {
        SampleIds = sampleIds;
        RefGenome = refGenome;
        RestrictedGeneIds = restrictedGeneIds;
        RequiredAminoAcids = requiredAminoAcids;
        OutputDir = "";
        WriteTransData = false;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE, true, "Sample - Id(s) separated by ';' or CSV file");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Ensembl data cache directory");
        options.addOption(GENE_ID_FILE, true, "Restrict to specific genes");
        options.addOption(REF_GENOME, true, "Ref genome");
        options.addOption(SV_FUSION_FILE, true, "SV fusion file (single sample or cohort)");
        options.addOption(WRITE_TRANS_DATA, false, "Write transcript data for each neo-epitope");
        options.addOption(REQ_AMINO_ACIDS, true, "Number of amino acids in neo-epitopes (default: 18)");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        DatabaseAccess.addDatabaseCmdLineArgs(options);
    }



}
