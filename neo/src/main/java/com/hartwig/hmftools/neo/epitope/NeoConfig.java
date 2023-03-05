package com.hartwig.hmftools.neo.epitope;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class NeoConfig
{
    public final List<String> Samples;

    public final RefGenomeInterface RefGenome;
    public final RefGenomeVersion RefGenVersion;

    public final List<String> RestrictedGeneIds;
    public final List<String> CommonHlaTypes;
    public final String MutationsFile;

    public final int RequiredAminoAcids;
    public final boolean WriteTransData;
    public final String SvFusionsDir;
    public final String SomaticVcf;
    public final String OutputDir;
    public final int Threads;

    public static final String SAMPLE = "sample";
    public static final String CANCER_TYPE = "cancer_type";
    public static final String PEPTIDE_LENGTHS = "peptide_lengths";
    public static final String PEPTIDE_FLANKS = "peptide_flanks";

    public static final String MUTATIONS_FILE = "mutations_file";
    public static final String SV_FUSION_DATA_DIR = "sv_fusion_data_dir";
    public static final String SOMATIC_VCF = "somatic_vcf";
    public static final String GENE_ID_FILE = "gene_id_file";
    public static final String REQ_AMINO_ACIDS = "req_amino_acids";

    public static final String WRITE_TRANS_DATA = "write_trans_data";

    public static final int DEFAULT_AMINO_ACID_REF_COUNT = 18;

    public NeoConfig(final CommandLine cmd)
    {
        Samples = Lists.newArrayList();

        if(cmd.hasOption(SAMPLE))
        {
            Samples.add(cmd.getOptionValue(SAMPLE));
        }
        else
        {
            Samples.addAll(loadSampleIdsFile(cmd));
        }

        final String refGenomeFilename = cmd.getOptionValue(REF_GENOME);
        RefGenome = loadRefGenome(refGenomeFilename);
        RefGenVersion = RefGenomeVersion.from(cmd);

        SvFusionsDir = cmd.getOptionValue(SV_FUSION_DATA_DIR);
        SomaticVcf = cmd.getOptionValue(SOMATIC_VCF);
        OutputDir = parseOutputDir(cmd);

        RequiredAminoAcids = Integer.parseInt(cmd.getOptionValue(REQ_AMINO_ACIDS, String.valueOf(DEFAULT_AMINO_ACID_REF_COUNT)));

        CommonHlaTypes = Lists.newArrayList();

        MutationsFile = cmd.getOptionValue(MUTATIONS_FILE);

        WriteTransData = cmd.hasOption(WRITE_TRANS_DATA);

        Threads = parseThreads(cmd);

        RestrictedGeneIds = Lists.newArrayList();
        if(cmd.hasOption(GENE_ID_FILE))
            RestrictedGeneIds.addAll(loadGeneIdsFile(cmd.getOptionValue(GENE_ID_FILE)));
    }

    public boolean isMultiSample() { return Samples.size() > 1; }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE, true, "Sample - Id(s) separated by ';' or CSV file");
        options.addOption(CANCER_TYPE, true, "Tumor cancer type (optional) - to retrieve cancer median TPM");
        options.addOption(MUTATIONS_FILE, true, "File with a list of point mutations");
        options.addOption(PEPTIDE_LENGTHS, true, "Peptide length min-max, separated by '-', eg 8-12");
        options.addOption(PEPTIDE_FLANKS, true, "Peptide flanking amino acids");
        EnsemblDataCache.addEnsemblDir(options);
        options.addOption(GENE_ID_FILE, true, "Restrict to specific genes");
        addRefGenomeConfig(options);
        options.addOption(SV_FUSION_DATA_DIR, true, "Linx neoepitope directory");
        options.addOption(SOMATIC_VCF, true, "Purple somatic VCF (use '*') as required");
        options.addOption(WRITE_TRANS_DATA, false, "Write transcript data for each neo-epitope");
        options.addOption(REQ_AMINO_ACIDS, true, "Number of amino acids in neo-epitopes (default: 18)");
        addLoggingOptions(options);
        addOutputOptions(options);
        addThreadOptions(options);
        DatabaseAccess.addDatabaseCmdLineArgs(options);
    }
}
