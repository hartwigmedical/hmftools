package com.hartwig.hmftools.neo.epitope;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class NeoConfig
{
    public final List<String> Samples;

    public final RefGenomeInterface RefGenome;
    public final RefGenomeVersion RefGenVersion;

    public final List<String> RestrictedGeneIds;

    public final int RequiredAminoAcids;
    public final boolean WriteTransData;
    public final String LinxDir;
    public final String SomaticVcf;
    public final String OutputDir;
    public final int Threads;

    public static final String CANCER_TYPE = "cancer_type";

    public static final String SOMATIC_VCF = "somatic_vcf";
    public static final String REQ_AMINO_ACIDS = "req_amino_acids";

    public static final String WRITE_TRANS_DATA = "write_trans_data";

    public static final int DEFAULT_AMINO_ACID_REF_COUNT = 15;

    public NeoConfig(final ConfigBuilder configBuilder)
    {
        Samples = Lists.newArrayList();

        if(configBuilder.hasValue(SAMPLE))
        {
            Samples.add(configBuilder.getValue(SAMPLE));
        }
        else
        {
            Samples.addAll(loadSampleIdsFile(configBuilder));
        }

        final String refGenomeFilename = configBuilder.getValue(REF_GENOME);
        RefGenome = loadRefGenome(refGenomeFilename);
        RefGenVersion = RefGenomeVersion.from(configBuilder);

        LinxDir = configBuilder.getValue(LINX_DIR_CFG);
        SomaticVcf = configBuilder.getValue(SOMATIC_VCF);
        OutputDir = parseOutputDir(configBuilder);

        RequiredAminoAcids = configBuilder.getInteger(REQ_AMINO_ACIDS);

        WriteTransData = configBuilder.hasFlag(WRITE_TRANS_DATA);

        Threads = parseThreads(configBuilder);

        RestrictedGeneIds = Lists.newArrayList();
        if(configBuilder.hasValue(GENE_ID_FILE))
            RestrictedGeneIds.addAll(loadGeneIdsFile(configBuilder.getValue(GENE_ID_FILE)));
    }

    public static void addCmdLineArgs(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, false, "Sample - Id(s) separated by ';' or CSV file");
        addSampleIdFile(configBuilder, false);
        configBuilder.addConfigItem(CANCER_TYPE, false, "Tumor cancer type (optional) - to retrieve cancer median TPM");
        configBuilder.addPath(LINX_DIR_CFG, true, "Linx neoepitope directory");
        configBuilder.addPath(SOMATIC_VCF, true, "Purple somatic VCF (use '*') as required");
        configBuilder.addFlag(WRITE_TRANS_DATA, "Write transcript data for each neo-epitope");

        configBuilder.addInteger(REQ_AMINO_ACIDS, "Number of amino acids in neo-epitopes", DEFAULT_AMINO_ACID_REF_COUNT);

        configBuilder.addPath(GENE_ID_FILE, false, GENE_ID_FILE_DESC);

        addRefGenomeConfig(configBuilder, true);
        EnsemblDataCache.addEnsemblDir(configBuilder);

        addLoggingOptions(configBuilder);
        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}
