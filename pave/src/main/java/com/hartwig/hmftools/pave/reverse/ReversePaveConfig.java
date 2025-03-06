package com.hartwig.hmftools.pave.reverse;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ReversePaveConfig
{
    public final String mVcfFile;
    public RefGenomeSource mRefGenome;
    public final RefGenomeVersion mRefGenVersion;
    public final EnsemblDataCache mEnsemblCache;

    public static final String VCF_FILE = "vcf_file";

    public static final Logger RPV_LOGGER = LogManager.getLogger(ReversePaveConfig.class);

    public ReversePaveConfig(final ConfigBuilder configBuilder)
    {
        mVcfFile = configBuilder.getValue(VCF_FILE);
        try
        {
            final String refGenomeFile = configBuilder.getValue(REF_GENOME);
            IndexedFastaSequenceFile refGenome = new IndexedFastaSequenceFile(new File(refGenomeFile));
            mRefGenome = new RefGenomeSource(refGenome);
        }
        catch(IOException e)
        {
            RPV_LOGGER.error("failed to load ref genome: {}", e.toString());
        }
        mRefGenVersion = RefGenomeVersion.from(configBuilder);
        mEnsemblCache = new EnsemblDataCache(configBuilder);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(VCF_FILE, true, "VCF input file");

        addRefGenomeConfig(configBuilder, true);
        addEnsemblDir(configBuilder, true);
        addThreadOptions(configBuilder);

        addOutputDir(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
