package com.hartwig.hmftools.sage.tinc;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.variant.pon.GnomadCache.GNOMAD_FREQUENCY_DIR;
import static com.hartwig.hmftools.common.variant.pon.GnomadCache.GNOMAD_FREQUENCY_FILE;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_FILE;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.pon.GnomadCache;

public class TincConfig
{
    public final RefGenomeVersion RefGenVersion;

    public final String InputVcf;
    public final String TumorId;
    public final String ReferenceId;

    public final String PonFilename;
    public final String GnomadDirectory;
    public final String GnomadFile;
    public final int Threads;

    public static final String INPUT_VCF = "input_vcf";

    public TincConfig(final ConfigBuilder configBuilder)
    {
        InputVcf = configBuilder.getValue(INPUT_VCF);
        TumorId = configBuilder.getValue(TUMOR);
        ReferenceId = configBuilder.getValue(REFERENCE);

        GnomadDirectory= configBuilder.getValue(GNOMAD_FREQUENCY_DIR);
        GnomadFile = configBuilder.getValue(GNOMAD_FREQUENCY_FILE);
        PonFilename = configBuilder.getValue(PON_FILE);
        RefGenVersion = RefGenomeVersion.from(configBuilder);

        Threads = parseThreads(configBuilder);
    }

    public static void registerFullConfig(final ConfigBuilder configBuilder)
    {
        registerConfig(configBuilder);

        addRefGenomeVersion(configBuilder);

        addOutputDir(configBuilder);
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(TUMOR, true, TUMOR_IDS_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_IDS_DESC);
        configBuilder.addPath(INPUT_VCF, true, "Input unfiltered VCF");
        configBuilder.addPath(PON_FILE, false, "PON entries");
        GnomadCache.addConfig(configBuilder);
    }
}
