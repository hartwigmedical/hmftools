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

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
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

    public final boolean RewriteVcf;
    public final boolean WriteFitVariants;
    public final String FitVariantsFile;

    public static final String INPUT_VCF = "input_vcf";
    public static final String WRITE_TINC_VCF = "write_tinc_vcf";
    public static final String FIT_VARIANTS_FILE = "fit_variants_file";
    public static final String WRITE_FIT_VARIANTS = "write_fit_variants";

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

        FitVariantsFile = configBuilder.getValue(FIT_VARIANTS_FILE);
        WriteFitVariants = !configBuilder.hasValue(FIT_VARIANTS_FILE) && configBuilder.hasFlag(WRITE_FIT_VARIANTS);
        RewriteVcf = configBuilder.hasFlag(WRITE_TINC_VCF);
    }

    public TincConfig(
            final RefGenomeVersion refGenVersion, final String inputVcf, final String tumorId, final String referenceId,
            final String ponFilename, final String gnomadDirectory, final String gnomadFile, final int threads, boolean overwriteVcf)
    {
        RefGenVersion = refGenVersion;
        InputVcf = inputVcf;
        TumorId = tumorId;
        ReferenceId = referenceId;
        PonFilename = ponFilename;
        GnomadDirectory = gnomadDirectory;
        GnomadFile = gnomadFile;
        Threads = threads;
        WriteFitVariants = false;
        FitVariantsFile = null;
        RewriteVcf = overwriteVcf;
    }

    public static void registerFullConfig(final ConfigBuilder configBuilder)
    {
        registerConfig(configBuilder);

        addRefGenomeConfig(configBuilder, false);

        addOutputDir(configBuilder);
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(TUMOR, true, TUMOR_IDS_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_IDS_DESC);
        configBuilder.addPath(INPUT_VCF, true, "Input unfiltered VCF");
        configBuilder.addPath(FIT_VARIANTS_FILE, false, "Cache file for fit variants");
        configBuilder.addFlag(WRITE_FIT_VARIANTS, "Write cache file for fit variants");
        configBuilder.addFlag(WRITE_TINC_VCF, "Write a new VCF with TINC changes if detected");
        configBuilder.addPath(PON_FILE, false, "PON entries");
        GnomadCache.addConfig(configBuilder);
    }
}
