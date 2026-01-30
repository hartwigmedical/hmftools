package com.hartwig.hmftools.geneutils.mapping;

import static com.hartwig.hmftools.common.blastn.BlastnRunner.registerBlastn;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

public class SeqTestConfig
{
    public final String InputFile;
    public final String OutputDir;
    public final String OutputId;

    public final int Threads;

    public final RefGenomeSource RefGenome;
    public final RefGenomeVersion RefGenomeVersion;

    private static final String INPUT_FILE = "input_file";

    public SeqTestConfig(final ConfigBuilder configBuilder)
    {
        String refGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenome = loadRefGenome(refGenomeFile);
        RefGenomeVersion = deriveRefGenomeVersion(RefGenome);

        InputFile = configBuilder.getValue(INPUT_FILE);
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        Threads = parseThreads(configBuilder);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(INPUT_FILE, true, "Input regions and sequences");

        registerBlastn(configBuilder, false);

        addOutputOptions(configBuilder, false);
        ConfigUtils.addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
        addRefGenomeFile(configBuilder, true);

    }
}
