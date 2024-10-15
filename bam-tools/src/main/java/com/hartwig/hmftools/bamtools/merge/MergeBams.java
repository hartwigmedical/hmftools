package com.hartwig.hmftools.bamtools.merge;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bamops.BamMerger;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class MergeBams
{
    private final BamMerger mBamMerger;

    private static final String INPUT_BAMS = "input_bams";
    private static final String OUTPUT_BAM = "output_bam";
    private static final String KEEP_INTERIM_BAMS = "keep_interim_bams";

    public MergeBams(final ConfigBuilder configBuilder)
    {
        String outputBam = configBuilder.getValue(OUTPUT_BAM);

        List<String> inputBams = Arrays.stream(configBuilder.getValue(INPUT_BAMS).split(",", -1)).collect(Collectors.toList());

        // check input BAMs exist
        boolean hasMissing = false;

        for(String inputBam : inputBams)
        {
            if(!Files.exists(Paths.get(inputBam)))
            {
                BT_LOGGER.error("missing input BAM: {}", inputBam);
                hasMissing = true;
            }
        }

        if(hasMissing)
            System.exit(1);

        String refGenome = configBuilder.getValue(REF_GENOME);

        String bamToolPath = configBuilder.getValue(BAMTOOL_PATH);
        boolean keepInterimBams = configBuilder.hasFlag(KEEP_INTERIM_BAMS);

        int threads = parseThreads(configBuilder);
        mBamMerger = new BamMerger(outputBam, inputBams, refGenome, bamToolPath, threads, keepInterimBams);
     }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        BT_LOGGER.info("starting BAM merge of {} BAMs", mBamMerger.inputBamCount());

        mBamMerger.merge();

        BT_LOGGER.info("BAM merge complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(INPUT_BAMS, true, "List of input BAMs to be merged, separated ','");
        configBuilder.addConfigItem(OUTPUT_BAM, true, "Output BAM filename");
        configBuilder.addFlag(KEEP_INTERIM_BAMS, "Do no delete per-thread BAMs");

        RefGenomeSource.addRefGenomeFile(configBuilder, true);
        BamToolName.addConfig(configBuilder);

        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        MergeBams mergeBams = new MergeBams(configBuilder);
        mergeBams.run();
    }
}
