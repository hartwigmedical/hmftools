package com.hartwig.hmftools.redux.utils;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.redux.write.BamMerger;

import org.jetbrains.annotations.NotNull;

public class MergeBams
{
    private final BamMerger mBamMerger;

    private static final String INPUT_BAMS = "input_bams";
    private static final String OUTPUT_BAM = "output_bam";

    public MergeBams(final ConfigBuilder configBuilder)
    {
        String outputBam = configBuilder.getValue(OUTPUT_BAM);

        List<String> inputBams = Arrays.stream(configBuilder.getValue(INPUT_BAMS).split(",", -1)).collect(Collectors.toList());

        String refGenome = configBuilder.getValue(REF_GENOME);

        mBamMerger = new BamMerger(outputBam, inputBams, refGenome);
     }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        RD_LOGGER.info("starting BAM merge");

        mBamMerger.merge();

        RD_LOGGER.info("BAM merge complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(INPUT_BAMS, true, "List of input BAMs to be merged, separated ','");
        configBuilder.addConfigItem(OUTPUT_BAM, true, "Output BAM filename");

        RefGenomeSource.addRefGenomeFile(configBuilder, true);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        MergeBams mergeBams = new MergeBams(configBuilder);
        mergeBams.run();
    }
}
