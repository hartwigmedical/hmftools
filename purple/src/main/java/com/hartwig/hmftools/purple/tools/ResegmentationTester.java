package com.hartwig.hmftools.purple.tools;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.purple.GermlineStatus.DIPLOID;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.purple.PurpleSegment;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.segment.Resegmentation;

import org.jetbrains.annotations.NotNull;

public class ResegmentationTester
{
    private final String mSampleId;
    private final String mPurpleDir;
    private final int mThreads;

    private final BufferedWriter mWriter;

    public ResegmentationTester(final ConfigBuilder configBuilder)
    {
        mSampleId = configBuilder.getValue(SAMPLE);
        mPurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);

        mThreads = parseThreads(configBuilder);

        mWriter = null; // initialiseWriter(parseOutputDir(configBuilder));
    }

    public void run()
    {
        PPL_LOGGER.info("running Purple resegmentation tester on sample({})", mSampleId);

        String samplePurpleDir = convertWildcardSamplePath(mPurpleDir, mSampleId);

        try
        {
            List<PurpleSegment> segments = PurpleSegment.read(PurpleSegment.generateFilename(mPurpleDir, mSampleId));
            List<ObservedRegion> observedRegions = segments.stream().map(x -> ObservedRegion.fromSegment(x)).collect(Collectors.toList());

            Resegmentation resegmentation = new Resegmentation(observedRegions);
            resegmentation.processSegments();
        }
        catch(Exception e)
        {
            PPL_LOGGER.error("failed to read Purple files: {}", e.toString());
            System.exit(1);
        }

        PPL_LOGGER.info("Purple resegmentation tester complete");
    }


    private BufferedWriter initialiseWriter(final String outputDir)
    {
        try
        {
            String fileName = outputDir + "resegmentation" + TSV_EXTENSION;

            /*
            String fileType = fileType();

            if(!fileType.isEmpty())
                fileName += "." + fileType;
            */

            // fileName += TSV_EXTENSION;

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("SampleId\tChromosome\tRegionStart\tRegionEnd\tGermlineStatus\tBafCount\tDepthWindows");
            writer.write("\tObsTumorRatio\tObsNormalRatio\tObsUnnormalisedNormalRatio\tRegionRefNormalisedCN");


            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to initialise output file output: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeSegmentData(final String sampleId, final ObservedRegion region)
    {
        try
        {
            mWriter.write(String.format("%s\t%s\t%d\t%d\t%s",
                    sampleId, region.chromosome(), region.start(), region.end(), region.germlineStatus()));

            mWriter.write(String.format("\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f",
                    region.bafCount(), region.depthWindowCount(), region.observedTumorRatio(),
                    region.observedNormalRatio(), region.unnormalisedObservedNormalRatio(), region.refNormalisedCopyNumber()));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to write germline gene deletion file output: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        // addSampleIdFile(configBuilder, true);

        configBuilder.addPath(PURPLE_DIR_CFG, true, PURPLE_DIR_CFG);
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);

        addLoggingOptions(configBuilder);
        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        ResegmentationTester resegmentationTester = new ResegmentationTester(configBuilder);
        resegmentationTester.run();
    }
}
