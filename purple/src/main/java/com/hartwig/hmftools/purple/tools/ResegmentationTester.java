package com.hartwig.hmftools.purple.tools;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleSegment;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.reseg.Resegmentation;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ResegmentationTester
{
    private final String mSampleId;
    private final String mPurpleDir;
    private final String mOutputDir;
    private final int mThreads;

    private final BufferedWriter mWriter;

    public ResegmentationTester(final ConfigBuilder configBuilder)
    {
        mSampleId = configBuilder.getValue(SAMPLE);
        mPurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        mOutputDir = parseOutputDir(configBuilder);

        mThreads = parseThreads(configBuilder);

        mWriter = null; // initialiseWriter(parseOutputDir(configBuilder));
    }

    public void run()
    {
        PPL_LOGGER.info("running Purple resegmentation tester on sample({})", mSampleId);

        // String samplePurpleDir = convertWildcardSamplePath(mPurpleDir, mSampleId);

        try
        {
            List<PurpleSegment> segments = PurpleSegment.read(PurpleSegment.generateFilename(mPurpleDir, mSampleId));
            List<ObservedRegion> observedRegions = segments.stream().map(x -> ObservedRegion.fromSegment(x)).collect(Collectors.toList());

            List<ObservedRegion> outputRegions = Resegmentation.run(observedRegions, mThreads);
            List<PurpleSegment> outputSegments = outputRegions.stream().map(x -> x.toSegment()).collect(Collectors.toList());

            String resegmentedFile = mOutputDir + mSampleId + ".purple.resegmented.tsv";

            PurpleSegment.write(resegmentedFile, outputSegments);

            writeSegmentDifferences(segments, outputSegments);

            // Resegmentation resegmentation = new Resegmentation(observedRegions);
            //resegmentation.processSegments();
        }
        catch(Exception e)
        {
            PPL_LOGGER.error("failed to read Purple files: {}", e.toString());
            System.exit(1);
        }

        PPL_LOGGER.info("Purple resegmentation tester complete");
    }

    private void writeSegmentDifferences(final List<PurpleSegment> oldSegments, final List<PurpleSegment> newSegments)
    {
        try
        {
            String filename = mOutputDir + mSampleId + ".purple.segmentation_diffs.tsv";
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("MismatchType").add("Diffs");
            sj.add(PurpleSegment.header());

            writer.write(sj.toString());
            writer.newLine();

            List<PurpleSegment> newSegmentsWorking = Lists.newArrayList(newSegments);

            for(PurpleSegment oldSegment : oldSegments)
            {
                PurpleSegment newSegment = newSegmentsWorking.stream()
                        .filter(x -> segmentsMatched(oldSegment, x, true)).findFirst().orElse(null);

                if(newSegment == null)
                {
                    newSegment = newSegmentsWorking.stream()
                            .filter(x -> segmentsMatched(oldSegment, x, false)).findFirst().orElse(null);
                }

                if(newSegment != null)
                {
                    writeSegmentDifference(writer, oldSegment, newSegment);
                    newSegmentsWorking.remove(newSegment);
                }
                else
                {
                    writeSegmentDifference(writer, oldSegment, null);
                }
            }

            for(PurpleSegment newSegment : newSegmentsWorking)
            {
                writeSegmentDifference(writer, null, newSegment);
            }

            writer.newLine();

            writer.close();
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to write segmentation diff output file: {}", e.toString());
        }
    }

    private static boolean segmentsMatched(final PurpleSegment first, final PurpleSegment second, boolean requireExactPositions)
    {
        if(!first.Chromosome.equals(second.Chromosome))
            return false;

        if(first.start() != second.start())
            return false;

        return !requireExactPositions || first.end() == second.end();
    }

    private static final int MAX_COUNT_DIFF = 10;
    private static final double MAX_RATIO_DIFF = 0.05;
    private static final int MAX_REGION_DIFF = 5000;

    private static boolean hasCountDifference(int value1, int value2)
    {
        return abs(value1 - value2) >= MAX_COUNT_DIFF;
    }

    private static boolean hasRatioDifference(double value1, double value2)
    {
        return abs(value1 - value2) >= MAX_RATIO_DIFF;
    }

    private static void writeSegmentDifference(
            final BufferedWriter writer, @Nullable final PurpleSegment oldSegment, @Nullable final PurpleSegment newSegment)
            throws IOException
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        // log diffs
        StringJoiner diffs = new StringJoiner(ITEM_DELIM);

        String mismatchType = "";
        if(oldSegment != null && newSegment != null)
        {
            // check specific fields for differences
            if(abs(oldSegment.PosEnd - newSegment.PosEnd) >= MAX_REGION_DIFF)
                diffs.add(format("RegionEnd(%d/%d)", oldSegment.PosEnd, newSegment.PosEnd));

            if(oldSegment.Support != newSegment.Support)
                diffs.add(format("Support(%s/%s)", oldSegment.Support, newSegment.Support));

            if(hasCountDifference(oldSegment.BafCount, newSegment.BafCount))
                diffs.add(format("BafCount(%d/%d)", oldSegment.BafCount, newSegment.BafCount));

            if(hasCountDifference(oldSegment.DepthWindowCount, newSegment.DepthWindowCount))
                diffs.add(format("DepthWindowCount(%d/%d)", oldSegment.DepthWindowCount, newSegment.DepthWindowCount));

            if(hasRatioDifference(oldSegment.ObservedTumorRatio, newSegment.ObservedTumorRatio))
                diffs.add(format("ObservedTumorRatio(%.3f/%.3f)", oldSegment.ObservedTumorRatio, newSegment.ObservedTumorRatio));

            if(oldSegment.GermlineState != newSegment.GermlineState)
                diffs.add(format("GermlineStatus(%s/%s)", oldSegment.GermlineState, newSegment.GermlineState));

            if(diffs.length() == 0)
                return;

            mismatchType = "VALUE";
        }
        else if(oldSegment != null)
        {
            mismatchType = "OLD_ONLY";
        }
        else
        {
            mismatchType = "NEW_ONLY";
        }

        sj.add(mismatchType);
        sj.add(diffs.toString());

        // log all other fields
        PurpleSegment refSegment = oldSegment != null ? oldSegment : newSegment;
        sj.add(refSegment.toTsv());

        writer.write(sj.toString());
        writer.newLine();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

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
