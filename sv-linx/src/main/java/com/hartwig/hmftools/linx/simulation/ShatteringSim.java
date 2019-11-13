package com.hartwig.hmftools.linx.simulation;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ShatteringSim
{
    // config
    private final ShatteringConfig mConfig;
    private final String mOutputDir;

    // state and results
    private final List<Segment> mSegments;
    private int mSegmentCount;
    private int mRemainingLinkCount;
    private Random mRandom;
    private int mRunIndex;
    private boolean mValidRun;
    private ShatteringResult mLatestResult;

    private BufferedWriter mResultsWriter;

    // testing
    private final List<Integer> mSpecifiedOrder;

    private static final Logger LOGGER = LogManager.getLogger(ShatteringSim.class);

    public ShatteringSim(final ShatteringConfig config, final String outputDir)
    {
        mConfig = config;
        mOutputDir = outputDir;
        initialiseWriter();

        mSegments = Lists.newArrayList();
        mSegmentCount = 0;
        mRemainingLinkCount = 0;
        mValidRun = true;
        mLatestResult = null;

        mRandom = new Random(123456);
        mSpecifiedOrder = Lists.newArrayList();
    }

    public boolean validRun() { return mValidRun; }

    public void run()
    {
        for(mSegmentCount = mConfig.SegmentCountMin; mSegmentCount <= mConfig.SegmentCountMax; ++mSegmentCount)
        {
            LOGGER.info("run with segmentCount({})", mSegmentCount);

            performRun();

            if(!mValidRun)
                break;
        }

        close();
    }

    private void performRun()
    {
        initialiseState();

        for(int i = 0; i < mConfig.Iterations; ++i)
        {
            runIteration();
            updateResults();
            writeResults();
            ++mRunIndex;

            if(mConfig.Iterations > 10000 && (i % 10000) == 1)
            {
                LOGGER.info("run index {}", i);
            }

            if(!mValidRun)
                break;
        }

    }

    private void initialiseState()
    {
        mRunIndex = 0;
        mSegments.clear();
        mValidRun = true;

        // create a set of segments with 2 unconnected ends, and 2 bounding segments with a single exposed end
        mSegments.add(new Segment(0, false, true));

        for(int i = 0; i < mSegmentCount; ++i)
        {
            mSegments.add(new Segment(i+1, true, true));
        }

        mSegments.add(new Segment(mSegments.size(), true, false));
    }

    private void clearRunState()
    {
        int breakendCount = mSegmentCount * 2 + 2; // each segment has 2 ends plus the start and end bounding ends
        mRemainingLinkCount = breakendCount;

        mSegments.stream().forEach(x -> x.clearLinks());
    }

    private void runIteration()
    {
        clearRunState();

        boolean foundLink = true;

        while(mRemainingLinkCount > 0 && foundLink)
        {
            // randomly find the next 2 ends to connect
            int[] nextIndices = getNextIndexPair(mRemainingLinkCount);
            int nextIndex1 = nextIndices[0];
            int nextIndex2 = nextIndices[1];

            Segment nextSegment1 = null;
            boolean seg1LinkOnStart = false;
            Segment nextSegment2 = null;
            boolean seg2LinkOnStart = false;

            foundLink = false;

            // use these randomly selected indices to find the next segments to be linked
            int unlinkedIndex = 0;
            for(int i = 0; i < mSegments.size(); ++i)
            {
                Segment segment = mSegments.get(i);

                if (segment.fullyLinked())
                    continue;

                for (int be = SE_START; be <= SE_END; ++be)
                {
                    boolean isStart = isStart(be);

                    if (!segment.isLinkOpen(isStart))
                        continue;

                    if (nextSegment1 == null && unlinkedIndex >= nextIndex1)
                    {
                        nextSegment1 = segment;
                        seg1LinkOnStart = isStart;
                    }
                    else if (nextSegment2 == null && unlinkedIndex >= nextIndex2)
                    {
                        nextSegment2 = segment;
                        seg2LinkOnStart = isStart;
                    }

                    ++unlinkedIndex;
                }

                if (nextSegment1 != null && nextSegment2 != null)
                {
                    foundLink = true;
                    nextSegment1.setLink(nextSegment2, seg1LinkOnStart);
                    nextSegment2.setLink(nextSegment1, seg2LinkOnStart);
                    mRemainingLinkCount -= 2;
                    break;
                }
            }

            if(!foundLink)
            {
                LOGGER.error("no link found: nextIndex({} & {}) remaining({})", nextIndex1, nextIndex2, mRemainingLinkCount);
                mValidRun = false;
                break;
            }

            if(!moreLinksPossible())
                break;
        }
    }

    private boolean moreLinksPossible()
    {
        Segment first = mSegments.get(0);

        if(first.isLinkOpen(false))
            return true;

        Segment last = mSegments.get(mSegments.size()-1);

        if(last.isLinkOpen(true))
            return true;

        // see a path can be traced from first to last
        Segment currentSegment = first;
        boolean nextLinkOnStart = false;

        int iterations = 0;

        while(true)
        {
            // check other side
            Segment nextSegment = currentSegment.getLink(nextLinkOnStart);

            if(nextSegment == null)
                return true;

            if(nextSegment == last)
                return false;

            boolean nextLinkedOnStart = nextSegment.getLink(true) == currentSegment;
            nextLinkOnStart = !nextLinkedOnStart;
            currentSegment = nextSegment;

            ++iterations;

            if(iterations > mSegments.size())
            {
                LOGGER.error("iterations({}) exceeds segCount({})", iterations, mSegments.size());
                mValidRun = false;
                return false;
            }
        }
    }

    private void updateResults()
    {
        /* record the following
            - number of segments (fixed for each test run)
            - number of links used
            - number of adjacent links preserved
        */

        int segmentsLinked = 0;
        int exactMatchCount = 0;

        Segment first = mSegments.get(0);
        Segment last = mSegments.get(mSegments.size()-1);

        Segment currentSegment = first;
        boolean nextLinkOnStart = false;

        List<Integer> linkedIndices = Lists.newLinkedList();
        linkedIndices.add(currentSegment.Id);

        int iterations = 0;

        while(true)
        {
            Segment nextSegment = currentSegment.getLink(nextLinkOnStart);

            if(nextSegment == null)
            {
                LOGGER.error("incomplete chain");
                mValidRun = false;
                break;
            }

            linkedIndices.add(nextSegment.Id);

            if(currentSegment.Id == (nextSegment.Id - 1)
            && currentSegment.getLink(false) == nextSegment && nextSegment.getLink(true) == currentSegment)
            {
                ++exactMatchCount;
            }
            else if(nextSegment.Id == (currentSegment.Id - 1)
            && nextSegment.getLink(false) == currentSegment && currentSegment.getLink(true) == nextSegment)
            {
                ++exactMatchCount;
            }

            if(nextSegment == last)
                break;

            ++segmentsLinked;

            boolean nextLinkedOnStart = nextSegment.getLink(true) == currentSegment;
            nextLinkOnStart = !nextLinkedOnStart;
            currentSegment = nextSegment;

            ++iterations;

            if(iterations > mSegments.size())
            {
                LOGGER.error("iterations({}) exceeds segCount({})", iterations, mSegments.size());
                mValidRun = false;
                return;
            }
        }

        // test now many adjacent segments are present, even if not directly next to each other
        // these are equivalent to deletion bridges (with zero loss)
        int adjacentPairs = 0;

        for(int i = 0; i < mSegments.size() - 1; ++i)
        {
            if(linkedIndices.contains(i) && linkedIndices.contains(i+1))
                ++adjacentPairs;
        }

        mLatestResult = ImmutableShatteringResult.builder()
                .testCount(mConfig.Iterations)
                .runIndex(mRunIndex)
                .segments(mSegmentCount)
                .linkedSegments(segmentsLinked)
                .exactMatches(exactMatchCount)
                .adjacentSegments(adjacentPairs)
                .build();

        LOGGER.debug("run({}) results: links({}) exact({}) adj({})", mRunIndex, segmentsLinked, exactMatchCount, adjacentPairs);
    }

    private int[] getNextIndexPair(int bounds)
    {
        int[] indices = new int[2];

        if(mSpecifiedOrder.size() >= 2)
        {
            indices[0] = min(mSpecifiedOrder.get(0), mRemainingLinkCount);
            mSpecifiedOrder.remove(0);
            indices[1] = min(mSpecifiedOrder.get(0), mRemainingLinkCount);
            mSpecifiedOrder.remove(0);
        }
        else
        {
            if(bounds <= 1)
            {
                mValidRun = false;
                LOGGER.error("request for random ints below bound of 1");
                return indices;
            }

            if(bounds <= 2)
            {
                indices[0] = 0;
                indices[1] = 1;
            }
            else
            {
                indices[0] = mRandom.nextInt(bounds);
                indices[1] = mRandom.nextInt(bounds);

                while(indices[1] == indices[0])
                {
                    indices[1] = mRandom.nextInt(bounds);
                }
            }
        }

        return indices;
    }

    public final ShatteringResult getLatestResults()
    {
        return mLatestResult;
    }

    private void initialiseWriter()
    {
        if(mOutputDir.isEmpty())
            return;

        try
        {
            String outputFileName = mOutputDir + "SIM_RESULTS.csv";

            mResultsWriter = createBufferedWriter(outputFileName, false);

            // definitional fields
            mResultsWriter.write("TestIterations,TestRun,SegCount,SegsLinked,ExactMatches,AdjacentSegs");

            mResultsWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile: {}", e.toString());
        }
    }

    private void writeResults()
    {
        if(mResultsWriter == null)
            return;

        try
        {
            mResultsWriter.write(String.format("%d,%d,%d,%d,%d,%d",
                    mLatestResult.testCount(), mLatestResult.runIndex(), mLatestResult.segments(),
                    mLatestResult.linkedSegments(), mLatestResult.exactMatches(), mLatestResult.adjacentSegments()));

            mResultsWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mResultsWriter);
    }

    public void setSpecifiedOrder(final List<Integer> order)
    {
        mSpecifiedOrder.clear();
        mSpecifiedOrder.addAll(order);
    }

}
