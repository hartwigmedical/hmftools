package com.hartwig.hmftools.svanalysis.simulation;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvSimShattering
{
    private static final Logger LOGGER = LogManager.getLogger(SvSimShattering.class);

    // config
    private int mIterations;
    private int mSegmentCount;
    private List<SvSimSegment> mSegments;
    private int mRemainingLinkCount;
    private String mOutputDir;

    // state and results
    private Random mRandom;
    private int mRunIndex;
    private boolean mValidRun;
    private int[] mLatestResults;


    private BufferedWriter mResultsWriter;

    // testing
    private List<Integer> mSpecifiedOrder;

    public SvSimShattering()
    {
        mSegments = Lists.newArrayList();
        mRemainingLinkCount = 0;
        mValidRun = true;
        mOutputDir = "";
        mResultsWriter = null;
        mLatestResults = new int[SS_RESULTS_ADJACENT_SEGS+1];

        mRandom = new Random(123456);
        mSpecifiedOrder = Lists.newArrayList();
    }

    public void setOutputDir(final String dir)
    {
        mOutputDir = dir;
        initialiseWriter();
    }

    public boolean validRun() { return mValidRun; }

    public void initialise(int segmentCount, int iterations)
    {
        mRunIndex = 0;
        mIterations = iterations;
        mSegments.clear();
        mSegmentCount = segmentCount;
        mValidRun = true;

        mSegments.add(new SvSimSegment(0, false, true));

        for(int i = 0; i < mSegmentCount; ++i)
        {
            mSegments.add(new SvSimSegment(i+1, true, true));
        }

        mSegments.add(new SvSimSegment(mSegments.size(), true, false));
    }

    public void run()
    {
        for(int i = 0; i < mIterations; ++i)
        {
            runIteration();
            updateResults();
            logResults();
            ++mRunIndex;

            if(mIterations > 10000 && (i % 10000) == 1)
            {
                LOGGER.info("run index {}", i);
            }

            if(!mValidRun)
                break;
        }
    }

    private void clearRunState()
    {
        mRemainingLinkCount = (mSegmentCount + 1 ) * 2;
        mSegments.stream().forEach(x -> x.clearLinks());
    }

    // private static int SPEC_RUN_INDEX = -1;
    private static int SPEC_RUN_INDEX = 23;

    private void runIteration()
    {
        clearRunState();

        boolean foundLink = true;

        if(mRunIndex == SPEC_RUN_INDEX)
        {
            LOGGER.debug("spec index({})", mRunIndex);
        }

        while(mRemainingLinkCount > 0 && foundLink)
        {
            int[] nextIndices = getNextIndexPair(mRemainingLinkCount);
            int nextIndex1 = nextIndices[0];
            int nextIndex2 = nextIndices[1];
            int index = 0;

            SvSimSegment nextSegment1 = null;
            boolean seg1LinkOnStart = false;
            SvSimSegment nextSegment2 = null;
            boolean seg2LinkOnStart = false;

            foundLink = false;

            int segIndex = 0;
            for(; segIndex < mSegments.size(); ++segIndex)
            {
                SvSimSegment segment = mSegments.get(segIndex);

                if(!segment.fullyLinked())
                {
                    for (int be = SVI_START; be <= SVI_END; ++be)
                    {
                        boolean isStart = isStart(be);
                        if (!segment.isLinkOpen(isStart))
                            continue;

                        if (nextSegment1 == null && index >= nextIndex1) // && nextSegment2 != segment
                        {
                            nextSegment1 = segment;
                            seg1LinkOnStart = isStart;
                        }
                        else if (nextSegment2 == null && index >= nextIndex2) //&& nextSegment1 != segment
                        {
                            nextSegment2 = segment;
                            seg2LinkOnStart = isStart;
                        }

                        ++index;
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

                if(segIndex == mSegments.size() - 1)
                {
                    // 2 free segments may not be found if the 2 random indices are in the last open segment with both links available
                    segIndex = 0;
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
        SvSimSegment first = mSegments.get(0);

        if(first.isLinkOpen(false))
            return true;

        SvSimSegment last = mSegments.get(mSegments.size()-1);

        if(last.isLinkOpen(true))
            return true;

        // see a path can be traced from first to last
        SvSimSegment currentSegment = first;
        boolean nextLinkOnStart = false;

        int iterations = 0;

        while(true)
        {
            // check other side
            SvSimSegment nextSegment = currentSegment.getLink(nextLinkOnStart);

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

    public static int SS_RESULTS_TEST_COUNT = 0;
    public static int SS_RESULTS_TEST_RUN = 1;
    public static int SS_RESULTS_SEGMENTS = 2;
    public static int SS_RESULTS_SEGS_LINKED = 3;
    public static int SS_RESULTS_EXACT_MATCHES = 4;
    public static int SS_RESULTS_ADJACENT_SEGS = 5;

    private void updateResults()
    {
        /* record the following
            - number of segments (fixed for each test run)
            - number of links used
            - number of adjacent links preserved
        */
        mLatestResults = new int[SS_RESULTS_ADJACENT_SEGS+1];

        int segmentsLinked = 0;
        int exactMatchCount = 0;

        SvSimSegment first = mSegments.get(0);
        SvSimSegment last = mSegments.get(mSegments.size()-1);

        SvSimSegment currentSegment = first;
        boolean nextLinkOnStart = false;

        List<Integer> linkedIndices = Lists.newLinkedList();
        linkedIndices.add(currentSegment.Id);

        int iterations = 0;

        while(true)
        {
            SvSimSegment nextSegment = currentSegment.getLink(nextLinkOnStart);

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

        // test now many consecutive segments are present
        int adjacentPairs = 0;

        for(int i = 0; i < mSegments.size() - 1; ++i)
        {
            if(linkedIndices.contains(i) && linkedIndices.contains(i+1))
                ++adjacentPairs;
        }

        mLatestResults[SS_RESULTS_TEST_COUNT] = mIterations;
        mLatestResults[SS_RESULTS_TEST_RUN] = mRunIndex;
        mLatestResults[SS_RESULTS_SEGMENTS] = mSegmentCount;
        mLatestResults[SS_RESULTS_SEGS_LINKED] = segmentsLinked;
        mLatestResults[SS_RESULTS_EXACT_MATCHES] = exactMatchCount;
        mLatestResults[SS_RESULTS_ADJACENT_SEGS] = adjacentPairs;

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

    public final int[] getLatestResults()
    {
        return mLatestResults;
    }

    public void logResults()
    {
        writeResults();
    }

    private void initialiseWriter()
    {
        try
        {
            String outputFileName = mOutputDir + "SIM_RESULTS.csv";

            mResultsWriter = createBufferedWriter(outputFileName, false);

            // definitional fields
            mResultsWriter.write("TestCount,TestRun,SegCount,SegsLinked,ExactMatches,AdjacentSegs");

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
                    mLatestResults[SS_RESULTS_TEST_COUNT], mLatestResults[SS_RESULTS_TEST_RUN],
                    mLatestResults[SS_RESULTS_SEGMENTS], mLatestResults[SS_RESULTS_SEGS_LINKED],
                    mLatestResults[SS_RESULTS_EXACT_MATCHES], mLatestResults[SS_RESULTS_ADJACENT_SEGS]));

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
