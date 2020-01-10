package com.hartwig.hmftools.linx.simulation;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

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
    private int mMaxLinkIndex;
    private final List<Integer> mRemainingLinks; // list of integers representing the unlinked segment start & end points
    private final Random mRandom;
    private int mRunIndex;
    private boolean mValidRun;
    private ShatteringResult mLatestResult;

    private Map<ShatteringResult,Integer> mGroupedResults;

    private BufferedWriter mResultsWriter;

    // for unit testing only - the order in which to add links
    private final List<Integer> mSpecifiedLinkOrder;

    private static final Logger LOGGER = LogManager.getLogger(ShatteringSim.class);

    public ShatteringSim(final ShatteringConfig config, final String outputDir)
    {
        mConfig = config;
        mOutputDir = outputDir;
        initialiseWriter();

        mSegments = Lists.newArrayList();
        mRemainingLinks = Lists.newArrayList();
        mSegmentCount = 0;
        mValidRun = true;
        mLatestResult = null;

        mRandom = new Random(123456);
        mSpecifiedLinkOrder = Lists.newArrayList();
        mGroupedResults = Maps.newHashMap();
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
        mGroupedResults.clear();
        initialiseState();

        for(int i = 0; i < mConfig.Iterations; ++i)
        {
            runIteration();
            updateResults();

            if(mConfig.GroupResults)
                registerResult(mLatestResult);
            else
                writeResults(mLatestResult, 1);

            ++mRunIndex;

            if(mConfig.Iterations > 10000 && i > 0 && (i % 10000) == 0)
            {
                LOGGER.info("run index {}", i);
            }

            if(!mValidRun)
                break;
        }

        if(mConfig.GroupResults)
        {
            for (Map.Entry<ShatteringResult, Integer> entry : mGroupedResults.entrySet())
            {
                writeResults(entry.getKey(), entry.getValue());
            }
        }
    }

    private void registerResult(final ShatteringResult result)
    {
        for(Map.Entry<ShatteringResult,Integer> entry : mGroupedResults.entrySet())
        {
            final ShatteringResult existingResult = entry.getKey();

            if(existingResult.equals(result))
            {
                entry.setValue(entry.getValue() + 1);
                return;
            }
        }

        mGroupedResults.put(result, 1);
    }

    private void initialiseState()
    {
        mRunIndex = 0;
        mSegments.clear();
        mValidRun = true;

        // create a set of N segments with 2 unconnected ends, and 2 bounding segments with a single exposed end
        // so in total there are N+2 segments
        mSegments.add(new Segment(0, false, true));

        for(int i = 0; i < mSegmentCount; ++i)
        {
            mSegments.add(new Segment(i+1, true, true));
        }

        mSegments.add(new Segment(mSegments.size(), true, false));
    }

    private void clearRunState()
    {
        mMaxLinkIndex = calcLinkCount(mSegments.size());
        mRemainingLinks.clear();

        for(int linkIndex = 0; linkIndex <= mMaxLinkIndex; ++linkIndex)
        {
            mRemainingLinks.add(linkIndex);
        }

        mSegments.stream().forEach(x -> x.clearLinks());
    }

    private void runIteration()
    {
        clearRunState();

        while(mRemainingLinks.size() >= 2)
        {
            // randomly find the next 2 ends to connect
            int[] nextIndices = getNextSegmentLinks();
            int nextIndex1 = nextIndices[0];
            int nextIndex2 = nextIndices[1];

            Segment nextSegment1 = getSegmentByLinkIndex(nextIndex1);
            boolean seg1LinkOnStart = isSegmentStartByLinkIndex(nextIndex1);
            Segment nextSegment2 = getSegmentByLinkIndex(nextIndex2);
            boolean seg2LinkOnStart = isSegmentStartByLinkIndex(nextIndex2);

            if(nextSegment1 == null || nextSegment2 == null)
            {
                LOGGER.error("invalid segment lookup with nextIndex({} & {})", nextIndex1, nextIndex2);
                mValidRun = false;
                break;
            }

            nextSegment1.setLink(nextSegment2, seg1LinkOnStart);
            nextSegment2.setLink(nextSegment1, seg2LinkOnStart);

            LOGGER.debug("linked indices({}-{} to {}-{}) remaining({})",
                    nextSegment1.Id, seg1LinkOnStart ? "start" : "end", nextSegment2.Id, seg2LinkOnStart ? "start" : "end",
                    mRemainingLinks.size());

            if(!moreLinksPossible())
            {
                LOGGER.debug("exiting with no more possible links, remainingLinks({})", nextIndex1, nextIndex2, mRemainingLinks.size());
                break;
            }
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

        // see if a path can be traced from first to last
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
        int inferredLinks = 0;
        boolean firstNonExactSeen = false;

        Segment first = mSegments.get(0);
        Segment last = mSegments.get(mSegments.size()-1);

        Segment currentSegment = first;
        boolean nextLinkOnStart = false;

        final List<Integer> linkedIndices = Lists.newLinkedList();
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

            boolean exactRepair = false;

            if(currentSegment.Id == (nextSegment.Id - 1)
            && currentSegment.getLink(false) == nextSegment && nextSegment.getLink(true) == currentSegment)
            {
                exactRepair = true;
            }
            else if(nextSegment.Id == (currentSegment.Id - 1)
            && nextSegment.getLink(false) == currentSegment && currentSegment.getLink(true) == nextSegment)
            {
                exactRepair = true;
            }

            if(exactRepair)
            {
                ++exactMatchCount;
            }
            else// if(nextSegment != last)
            {
                if(!firstNonExactSeen)
                    firstNonExactSeen = true;
                else
                    ++inferredLinks;
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

        final List<Segment> lostSegments = mSegments.stream()
                .filter(x -> !x.EndSegment)
                .filter(x -> !linkedIndices.contains(x.Id))
                .collect(Collectors.toList());

        // determine contiguous lost sections (ie deletion bridges or DELs)
        int inferredLost = 0;

        if(!lostSegments.isEmpty())
        {
            if(lostSegments.size() == mSegmentCount)
            {
                // all segment lost can be just considered a single DEL event
                inferredLost = 1;
            }
            else
            {
                int index = 0;
                while (index < lostSegments.size())
                {
                    ++inferredLost;

                    int nextIndex = index + 1;

                    while (nextIndex < lostSegments.size())
                    {
                        if (lostSegments.get(nextIndex-1).Id + 1 == lostSegments.get(nextIndex).Id)
                            ++nextIndex;
                        else
                            break;
                    }

                    index = nextIndex;
                }
            }
        }

        if(segmentsLinked == 1 && inferredLinks == 1 && inferredLost != 2)
        {

        }

        mLatestResult = ImmutableShatteringResult.builder()
                .runIndex(mRunIndex)
                .segments(mSegmentCount)
                .linkedSegments(segmentsLinked)
                .exactRepairs(exactMatchCount)
                .adjacentSegments(adjacentPairs)
                .inferredLinks(inferredLinks)
                .inferredLost(inferredLost)
                .build();

        LOGGER.debug("run({}) results: links(kept={} lost={} exact={} adj={}) inferred(links={} lost={})",
                mRunIndex, segmentsLinked, lostSegments.size(), exactMatchCount, adjacentPairs, inferredLinks, inferredLost);
    }

    private int[] getNextSegmentLinks()
    {
        int[] indices = new int[2];

        if(mSpecifiedLinkOrder.size() >= 2)
        {
            indices[0] = mSpecifiedLinkOrder.get(0);
            mSpecifiedLinkOrder.remove(0);
            indices[1] = mSpecifiedLinkOrder.get(0);
            mSpecifiedLinkOrder.remove(0);
        }
        else
        {
            if (mRemainingLinks.size() == 2)
            {
                indices[0] = mRemainingLinks.get(0);
                indices[1] = mRemainingLinks.get(1);
                mRemainingLinks.clear();
            }
            else
            {
                int index = mRandom.nextInt(mRemainingLinks.size());
                Integer index0 = mRemainingLinks.get(index);
                mRemainingLinks.remove(index0);
                indices[0] = index0;

                index = mRandom.nextInt(mRemainingLinks.size());
                Integer index1 = mRemainingLinks.get(index);
                mRemainingLinks.remove(index1);
                indices[1] = index1;
            }
        }

        return indices;
    }

    public static int calcLinkCount(int segmentCount)
    {
        // for segments 0, 1, 2 and 3, there will be index 0, 1-2, 3-4 and 5
        return (segmentCount - 1) * 2 - 1;
    }

    private Segment getSegmentByLinkIndex(int linkIndex)
    {
        if(linkIndex < 0 || linkIndex > mMaxLinkIndex)
        {
            LOGGER.error("invalid linkIndex({}) vs segmentCount({})", linkIndex, mSegments.size());
            return null;
        }

        if(linkIndex == 0)
            return mSegments.get(0);

        int segIndex = linkIndex / 2;

        if((linkIndex % 2) == 1)
            ++segIndex;

        return mSegments.get(segIndex);
    }

    private boolean isSegmentStartByLinkIndex(int linkIndex)
    {
        if(linkIndex == 0)
            return false;
        else if(linkIndex == mMaxLinkIndex)
            return true;
        else
            return (linkIndex % 2) == 1;
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
            String outputFileName = mOutputDir + "LNX_SIM_RESULTS.csv";

            mResultsWriter = createBufferedWriter(outputFileName, false);

            // definitional fields
            mResultsWriter.write("TestIterations,SegCount,SegsLinked,ExactRepairs,AdjacentSegs,InfLinks,InfLost");

            if(mConfig.GroupResults)
                mResultsWriter.write(",RepeatCount");
            else
                mResultsWriter.write(",TestRun");

            mResultsWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile: {}", e.toString());
        }
    }

    private void writeResults(final ShatteringResult result, int repeatCount)
    {
        if(mResultsWriter == null)
            return;

        try
        {
            mResultsWriter.write(String.format("%d,%d,%d,%d,%d,%d,%d",
                    mConfig.Iterations, result.segments(), result.linkedSegments(),
                    result.exactRepairs(), result.adjacentSegments(), result.inferredLinks(), result.inferredLost()));

            if(mConfig.GroupResults)
            {
                mResultsWriter.write(String.format(",%d", repeatCount));
            }
            else
            {
                mResultsWriter.write(String.format(",%d", result.runIndex()));
            }

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
        mSpecifiedLinkOrder.clear();
        mSpecifiedLinkOrder.addAll(order);
    }

}
