package com.hartwig.hmftools.svtools.simulation;

import static java.lang.Math.floor;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;

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
    private final List<Integer> mRemainingLinks;
    private final List<int[]> mRemainingLinkPairs;
    private final Random mRandom;
    private int mRunIndex;
    private boolean mValidRun;
    private ShatteringResult mLastResult;
    private String mLastLinkStr;
    private final List<int[]> mIndexSelector;

    private final Map<ShatteringResult,Integer> mGroupedResults;

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
        mRemainingLinkPairs = Lists.newArrayList();
        mSegmentCount = 0;
        mValidRun = true;
        mLastResult = null;
        mLastLinkStr = "";

        mRandom = new Random();
        mGroupedResults = Maps.newHashMap();

        mIndexSelector = Lists.newArrayList();
        mSpecifiedLinkOrder = Lists.newArrayList();
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

        if(!mConfig.ExhuastiveSearch)
        {
            for (int i = 0; i < mConfig.Iterations; ++i)
            {
                runIteration();

                if (!mValidRun)
                    break;

                mLastResult = generateResults(mSegments, mLastLinkStr);

                if (mConfig.GroupResults)
                    registerResult(mLastResult);
                else
                    writeResults(mLastResult, 1);

                ++mRunIndex;

                if (mConfig.Iterations > 10000 && i > 0 && (i % 10000) == 0)
                {
                    LOGGER.info("run index {}", i);
                }
            }
        }
        else
        {
            runIterationRecursively();
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

        for(int i = 0; i <= mSegmentCount; ++i)
        {
            mIndexSelector.add(new int[]{0,0});
        }

        // create a set of N segments with 2 unconnected ends, and 2 bounding segments with a single exposed end
        // so in total there are N+2 segments
        mSegments.add(new Segment(0, false, true));

        for(int i = 0; i < mSegmentCount; ++i)
        {
            mSegments.add(new Segment(i+1, true, true));
        }

        mSegments.add(new Segment(mSegments.size(), true, false));
        mMaxLinkIndex = calcLinkCount(mSegments.size());
    }

    private void clearRunState()
    {
        mRemainingLinks.clear();

        for (int i = 0; i <= mMaxLinkIndex; ++i)
        {
            mRemainingLinks.add(i);
        }

        mSegments.stream().forEach(x -> x.clearLinks());
    }

    private void runIteration()
    {
        clearRunState();

        String linksStr = "";
        int roundIndex = 0;

        while(!mRemainingLinks.isEmpty())
        {
            // randomly find the next 2 ends to connect
            int[] nextIndices = getNextSegmentLinks(roundIndex);
            int nextIndex1 = nextIndices[0];
            int nextIndex2 = nextIndices[1];

            Segment nextSegment1 = getSegmentByLinkIndex(nextIndex1, mSegments);
            boolean seg1LinkOnStart = isSegmentStartByLinkIndex(nextIndex1);
            Segment nextSegment2 = getSegmentByLinkIndex(nextIndex2, mSegments);
            boolean seg2LinkOnStart = isSegmentStartByLinkIndex(nextIndex2);

            if(nextSegment1 == null || nextSegment2 == null)
            {
                LOGGER.error("invalid segment lookup with nextIndex({} & {})", nextIndex1, nextIndex2);
                mValidRun = false;
                break;
            }

            nextSegment1.setLink(nextSegment2, seg1LinkOnStart);
            nextSegment2.setLink(nextSegment1, seg2LinkOnStart);

            if(LOGGER.isDebugEnabled() || mSegmentCount <= 5)
            {
                String link = String.format("%d:%s-%d:%s",
                        nextSegment1.Id, seg1LinkOnStart ? "s" : "e", nextSegment2.Id, seg2LinkOnStart ? "s" : "e");

                LOGGER.debug("{}: linked({}) remaining links({})", roundIndex, link, mRemainingLinks.size());

                linksStr = appendStr(linksStr, link, ';');
            }

            if(!moreLinksPossible(mSegments))
            {
                LOGGER.debug("exiting with no more possible links, remaining links({})", mRemainingLinks.size());
                break;
            }

            ++roundIndex;
        }

        mLastLinkStr = linksStr;
    }

    private void runIterationRecursively()
    {
        for (int i = 0; i <= mMaxLinkIndex - 1; ++i)
        {
            for (int j = i + 1; j <= mMaxLinkIndex; ++j)
            {
                mRemainingLinkPairs.add(new int[] { i, j });
            }
        }

        mSegments.stream().forEach(x -> x.clearLinks());

        String linkStr = "";
        findNextLink(mSegments, mRemainingLinkPairs, linkStr);
    }

    private void findNextLink(final List<Segment> segments, final List<int[]> remainingLinks, final String linkStr)
    {
        // take the next possible link off and call again recursively
        if(remainingLinks.isEmpty() || !mValidRun)
            return;

        for(int linkIndex = 0; linkIndex < remainingLinks.size(); ++linkIndex)
        {
            List<int[]> newRemainingLinks = Lists.newArrayList();
            newRemainingLinks.addAll(remainingLinks);

            int[] nextIndices = newRemainingLinks.get(linkIndex);
            int nextIndex1 = nextIndices[0];
            int nextIndex2 = nextIndices[1];

            final List<Segment> newSegments = segments.stream().map(x -> Segment.from(x)).collect(Collectors.toList());
            String newLinkStr = linkStr;

            // recreate existing links
            for(int i = 0; i < segments.size(); ++i)
            {
                final Segment refSegment = segments.get(i);
                Segment copySegment = newSegments.get(i);

                if(refSegment.getLink(true) != null)
                {
                    copySegment.setLink(newSegments.get(refSegment.getLink(true).Id), true);
                }

                if(refSegment.getLink(false) != null)
                {
                    copySegment.setLink(newSegments.get(refSegment.getLink(false).Id), false);
                }
            }

            Segment nextSegment1 = getSegmentByLinkIndex(nextIndex1, newSegments);
            boolean seg1LinkOnStart = isSegmentStartByLinkIndex(nextIndex1);
            Segment nextSegment2 = getSegmentByLinkIndex(nextIndex2, newSegments);
            boolean seg2LinkOnStart = isSegmentStartByLinkIndex(nextIndex2);

            if(nextSegment1 == null || nextSegment2 == null)
            {
                LOGGER.error("invalid segment lookup with nextIndex({} & {})", nextIndex1, nextIndex2);
                mValidRun = false;
                return;
            }

            nextSegment1.setLink(nextSegment2, seg1LinkOnStart);
            nextSegment2.setLink(nextSegment1, seg2LinkOnStart);

            if(LOGGER.isDebugEnabled() || mSegmentCount <= 5)
            {
                String link = String.format("%d:%s-%d:%s",
                        nextSegment1.Id, seg1LinkOnStart ? "s" : "e", nextSegment2.Id, seg2LinkOnStart ? "s" : "e");

                newLinkStr = appendStr(linkStr, link, ';');

                LOGGER.debug("linked({}) remaining linkPairs({})", link, remainingLinks.size());
            }

            if(!moreLinksPossible(newSegments))
            {
                LOGGER.debug("no more possible links, remaining linkPairs({})", remainingLinks.size());

                // register the result
                ShatteringResult result = generateResults(newSegments, newLinkStr);

                if (mConfig.GroupResults)
                    registerResult(result);
                else
                    writeResults(result, 1);

                continue;
            }

            // prepare state for the next call - first remove conflicting links
            purgeConflictingLinks(newRemainingLinks, nextIndices);

            if(newRemainingLinks.isEmpty())
            {
                LOGGER.error("no more possible links");
                continue;
            }

            findNextLink(newSegments, newRemainingLinks, newLinkStr);
        }
    }

    private boolean moreLinksPossible(final List<Segment> segments)
    {
        Segment first = segments.get(0);

        if(first.isLinkOpen(false))
            return true;

        Segment last = segments.get(segments.size()-1);

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

            if(iterations > segments.size())
            {
                LOGGER.error("iterations({}) exceeds segCount({})", iterations, segments.size());
                mValidRun = false;
                return false;
            }
        }
    }

    private ShatteringResult generateResults(final List<Segment> segments, final String linkStr)
    {
        /* record the following
            - number of segments (fixed for each test run)
            - number of links used
            - number of adjacent links preserved
            - inferred TIs and DELs/DBs
        */

        int segmentsLinked = 0;
        int exactMatchCount = 0;
        int inferredLinks = 0;
        boolean firstNonExactSeen = false;

        Segment first = segments.get(0);
        Segment last = segments.get(segments.size()-1);

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

            if(iterations > segments.size())
            {
                LOGGER.error("iterations({}) exceeds segCount({})", iterations, segments.size());
                mValidRun = false;
                return null;
            }
        }

        // test now many adjacent segments are present, even if not directly next to each other
        // these are equivalent to deletion bridges (with zero loss)
        int adjacentPairs = 0;

        for(int i = 0; i < segments.size() - 1; ++i)
        {
            if(linkedIndices.contains(i) && linkedIndices.contains(i+1))
                ++adjacentPairs;
        }

        final List<Segment> lostSegments = segments.stream()
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

        LOGGER.debug("run({}) results: links(kept={} lost={} exact={} adj={}) inferred(links={} lost={}) linkStr({})",
                mRunIndex, segmentsLinked, lostSegments.size(), exactMatchCount, adjacentPairs,
                inferredLinks, inferredLost, linkStr);

        return ImmutableShatteringResult.builder()
                .runIndex(mRunIndex)
                .segments(mSegmentCount)
                .linkedSegments(segmentsLinked)
                .exactRepairs(exactMatchCount)
                .adjacentSegments(adjacentPairs)
                .inferredLinks(inferredLinks)
                .inferredLost(inferredLost)
                .linkStr(linkStr)
                .build();
    }

    private int[] getNextSegmentLinks(int roundIndex)
    {
        final int[] nextIndices = new int[2];

        if(mSpecifiedLinkOrder.size() >= 2)
        {
            nextIndices[0] = mSpecifiedLinkOrder.get(0);
            nextIndices[1] = mSpecifiedLinkOrder.get(1);
            mSpecifiedLinkOrder.remove(0);
            mSpecifiedLinkOrder.remove(0);
        }
        else
        {
            if(mRemainingLinks.size() == 2)
            {
                nextIndices[0] = mRemainingLinks.get(0);
                nextIndices[1] = mRemainingLinks.get(1);
                mRemainingLinks.clear();
            }
            else
            {
                int randIndex = mRandom.nextInt(mRemainingLinks.size());
                Integer nextIndex = mRemainingLinks.get(randIndex);
                nextIndices[0] = nextIndex;
                mRemainingLinks.remove(nextIndex);

                randIndex = mRandom.nextInt(mRemainingLinks.size());
                nextIndex = mRemainingLinks.get(randIndex);
                nextIndices[1] = nextIndex;
                mRemainingLinks.remove(nextIndex);
            }
        }

        return nextIndices;
    }

    public void purgeConflictingLinks(final List<int[]> existingLinks, final int[] newLinks)
    {
        int lpIndex = 0;
        while(lpIndex < existingLinks.size())
        {
            int[] pair = existingLinks.get(lpIndex);
            if(pair[0] == newLinks[0] || pair[1] == newLinks[1] || pair[0] == newLinks[1] || pair[1] == newLinks[0])
                existingLinks.remove(lpIndex);
            else
                ++lpIndex;
        }
    }

    private int getNextRandom()
    {
        double rand = mRandom.nextDouble();
        return (int)floor(mRemainingLinks.size() * rand);
    }

    public static int calcLinkCount(int segmentCount)
    {
        // for segments 0, 1, 2 and 3, there will be index 0, 1-2, 3-4 and 5
        return (segmentCount - 1) * 2 - 1;
    }

    private Segment getSegmentByLinkIndex(int linkIndex, final List<Segment> segments)
    {
        if(linkIndex < 0 || linkIndex > mMaxLinkIndex)
        {
            LOGGER.error("invalid linkIndex({}) vs segmentCount({})", linkIndex, segments.size());
            return null;
        }

        if(linkIndex == 0)
            return segments.get(0);

        int segIndex = linkIndex / 2;

        if((linkIndex % 2) == 1)
            ++segIndex;

        return segments.get(segIndex);
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
        return mLastResult;
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
                mResultsWriter.write(",TestRun,LinkStr");

            mResultsWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile: {}", e.toString());
        }
    }

    private void writeResults(final ShatteringResult result, int repeatCount)
    {
        if(mResultsWriter == null || !mValidRun)
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
                mResultsWriter.write(String.format(",%d,%s", result.runIndex(), result.linkStr()));
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
