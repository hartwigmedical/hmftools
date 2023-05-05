package com.hartwig.hmftools.markdups;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.markdups.common.Fragment;

import htsjdk.samtools.SAMRecord;

public class ReadPositionsCache
{
    // a ring buffer to store reads at each read starting position
    private final String mChromosome;
    private final FragmentGroup[] mForwardPositions;
    private final Map<Integer,FragmentGroup> mReversePositions;
    private final Map<String, Fragment> mFragments;
    private final Consumer<List<Fragment>> mReadGroupHandler;
    private int mMinPosition;
    private int mMinPositionIndex;
    private final int mCapacity;
    private final boolean mUseMateCigar;

    private int mLastFragmentLogCount;
    private int mLastLogReadCount;

    private long mFragmemtCacheCount;
    private long mFragmemtUnmatchedCount;
    private long mFragmemtNoCacheCount;

    private class FragmentGroup
    {
        // fragments with a matching start position
        public final List<Fragment> Fragments;

        public FragmentGroup(final Fragment fragment)
        {
            Fragments = Lists.newArrayList(fragment);
        }
    }

    public ReadPositionsCache(
            final String chromosome, int capacity, boolean useMateCigar, final Consumer<List<Fragment>> evictionHandler)
    {
        mChromosome = chromosome;
        mReadGroupHandler = evictionHandler;
        mCapacity = capacity;
        mForwardPositions = new FragmentGroup[mCapacity];
        mReversePositions = Maps.newHashMap();
        mFragments = Maps.newHashMap();
        mMinPosition = 0;
        mMinPositionIndex = 0;
        mUseMateCigar = useMateCigar;

        mLastFragmentLogCount = 0;
        mLastLogReadCount = 0;
        mFragmemtCacheCount = 0;
        mFragmemtUnmatchedCount = 0;
        mFragmemtNoCacheCount = 0;
    }

    public boolean processRead(final SAMRecord read)
    {
        // supplementaries just check for a fragment match otherwise no further processing
        if(read.getSupplementaryAlignmentFlag())
        {
            Fragment fragment = mFragments.get(read.getReadName());

            if(fragment != null) // add to fragment if in a current group
            {
                fragment.addRead(read);
                return true;
            }

            return false;
        }

        // check for an existing mate if on the same chromosome
        // store in a group of fragments with a matching first fragment coordinate
        // if the mate has a lower position or is on a lower chromosome, don't add it to a position group
        boolean mateUnmapped = read.getMateUnmappedFlag();
        boolean readUnmapped = read.getReadUnmappedFlag();

        boolean sameChromosome = !mateUnmapped && read.getMateReferenceName().equals(mChromosome);

        // skip if mate is on a lower chromosome
        if(!sameChromosome && read.getReferenceIndex() > read.getMateReferenceIndex())
            return false;

        Fragment fragment = mFragments.get(read.getReadName());

        if(fragment != null) // add to fragment if in a current group
        {
            fragment.addRead(read);
            return true;
        }

        if(readUnmapped)
            return false;

        if(sameChromosome && read.getAlignmentStart() > read.getMateAlignmentStart()) // mate already processed and evicted
            return false;

        storeInitialRead(read);
        return true;
    }

    private void storeInitialRead(final SAMRecord read)
    {
        ++mLastLogReadCount;

        Fragment fragment = new Fragment(read);
        fragment.intialiseCoordinates(mUseMateCigar);

        int fragmentPosition = fragment.initialPosition();

        if(fragmentPosition > 0)
        {
            if(mMinPosition == 0)
                resetMinPosition(fragmentPosition);
            else
                checkFlush(fragmentPosition);
        }

        // only store the read if its mate is local and expected to be added
        if(fragment.hasLocalMate())
        {
            mFragments.put(read.getReadName(), fragment);
            ++mFragmemtCacheCount;
        }
        else
        {
            ++mFragmemtNoCacheCount;
        }

        if(fragmentPosition > 0)
        {
            int index = calcIndex(fragmentPosition);

            if(index < 0 || index >= mCapacity)
            {
                MD_LOGGER.error("fragment({}) outside forward strand array bounds", fragment);
                return;
            }

            FragmentGroup element = mForwardPositions[index];
            if(element == null)
            {
                element = new FragmentGroup(fragment);
                mForwardPositions[index] = element;
            }
            else
            {
                element.Fragments.add(fragment);
            }
        }
        else
        {
            // store in reverse strand map
            FragmentGroup element = mReversePositions.get(fragmentPosition);
            if(element == null)
            {
                element = new FragmentGroup(fragment);
                mReversePositions.put(fragmentPosition, element);
            }
            else
            {
                element.Fragments.add(fragment);
            }
        }
    }

    private int calcIndex(int position)
    {
        // capacity = 10, min position = 1, min index = 0, position of 10 is index 9
        // capacity = 10, min position = 2, min index = 1, position of 10 is index 9, position of 11 is 0
        int distanceFromMinPosition = position - mMinPosition;

        if(mMinPositionIndex + distanceFromMinPosition < mCapacity)
            return mMinPositionIndex + distanceFromMinPosition;

        // index is from start of ring buffer
        return distanceFromMinPosition + mMinPositionIndex - mCapacity;
    }

    private void checkFlush(int position)
    {
        if(mMinPosition == 0)
        {
            resetMinPosition(position);
            return;
        }

        int distanceFromMinPosition = position - mMinPosition;

        if(distanceFromMinPosition < mCapacity)
            return;

        int flushCount = position - mMinPosition - mCapacity + 1;

        int flushedElements = 0;

        // only iterate at most once through the array
        for(int i = 0; i < min(flushCount, mCapacity); i++)
        {
            FragmentGroup element = mForwardPositions[mMinPositionIndex];

            // clear and process each element and depth
            if(element != null)
            {
                ++flushedElements;

                removeMappedFragments(element.Fragments); // need remove frags first since the processing can remove elements

                mReadGroupHandler.accept(element.Fragments);

                mForwardPositions[mMinPositionIndex] = null;
            }

            mMinPosition++;

            if(mMinPositionIndex + 1 >= mForwardPositions.length)
                mMinPositionIndex = 0;
            else
                ++mMinPositionIndex;
        }

        if(flushCount >= mCapacity)
            resetMinPosition(position);

        if(flushedElements == 0)
            return;

        // flush out any reverse strand position which is now earlier than the current forward strand read start position
        Set<Integer> flushedPositions = Sets.newHashSet();
        for(Map.Entry<Integer,FragmentGroup> entry : mReversePositions.entrySet())
        {
            int reversePosition = abs(entry.getKey());
            if(reversePosition < position)
            {
                flushedPositions.add(entry.getKey());

                removeMappedFragments(entry.getValue().Fragments);
                mReadGroupHandler.accept(entry.getValue().Fragments);
            }
        }

        flushedPositions.forEach(x -> mReversePositions.remove(x));

        checkFragmentLog();
    }

    private void removeMappedFragments(final List<Fragment> fragments)
    {
        for(Fragment fragment : fragments)
        {
            if(fragment.hasLocalMate())
            {
                if(fragment.readCount() == 1)
                    ++mFragmemtUnmatchedCount;

                mFragments.remove(fragment.id());
            }
        }

        // fragments.stream().filter(x -> x.hasLocalMate()).forEach(x -> mFragments.remove(x.id()));
    }

    public void evictAll()
    {
        if(mFragmemtNoCacheCount > LOG_FRAG_COUNT || mFragmemtCacheCount > LOG_FRAG_COUNT)
        {
            MD_LOGGER.debug("read cache eviction: chr({}:{}) fragments(forward={} reverse={}) cache(none={} cache={} unmatched={})",
                    mChromosome, mMinPosition, mFragments.size(), mReversePositions.size(),
                    mFragmemtNoCacheCount, mFragmemtCacheCount, mFragmemtUnmatchedCount);
        }

        for(int i = 0; i < mCapacity; i++)
        {
            FragmentGroup element = mForwardPositions[i];

            // clear and process each element and depth
            if(element != null)
            {
                mReadGroupHandler.accept(element.Fragments);
                mForwardPositions[i] = null;
            }
        }

        mReversePositions.values().forEach(x -> mReadGroupHandler.accept(x.Fragments));

        mReversePositions.clear();
        mFragments.clear();
        mLastLogReadCount = 0;
        mLastFragmentLogCount = 0;
        mFragmemtCacheCount = 0;
        mFragmemtUnmatchedCount = 0;
        mFragmemtNoCacheCount = 0;
    }

    private void resetMinPosition(int position)
    {
        mMinPositionIndex = 0;
        mMinPosition = max(1, position - (int)round(mCapacity * 0.5));
    }

    private static final int LOG_FRAG_COUNT = 10000;

    private void checkFragmentLog()
    {
        if(mLastLogReadCount < LOG_FRAG_COUNT)
            return;

        mLastLogReadCount = 0;

        int forwardPositions = 0;
        int forwardFrags = 0;

        for(int i = 0; i < mCapacity; i++)
        {
            if(mForwardPositions[i] != null)
            {
                forwardFrags += mForwardPositions[i].Fragments.size();
                ++forwardPositions;
            }
        }

        int reverseFrags = mReversePositions.values().stream().mapToInt(x -> x.Fragments.size()).sum();

        int fragmentSize = forwardFrags + reverseFrags;

        if(abs(fragmentSize - mLastFragmentLogCount) < LOG_FRAG_COUNT)
            return;

        mLastFragmentLogCount = fragmentSize;

        MD_LOGGER.debug("read cache: chr({} minPos={}) fragments({}) forward({} frags={}) reverse({} frags={})",
                mChromosome, mMinPosition, fragmentSize, forwardPositions, forwardFrags, mReversePositions.size(), reverseFrags);
    }
}
