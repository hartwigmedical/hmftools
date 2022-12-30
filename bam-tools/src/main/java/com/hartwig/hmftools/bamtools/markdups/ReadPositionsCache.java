package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentCoordinates.NO_COORDS;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.NONE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.checkDuplicateFragments;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.SAMRecord;

public class ReadPositionsCache
{
    // a ring buffer to store reads at each read starting position
    private final String mChromosome;
    private final CandidateDuplicates[] mForwardPositions;
    private final Map<Integer,CandidateDuplicates> mReversePositions;
    private final Map<String,Fragment> mFragments;
    private final Consumer<List<Fragment>> mReadGroupHandler;
    private final List<Fragment> mResolvedFragments;
    private int mMinPosition;
    private int mMinPositionIndex;
    private int mLastFragmentLogCount;

    private final int mCapacity;

    public ReadPositionsCache(final String chromosome, int capacity, final Consumer<List<Fragment>> evictionHandler)
    {
        mChromosome = chromosome;
        mReadGroupHandler = evictionHandler;
        mCapacity = capacity;
        mForwardPositions = new CandidateDuplicates[mCapacity];
        mReversePositions = Maps.newHashMap();
        mFragments = Maps.newHashMap();
        mResolvedFragments = Lists.newArrayList();
        mMinPosition = 0;
        mMinPositionIndex = 0;
        mLastFragmentLogCount = 0;
    }

    public boolean processRead(final SAMRecord read)
    {
        // supplementaries aren't processed
        // check for an existing mate if on the same chromosome
        // store in a group of fragments with a matching first fragment coordinate
        // check for duplicates with existing fragments, and if a duplicate the don't store and write immediately
        if(!read.getReadPairedFlag() || read.getMateUnmappedFlag())
        {
            storeInitialRead(read);
            return true;
        }

        // mate is elsewhere so store this primary read
        if(!read.getMateReferenceName().equals(mChromosome))
        {
            storeInitialRead(read);
            return true;
        }

        // check for the mate already cached and nearby
        Fragment fragment = mFragments.get(read.getReadName());

        if(fragment != null)
        {
            fragment.addRead(read);
            return true;
        }

        // if the mate was already processed earlier, then could skip checking its duplicate status again
        // but for the sake of UMIs, they still need to be grouped with their duplicate reads
        /*
        if(read.getAlignmentStart() > read.getMateAlignmentStart())
            return false;
        */

        storeInitialRead(read);
        return true;
    }

    private void storeInitialRead(final SAMRecord read)
    {
        Fragment fragment = new Fragment(read);

        if(fragment.coordinates() == NO_COORDS)
        {
            // need to understand why this isn't set for some fragments
            if(fragment.reads().get(0).getDuplicateReadFlag())
                fragment.setStatus(DUPLICATE);
            else
                fragment.setStatus(NONE);

            mResolvedFragments.add(fragment);

            // store this fragment if its mate is in the same partition
            if(read.getMateReferenceName().equals(mChromosome))
                mFragments.put(read.getReadName(), fragment);

            return;
        }

        int fragmentPosition = fragment.initialPosition();

        if(fragmentPosition > 0)
        {
            if(mMinPosition == 0)
                resetMinPosition(fragmentPosition);
            else
                checkFlush(fragmentPosition);
        }

        Fragment duplicateFragment = handleInitialFragment(fragment);

        if(duplicateFragment != null)
            mResolvedFragments.add(duplicateFragment);

        // store this fragment if its mate is in the same partition
        if(read.getMateReferenceName().equals(mChromosome))
        {
            mFragments.put(read.getReadName(), fragment);
        }
    }

    private Fragment handleInitialFragment(final Fragment fragment)
    {
        int fragmentPosition = fragment.initialPosition();

        CandidateDuplicates element = null;

        if(fragmentPosition > 0)
        {
            int index = calcIndex(fragmentPosition);

            if(index < 0 || index >= mCapacity)
            {
                BM_LOGGER.error("fragment({}) outside forward strand array bounds", fragment);
                return null;
            }

            element = mForwardPositions[index];
            if(element == null)
            {
                element = new CandidateDuplicates(fragmentPosition, fragment);
                mForwardPositions[index] = element;
                return null;
            }
        }
        else
        {
            // store in reverse strand map
            element = mReversePositions.get(fragmentPosition);
            if(element == null)
            {
                element = new CandidateDuplicates(fragmentPosition, fragment);
                mReversePositions.put(fragmentPosition, element);
                return null;
            }
        }

        // check new fragment against existing for this position
        return checkDuplicateFragments(fragment, element.Fragments);
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
            CandidateDuplicates element = mForwardPositions[mMinPositionIndex];

            // clear and process each element and depth
            if(element != null)
            {
                ++flushedElements;

                mResolvedFragments.addAll(element.Fragments);

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
        for(Map.Entry<Integer,CandidateDuplicates> entry : mReversePositions.entrySet())
        {
            int reversePosition = abs(entry.getKey());
            if(reversePosition < position)
            {
                flushedPositions.add(entry.getKey());

                List<Fragment> fragments = entry.getValue().Fragments;

                mResolvedFragments.addAll(fragments);
            }
        }

        flushedPositions.forEach(x -> mReversePositions.remove(x));

        mResolvedFragments.forEach(x -> mFragments.remove(x.id()));

        mReadGroupHandler.accept(mResolvedFragments);

        mResolvedFragments.clear();

        checkFragmentLog();
    }

    public void evictAll()
    {
        for(int i = 0; i < mCapacity; i++)
        {
            CandidateDuplicates element = mForwardPositions[i];

            // clear and process each element and depth
            if(element != null)
            {
                mResolvedFragments.addAll(element.Fragments);
                mForwardPositions[i] = null;
            }
        }

        mReversePositions.values().forEach(x -> mResolvedFragments.addAll(x.Fragments));

        mReadGroupHandler.accept(mResolvedFragments);

        mResolvedFragments.clear();
        mReversePositions.clear();
        mFragments.clear();
    }

    private void resetMinPosition(int position)
    {
        mMinPositionIndex = 0;
        mMinPosition = max(1, position - (int)round(mCapacity * 0.5));
    }

    private static final int LOG_FRAG_COUNT = 10000;

    private void checkFragmentLog()
    {
        if(abs(mFragments.size() - mLastFragmentLogCount) < LOG_FRAG_COUNT)
            return;

        mLastFragmentLogCount = mFragments.size();
        BM_LOGGER.debug("{}", cacheStatsStr());
    }

    public String cacheStatsStr()
    {
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

        return format("read cache: fragments(%d) forward(%d frags=%d) reverse(%d frags=%d) minPosition(%d)",
                mFragments.size(), forwardPositions, forwardFrags, mReversePositions.size(), reverseFrags, mMinPosition);
    }
}
