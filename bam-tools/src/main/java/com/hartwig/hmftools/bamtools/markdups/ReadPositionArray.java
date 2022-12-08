package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;

import java.util.function.Consumer;

import com.hartwig.hmftools.bamtools.ReadGroup;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

public class ReadPositionArray
{
    // a ring buffer to store reads at each read starting position
    private final String mChromosome;
    private final PositionReadGroups[] mElements;
    private final Consumer<PositionReadGroups> mReadGroupHandler;
    private int mMinPosition;
    private int mMinPositionIndex;

    private final int mCapacity;

    public ReadPositionArray(final String chromosome, int capacity, final Consumer<PositionReadGroups> evictionHandler)
    {
        mChromosome = chromosome;
        mReadGroupHandler = evictionHandler;
        mCapacity = capacity;
        mElements = new PositionReadGroups[mCapacity];
        mMinPosition = 0;
        mMinPositionIndex = 0;
    }

    public boolean withinRange(int position) { return position >= mMinPosition && position < mMinPosition + mCapacity; }

    public boolean processRead(final SAMRecord record)
    {
        int readPosStart = record.getAlignmentStart();

        if(mMinPosition == 0)
            mMinPosition = readPosStart;

        if(!isValidPosition(readPosStart))
            return false;

        boolean readStored = handleRead(record);

        checkFlush(readPosStart);

        return readStored;
    }

    private boolean handleRead(final SAMRecord record)
    {
        /* scenarios:
            - unpaired or mate is unmapped
                - store as potential duplicate
                - supplementaries - store if within range, to avoid processing later
            - read's mate is on another chromosome - store since can still evaluate if duplicate (from both chromosomes)
            - lower pos with mate higher on same chromosome - store regardless
            - higher pos with mate lower on same chromosome
                - mate within current range - look up and add
                - mate already flushed, return unhandled
            - supplementary
                - with mate within current range - look up and add
                - otherwise return unhandled
        */
        SupplementaryReadData suppData = record.getSupplementaryAlignmentFlag() ? SupplementaryReadData.from(record) : null;

        if(!record.getReadPairedFlag() || record.getMateUnmappedFlag())
        {
            // could be a secondary or supplementary, otherwise store
            if(record.getSupplementaryAlignmentFlag())
            {
                if(suppData.Chromosome.equals(mChromosome) && withinRange(suppData.Position))
                {
                    return storeAdditionalRead(record, suppData.Position);
                }
                else
                {
                    return false;
                }
            }

            storeLowerRead(record);
            return true;
        }

        int matePosStart = record.getMateAlignmentStart();
        String mateStr = record.getMateReferenceName();

        if(record.getSupplementaryAlignmentFlag())
        {
            if(mateStr.equals(mChromosome) && withinRange(matePosStart))
            {
                return storeAdditionalRead(record, matePosStart);
            }
            else if(suppData.Chromosome.equals(mChromosome) && withinRange(suppData.Position))
            {
                return storeAdditionalRead(record, suppData.Position);
            }
            else
            {
                return false;
            }
        }

        // mate is elsewhere so store this primary read
        if(!mateStr.equals(mChromosome))
        {
            storeLowerRead(record);
            return true;
        }

        int readPosStart = record.getAlignmentStart();
        if(readPosStart < matePosStart)
        {
            storeLowerRead(record);
            return true;
        }

        // higher of 2 reads on the same chromosome
        if(withinRange(matePosStart))
        {
            return storeAdditionalRead(record, matePosStart);
        }
        else
        {
            // will be linked via the group combiner cache
            return false;
        }
    }

    private void storeLowerRead(final SAMRecord record)
    {
        int index = calcIndex(record.getAlignmentStart());

        PositionReadGroups element = mElements[index];
        if(element == null)
        {
            element = new PositionReadGroups(record);
            mElements[index] = element;
        }
        else
        {
            element.ReadGroups.add(new ReadGroup(record));
        }
    }

    private boolean storeAdditionalRead(final SAMRecord record, int lowerPosStart)
    {
        int index = calcIndex(lowerPosStart);

        PositionReadGroups element = mElements[index];
        if(element != null)
        {
            ReadGroup readGroup = element.ReadGroups.stream().filter(x -> x.id().equals(record.getReadName())).findFirst().orElse(null);

            if(readGroup != null)
            {
                readGroup.addRead(record);
                return true;
            }
        }

        return false;
    }

    public void evictAll()
    {
        checkFlush(-1);
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

    private boolean isValidPosition(int position)
    {
        if(mMinPosition > 0 && position < mMinPosition)
        {
            BM_LOGGER.warn("ignoring read with position({}) before prior position({})", position, mMinPosition);
            return false;
        }

        return true;
    }

    private void checkFlush(int position)
    {
        int flushCount = 0;

        if(position > 0)
        {
            if(mMinPosition == 0)
            {
                resetMinPosition(position);
                return;
            }

            int distanceFromMinPosition = position - mMinPosition;

            if(distanceFromMinPosition < mCapacity)
                return;

            flushCount = position - mMinPosition - mCapacity + 1;
        }
        else
        {
            flushCount = mCapacity;
        }

        // only iterate at most once through the array
        for(int i = 0; i < min(flushCount, mCapacity); i++)
        {
            PositionReadGroups element = mElements[mMinPositionIndex];

            // clear and process each element and depth
            if(element != null)
            {
                mReadGroupHandler.accept(element);
                mElements[mMinPositionIndex] = null;
            }

            mMinPosition++;

            if(mMinPositionIndex + 1 >= mElements.length)
                mMinPositionIndex = 0;
            else
                ++mMinPositionIndex;
        }

        if(flushCount >= mCapacity)
            resetMinPosition(position);
    }

    private void resetMinPosition(int position)
    {
        mMinPositionIndex = 0;
        mMinPosition = max(1, position - (int)round(mCapacity * 0.5));
    }
}
