package com.hartwig.hmftools.sage.filter;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractUmiType;

import com.hartwig.hmftools.common.bam.UmiReadType;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.evidence.RawContext;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.sync.FragmentData;

import htsjdk.samtools.SAMRecord;

public class StrandBiasData
{
    private final boolean mIsAltBias;
    private int mForwardCount;
    private int mReverseCount;

    public StrandBiasData(final boolean isAltBias)
    {
        mIsAltBias = isAltBias;
        mForwardCount = 0;
        mReverseCount = 0;
    }

    public void add(boolean isForward)
    {
        if(isForward)
            ++mForwardCount;
        else
            ++mReverseCount;
    }

    public int forward() { return mForwardCount; }
    public int reverse() { return mReverseCount; }
    public int depth() { return mForwardCount + mReverseCount; }

    public double bias()
    {
        double depth = mForwardCount + mReverseCount;
        return depth > 0 ? mForwardCount / depth : 0.5;
    }

    public double minBias()
    {
        double bias = bias();
        return min(bias, 1 - bias);
    }

    public void registerFragment(final SAMRecord record)
    {
        boolean readIsForward = !record.getReadNegativeStrandFlag();

        if(!record.getReadPairedFlag())
        {
            add(readIsForward);
            return;
        }

        if(extractUmiType(record) == UmiReadType.DUAL)
        {
            ++mForwardCount;
            ++mReverseCount;
            return;
        }

        // make the distinction between F1R2 and F2R1
        boolean firstIsForward = record.getFirstOfPairFlag() ? readIsForward : !record.getMateNegativeStrandFlag();
        boolean secondIsForward = !record.getFirstOfPairFlag() ? readIsForward : !record.getMateNegativeStrandFlag();

        if(firstIsForward != secondIsForward)
        {
            add(firstIsForward);
        }
    }

    private void registerRead(final SAMRecord record)
    {
        add(!record.getReadNegativeStrandFlag());

        /*
        SG_LOGGER.debug("read({}) coords({}-{}) on {}} strand",
                record.getReadName(), record.getAlignmentStart(), record.getAlignmentEnd(),
                record.getReadNegativeStrandFlag() ? "reverse" : "forward");
        */
    }

    public void registerRead(
            final SAMRecord record, final FragmentData fragment, final ReadContextCounter rcCounter)
    {
        // no fragment sync - take the read's strand
        if(fragment == null)
        {
            registerRead(record);
            return;
        }

        if(fragment.First.getReadNegativeStrandFlag() == fragment.Second.getReadNegativeStrandFlag())
            return; // ignore the inverted case

        // read falls only within one or the other read's bounds - use that read's strand
        boolean withinFirst = positionWithin(rcCounter.position(), fragment.First.getAlignmentStart(), fragment.First.getAlignmentEnd());
        boolean withinSecond = positionWithin(rcCounter.position(), fragment.Second.getAlignmentStart(), fragment.Second.getAlignmentEnd());

        /*
        SG_LOGGER.debug("read({}) var({}) withinRead(first={} second={})",
                record.getReadName(), variant.position(), withinFirst, withinSecond);
        */

        SimpleVariant variant = rcCounter.variant();

        if(withinFirst && !withinSecond)
        {
            ReadContextMatch matchType = determineReadContextMatch(rcCounter, fragment.First);

            if(hasSupport(matchType))
                registerRead(fragment.First);

            return;
        }
        else if(!withinFirst && withinSecond)
        {
            ReadContextMatch matchType = determineReadContextMatch(rcCounter, fragment.Second);

            if(hasSupport(matchType))
                registerRead(fragment.Second);

            return;
        }

        // look at the raw alt support from each read to determine which to count
        ReadContextMatch firstMatchType = determineReadContextMatch(rcCounter, fragment.First);
        ReadContextMatch secondMatchType = determineReadContextMatch(rcCounter, fragment.Second);

        if(hasSupport(firstMatchType))
            registerRead(fragment.First);

        if(hasSupport(secondMatchType))
            registerRead(fragment.Second);

        /*
        SG_LOGGER.debug("read({}) var({}) isAltBias({}) altSupport(first={} second={})",
                record.getReadName(), variant.position(), mIsAltBias, firstMatchType, secondMatchType);
        */
    }

    private ReadContextMatch determineReadContextMatch(final ReadContextCounter rcCounter, final SAMRecord read)
    {
        RawContext rawContext = RawContext.createFromRead(rcCounter.variant(), read);

        if(rawContext.ReadVariantIndex < 0)
            return ReadContextMatch.NONE;

        if(!rcCounter.matcher().coversVariant(read, rawContext.ReadVariantIndex))
            return ReadContextMatch.NONE;

        return rcCounter.matcher().determineReadMatch(read, rawContext.ReadVariantIndex);
    }

    private boolean hasSupport(final ReadContextMatch matchType)
    {
        return mIsAltBias ? matchType.SupportsAlt : matchType.SupportsRef;
    }

    public String toString() { return format("fwd=%d rev=%d total=%d bias=%.3f", mForwardCount, mReverseCount, depth(), bias()); }
}
