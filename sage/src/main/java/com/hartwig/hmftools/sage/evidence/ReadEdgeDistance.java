package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;

import com.hartwig.hmftools.sage.sync.FragmentData;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ReadEdgeDistance
{
    private final int mVariantPosition;

    // all distances are recorded as a percentage of the read length from which they were taken

    private double mMaxDistanceFromEdge; // includes any soft-clipped bases
    private double mMaxDistanceFromEdgeAlt;

    // captures edge distance percentage in a sided-way, where 0.1 is 10% from left edge, nad 0.9 is 10% from the right
    private double mTotalDistanceFromEdge;
    private double mTotalDistanceFromEdgeAlt;
    private int mUpdates;
    private int mUpdatesAlt;

    public ReadEdgeDistance(final int variantPosition)
    {
        mVariantPosition = variantPosition;
        mMaxDistanceFromEdge = 0;
        mMaxDistanceFromEdgeAlt = 0;
        mTotalDistanceFromEdge = 0;
        mTotalDistanceFromEdgeAlt = 0;
        mUpdates = 0;
        mUpdatesAlt = 0;
    }

    public static int calcAdjustedVariantPosition(final int variantPosition, final int indelLength)
    {
        // INDELs measure MED from their mid-point
        return indelLength == 0 ? variantPosition : variantPosition +  indelLength / 2;
    }

    public double maxAltDistanceFromEdge() { return mMaxDistanceFromEdgeAlt; }
    public double maxDistanceFromEdge() { return mMaxDistanceFromEdge; }

    public double avgDistanceFromEdge() { return calcAdjustedAvgEdge(false); }
    public double avgAltDistanceFromEdge() { return calcAdjustedAvgEdge(true); }

    private double calcAdjustedAvgEdge(boolean useAlt)
    {
        double totalDistance = useAlt ? mTotalDistanceFromEdgeAlt : mTotalDistanceFromEdge;
        int updates = useAlt ? mUpdatesAlt : mUpdates;
        double average = updates > 0 ? totalDistance / updates : 0;
        return average <= 0.5 ? average : 1 - average;
    }

    public void update(final SAMRecord record, final FragmentData fragmentData, boolean altSupport)
    {
        if(record.getCigar().containsOperator(CigarOperator.N)) // unnecessary in append mode
            return;

        // determine how far from the edge of the read the variant is, ignoring soft-clip bases and realigned reads
        // and for INDELs use the position of the middle base of INDEL
        // take the lowest of the 2 distances for each read, eg 50 and 100 bases, then take 50

        Integer minDistance = null;

        if(fragmentData != null)
        {
            Integer minDistanceFirst = calcDistanceFromReadEdge(mVariantPosition, fragmentData.First);
            Integer minDistanceSecond = calcDistanceFromReadEdge(mVariantPosition, fragmentData.Second);

            // use the fragment bounds if the variant falls within the overlapping section, otherwise just the read it relates to
            if(minDistanceFirst != null && minDistanceSecond != null)
                minDistance = calcDistanceFromReadEdge(mVariantPosition, record);
            else if(minDistanceFirst != null)
                minDistance = minDistanceFirst;
            else if(minDistanceSecond != null)
                minDistance = minDistanceSecond;
        }
        else
        {
            minDistance = calcDistanceFromReadEdge(mVariantPosition, record);
        }

        if(minDistance == null)
            return;

        double readLengthPerc = minDistance / (double)record.getReadLength();
        double readLengthPercAbs = abs(readLengthPerc);
        double readLengthPercSided = readLengthPerc > 0 ? readLengthPerc : 1 - readLengthPercAbs;

        ++mUpdates;
        mMaxDistanceFromEdge = max(readLengthPercAbs, mMaxDistanceFromEdge);

        if(altSupport)
        {
            mMaxDistanceFromEdgeAlt = max(readLengthPercAbs, mMaxDistanceFromEdgeAlt);
            mTotalDistanceFromEdgeAlt += readLengthPercSided;
            ++mUpdatesAlt;
        }

        mTotalDistanceFromEdge += readLengthPercSided;
    }

    private static Integer calcDistanceFromReadEdge(final int variantPosition, final SAMRecord record)
    {
        // returns null if invalid, otherwise +ve distance if closer to the left, or -ve distance if closer to right read edge
        if(!positionWithin(variantPosition, record.getAlignmentStart(), record.getAlignmentEnd()))
            return null;

        int distFromStart = variantPosition - record.getAlignmentStart();
        int distFromEnd = record.getAlignmentEnd() - variantPosition;

        return distFromStart < distFromEnd ? distFromStart : -distFromEnd;
    }
}
