package com.hartwig.hmftools.sage.evidence;

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

    public double avgDistanceFromEdge() { return mUpdates > 0 ? mTotalDistanceFromEdge / mUpdates : 0; }
    public double avgAltDistanceFromEdge() { return mUpdatesAlt > 0 ? mTotalDistanceFromEdgeAlt / mUpdatesAlt : 0; }

    public void update(final SAMRecord record, final FragmentData fragmentData, boolean altSupport)
    {
        if(record.getCigar().containsOperator(CigarOperator.N)) // unnecessary in append mode
            return;

        // determine how far from the edge of the read the variant is, ignoring soft-clip bases and realigned reads
        // and for INDELs use the position of the middle base of INDEL
        // take the lowest of the 2 distances for each read, eg 50 and 100 bases, then take 50

        int minDistance = NO_READ_EDGE_DISTANCE;

        if(fragmentData != null)
        {
            int minDistanceFirst = calcDistanceFromReadEdge(mVariantPosition, fragmentData.First);
            int minDistanceSecond = calcDistanceFromReadEdge(mVariantPosition, fragmentData.Second);

            // use the fragment bounds if the variant falls within the overlapping section, otherwise just the read it relates to
            if(minDistanceFirst > NO_READ_EDGE_DISTANCE && minDistanceSecond > NO_READ_EDGE_DISTANCE)
                minDistance = calcDistanceFromReadEdge(mVariantPosition, record);
            else if(minDistanceFirst > NO_READ_EDGE_DISTANCE)
                minDistance = minDistanceFirst;
            else if(minDistanceSecond > NO_READ_EDGE_DISTANCE)
                minDistance = minDistanceSecond;
        }
        else
        {
            minDistance = calcDistanceFromReadEdge(mVariantPosition, record);
        }

        double readLengthPercentage = minDistance / (double)record.getReadLength();
        ++mUpdates;
        mMaxDistanceFromEdge = max(readLengthPercentage, mMaxDistanceFromEdge);

        if(altSupport)
        {
            mMaxDistanceFromEdgeAlt = max(readLengthPercentage, mMaxDistanceFromEdgeAlt);
            mTotalDistanceFromEdgeAlt += readLengthPercentage;
            ++mUpdatesAlt;
        }

        mTotalDistanceFromEdge += readLengthPercentage;
    }

    private static final int NO_READ_EDGE_DISTANCE = -1;

    private static int calcDistanceFromReadEdge(final int variantPosition, final SAMRecord record)
    {
        if(!positionWithin(variantPosition, record.getAlignmentStart(), record.getAlignmentEnd()))
            return NO_READ_EDGE_DISTANCE;

        int distFromStart = variantPosition - record.getAlignmentStart();
        int distFromEnd = record.getAlignmentEnd() - variantPosition;

        return min(distFromStart, distFromEnd);
    }
}
