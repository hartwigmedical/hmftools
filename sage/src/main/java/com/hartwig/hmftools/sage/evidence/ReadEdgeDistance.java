package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.samtools.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.samtools.CigarUtils.rightSoftClipLength;

import com.hartwig.hmftools.sage.sync.FragmentData;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ReadEdgeDistance
{
    private final int mVariantPosition;
    
    private int mMaxDistanceFromUnclippedEdge; // includes any soft-clipped bases
    private int mMaxDistanceFromUnclippedEdgeAlt;
    private int mMaxDistanceFromAlignedEdgeAlt;
    
    public ReadEdgeDistance(final int variantPosition)
    {
        mVariantPosition = variantPosition;
        mMaxDistanceFromUnclippedEdgeAlt = 0;
        mMaxDistanceFromUnclippedEdge = 0;
        mMaxDistanceFromAlignedEdgeAlt = 0;
    }

    public static int calcAdjustedVariantPosition(final int variantPosition, final int indelLength)
    {
        // INDELs measure MED from their mid-point
        return indelLength == 0 ? variantPosition : variantPosition +  indelLength / 2;
    }

    public int maxAltDistanceFromAlignedEdge() { return mMaxDistanceFromAlignedEdgeAlt; }
    public int maxDistanceFromUnclippedEdge() { return mMaxDistanceFromUnclippedEdge; }
    public int maxAltDistanceFromUnclippedEdge() { return mMaxDistanceFromUnclippedEdgeAlt; }

    public void update(final SAMRecord record, final FragmentData fragmentData, boolean altSupport)
    {
        if(mMaxDistanceFromUnclippedEdgeAlt >= record.getReadBases().length / 2)
            return;

        if(record.getCigar().containsOperator(CigarOperator.N)) // unnecessary in append mode
            return;

        // determine how far from the edge of the read the variant is, ignoring soft-clip bases and realigned reads
        // and for INDELs use the position of the middle base of INDEL
        // take the lower of the 2 distances for each read, eg 50 and 100 bases, then take 50

        int minDistance = NO_READ_EDGE_DISTANCE;
        int minDistanceUnclipped = NO_READ_EDGE_DISTANCE;

        if(fragmentData != null)
        {
            int minDistanceFirst = calcDistanceFromReadEdge(mVariantPosition, fragmentData.First, false);
            int minDistanceSecond = calcDistanceFromReadEdge(mVariantPosition, fragmentData.Second, false);

            if(minDistanceFirst > NO_READ_EDGE_DISTANCE && minDistanceSecond > NO_READ_EDGE_DISTANCE)
                minDistance = min(minDistanceFirst, minDistanceSecond);
            else if(minDistanceFirst > NO_READ_EDGE_DISTANCE)
                minDistance = minDistanceFirst;
            else if(minDistanceSecond > NO_READ_EDGE_DISTANCE)
                minDistance = minDistanceSecond;

            // again for the unclipped distances
            minDistanceFirst = calcDistanceFromReadEdge(mVariantPosition, fragmentData.First, true);
            minDistanceSecond = calcDistanceFromReadEdge(mVariantPosition, fragmentData.Second, true);

            if(minDistanceFirst > NO_READ_EDGE_DISTANCE && minDistanceSecond > NO_READ_EDGE_DISTANCE)
                minDistanceUnclipped = min(minDistanceFirst, minDistanceSecond);
            else if(minDistanceFirst > NO_READ_EDGE_DISTANCE)
                minDistanceUnclipped = minDistanceFirst;
            else if(minDistanceSecond > NO_READ_EDGE_DISTANCE)
                minDistanceUnclipped = minDistanceSecond;
        }
        else
        {
            minDistance = calcDistanceFromReadEdge(mVariantPosition, record, false);
            minDistanceUnclipped = calcDistanceFromReadEdge(mVariantPosition, record, true);
        }

        mMaxDistanceFromUnclippedEdge = max(minDistanceUnclipped, mMaxDistanceFromUnclippedEdge);

        if(altSupport)
        {
            mMaxDistanceFromAlignedEdgeAlt = max(minDistance, mMaxDistanceFromAlignedEdgeAlt);
            mMaxDistanceFromUnclippedEdgeAlt = max(minDistanceUnclipped, mMaxDistanceFromUnclippedEdgeAlt);
        }
    }

    private static final int NO_READ_EDGE_DISTANCE = -1;

    private static int calcDistanceFromReadEdge(final int variantPosition, final SAMRecord record, boolean includeSoftClipBases)
    {
        if(!positionWithin(variantPosition, record.getAlignmentStart(), record.getAlignmentEnd()))
            return NO_READ_EDGE_DISTANCE;

        int distFromStart = variantPosition - record.getAlignmentStart();
        int distFromEnd = record.getAlignmentEnd() - variantPosition;

        if(includeSoftClipBases)
        {
            distFromStart += leftSoftClipLength(record);
            distFromEnd += rightSoftClipLength(record);
        }

        return min(distFromStart, distFromEnd);
    }
}
