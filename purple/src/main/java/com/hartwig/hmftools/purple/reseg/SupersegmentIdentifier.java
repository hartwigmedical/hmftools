package com.hartwig.hmftools.purple.reseg;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.purple.region.ObservedRegion;

public final class SupersegmentIdentifier
{
    // identifies the 'supersegments' - contiguous runs of SV-unsupported diploid segments,
    // potentially bridging over isolated skippable germline deletions or short SVs
    private SupersegmentIdentifier() {}

    public static SupersegmentResult identify(final List<ObservedRegion> segmentsInOrder)
    {
        int n = segmentsInOrder.size();

        boolean[] bothNone = new boolean[n];
        boolean[] canBeSkipped = new boolean[n];

        for(int i = 0; i < n; i++)
        {
            ObservedRegion segment = segmentsInOrder.get(i);
            SegmentSupport otherSupport = (i + 1 < n) ? segmentsInOrder.get(i + 1).support() : null;

            bothNone[i] = segment.support() == SegmentSupport.NONE && otherSupport == SegmentSupport.NONE
                    && segment.germlineStatus() == GermlineStatus.DIPLOID;

            boolean isDeletion = segment.germlineStatus() == GermlineStatus.HOM_DELETION
                    || segment.germlineStatus() == GermlineStatus.HET_DELETION;

            boolean isUnknownBridge = segment.germlineStatus() == GermlineStatus.UNKNOWN && otherSupport != null
                    && segment.support() == otherSupport
                    && (segment.support() == SegmentSupport.DEL || segment.support() == SegmentSupport.DUP);

            canBeSkipped[i] = isDeletion || isUnknownBridge;

            /* TODO: check if matters
            if(i > 0 && bothNone[i] && bothNone[i - 1] && !segment.chromosome().equals(segmentsInOrder.get(i - 1).chromosome()))
            {
                PPL_LOGGER.debug("bothNone run spans chromosome boundary at index({}): {} -> {}",
                        i, segmentsInOrder.get(i - 1).chromosome(), segment.chromosome());
            }
            */
        }

        boolean[] rawNewSuperSegment = new boolean[n];

        for(int i = 0; i < n; i++)
        {
            boolean prevBothNone = i > 0 && bothNone[i - 1];
            rawNewSuperSegment[i] = bothNone[i] && !prevBothNone;
        }

        boolean[] skipOver = new boolean[n];

        for(int i = 0; i < n; i++)
        {
            boolean prevBothNone = i > 0 && bothNone[i - 1];
            boolean nextIsRawStart = (i + 1 < n) && rawNewSuperSegment[i + 1];
            skipOver[i] = canBeSkipped[i] && prevBothNone && nextIsRawStart;
        }

        boolean[] newSuperSegment = rawNewSuperSegment.clone();

        for(int i = 0; i < n; i++)
        {
            if(skipOver[i] && i + 1 < n)
                newSuperSegment[i + 1] = false;
        }

        List<Supersegment> supersegments = new ArrayList<>();
        List<ObservedRegion> excludedSegments = new ArrayList<>();

        List<ObservedRegion> currentBothNone = null;
        List<ObservedRegion> currentSkippable = null;

        for(int i = 0; i < n; i++)
        {
            ObservedRegion segment = segmentsInOrder.get(i);

            if(newSuperSegment[i])
            {
                if(currentBothNone != null)
                    supersegments.add(new Supersegment(supersegments.size(), currentBothNone, currentSkippable));

                currentBothNone = new ArrayList<>();
                currentSkippable = new ArrayList<>();
            }

            if(bothNone[i])
            {
                currentBothNone.add(segment);
            }
            else
            {
                excludedSegments.add(segment);

                if(skipOver[i])
                    currentSkippable.add(segment);
            }
        }

        if(currentBothNone != null)
            supersegments.add(new Supersegment(supersegments.size(), currentBothNone, currentSkippable));

        return new SupersegmentResult(supersegments, excludedSegments);
    }
}
