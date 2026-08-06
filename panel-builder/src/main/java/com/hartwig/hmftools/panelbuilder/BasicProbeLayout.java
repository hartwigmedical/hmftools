package com.hartwig.hmftools.panelbuilder;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

// The start + insert + end decomposition of a probe sequence. Valid only for probes with at most two regions (single region, single
// breakend, or a two-region variant), keeping that assumption out of SequenceDefinition.
public record BasicProbeLayout(
        @Nullable ChrBaseRegion startRegion,
        @Nullable Orientation startOrientation,
        // Empty if there is no insert.
        String insertSequence,
        @Nullable ChrBaseRegion endRegion,
        @Nullable Orientation endOrientation
)
{
    // Matches the sequence definition against each valid basic shape (single region, single breakend, two-region variant) and rejects
    // anything else (e.g. more than two regions, multiple inserts, insert only).
    public static BasicProbeLayout from(final SequenceDefinition definition)
    {
        List<SequenceSegment> segments = definition.segments();
        SequenceSegment first = segments.get(0);

        if(segments.size() == 1 && first instanceof RefSegment only)
        {
            // Single region.
            return new BasicProbeLayout(only.region(), only.orientation(), "", null, null);
        }
        if(segments.size() == 2)
        {
            SequenceSegment second = segments.get(1);
            if(first instanceof RefSegment start && second instanceof InsertSeqSegment insert)
            {
                // Single breakend: region then insert.
                return new BasicProbeLayout(start.region(), start.orientation(), insert.sequence(), null, null);
            }
            if(first instanceof InsertSeqSegment insert && second instanceof RefSegment end)
            {
                // Single breakend: insert then region.
                return new BasicProbeLayout(null, null, insert.sequence(), end.region(), end.orientation());
            }
            if(first instanceof RefSegment start && second instanceof RefSegment end)
            {
                // Two regions, no insert.
                return new BasicProbeLayout(start.region(), start.orientation(), "", end.region(), end.orientation());
            }
        }
        if(segments.size() == 3
                && first instanceof RefSegment start
                && segments.get(1) instanceof InsertSeqSegment insert
                && segments.get(2) instanceof RefSegment end)
        {
            // Two regions with an insert between.
            return new BasicProbeLayout(start.region(), start.orientation(), insert.sequence(), end.region(), end.orientation());
        }
        throw new IllegalArgumentException("Sequence definition is not expressible as a legacy start/insert/end layout: " + definition);
    }
}
