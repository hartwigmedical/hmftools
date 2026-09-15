package com.hartwig.hmftools.panelbuilder;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

// Defines the region(s) and base sequence of a probe as an ordered list of segments contiguous in probe-sequence space.
// A segment is either a reference genome region (optionally reverse complemented) or a literal insert sequence.
// Common shapes:
//   - A single reference genome region (most probes).
//   - 1-2 regions with an optional insert (variant probes).
//   - Multiple disjoint regions (spliced probes, e.g. RNA across exon junctions).
public record SequenceDefinition(List<SequenceSegment> segments) implements Comparable<SequenceDefinition>
{
    public SequenceDefinition
    {
        segments = List.copyOf(segments);

        if(segments.isEmpty())
        {
            throw new IllegalArgumentException("Sequence definition has no segments");
        }
        if(segments.stream().noneMatch(segment -> segment instanceof RefSegment))
        {
            throw new IllegalArgumentException("Sequence definition must contain at least one region");
        }
        for(int i = 1; i < segments.size(); ++i)
        {
            SequenceSegment prev = segments.get(i - 1);
            SequenceSegment next = segments.get(i);
            // Consecutive inserts should have been a single insert segment.
            if(prev instanceof InsertSeqSegment && next instanceof InsertSeqSegment)
            {
                throw new IllegalArgumentException("Consecutive insert segments should be a single segment");
            }
            // Consecutive regions that are directly adjacent in the genome with the same orientation should have been a single region.
            if(prev instanceof RefSegment prevRegion && next instanceof RefSegment nextRegion && prevRegion.isGenomeAdjacentTo(nextRegion))
            {
                throw new IllegalArgumentException("Adjacent regions with the same orientation should be a single region");
            }
        }
    }

    public static SequenceDefinition singleRegion(final ChrBaseRegion region)
    {
        return new SequenceDefinition(List.of(new RefSegment(region, Orientation.FORWARD)));
    }

    // Single genome breakend with a novel insert sequence following it.
    public static SequenceDefinition forwardSgl(final ChrBaseRegion startRegion, final String insertSequence)
    {
        return new SequenceDefinition(List.of(new RefSegment(startRegion, Orientation.FORWARD), new InsertSeqSegment(insertSequence)));
    }

    // Novel insert sequence followed by a single genome breakend.
    public static SequenceDefinition reverseSgl(final String insertSequence, final ChrBaseRegion endRegion)
    {
        return new SequenceDefinition(List.of(new InsertSeqSegment(insertSequence), new RefSegment(endRegion, Orientation.FORWARD)));
    }

    // Two genome regions with an optional insert between them (SNV/INDEL/SV probes).
    public static SequenceDefinition variant(final ChrBaseRegion startRegion, final Orientation startOrientation,
            final String insertSequence, final ChrBaseRegion endRegion, final Orientation endOrientation)
    {
        List<SequenceSegment> segments = insertSequence.isEmpty()
                ? List.of(new RefSegment(startRegion, startOrientation), new RefSegment(endRegion, endOrientation))
                : List.of(
                        new RefSegment(startRegion, startOrientation), new InsertSeqSegment(insertSequence),
                        new RefSegment(endRegion, endOrientation));
        return new SequenceDefinition(segments);
    }

    // A probe spanning one or more genome regions, all genome-forward (e.g. a spliced RNA probe from a region mapping).
    // The regions must not be directly adjacent with the same orientation (they should already be merged); see the canonical constructor.
    public static SequenceDefinition spliced(final List<ChrBaseRegion> regions)
    {
        List<SequenceSegment> segments = regions.stream()
                .map(region -> (SequenceSegment) new RefSegment(region, Orientation.FORWARD))
                .toList();
        return new SequenceDefinition(segments);
    }

    public List<ChrBaseRegion> regions()
    {
        return segments.stream()
                .filter(segment -> segment instanceof RefSegment)
                .map(segment -> ((RefSegment) segment).region())
                .toList();
    }

    public int baseLength()
    {
        return segments.stream().mapToInt(SequenceSegment::baseLength).sum();
    }

    // Checks if the probe is defined by a single reference genome region.
    public boolean isSingleRegion()
    {
        return segments.size() == 1 && segments.get(0) instanceof RefSegment;
    }

    // Checks if the probe spans more than one disjoint reference genome region (e.g. SV probe or spliced exons probe).
    public boolean isMultiRegion()
    {
        return regions().size() > 1;
    }

    // Gets the single region that the probe is defined by, or throws an exception.
    public ChrBaseRegion singleRegion()
    {
        ChrBaseRegion region = singleRegionOrNull();
        if(region == null)
        {
            throw new IllegalArgumentException("Probe has multiple regions");
        }
        return region;
    }

    @Nullable
    public ChrBaseRegion singleRegionOrNull()
    {
        return isSingleRegion() ? ((RefSegment) segments.get(0)).region() : null;
    }

    @Override
    public int compareTo(final SequenceDefinition other)
    {
        // Consistent structural ordering for output and debugging purposes.
        int common = Math.min(segments.size(), other.segments.size());
        for(int i = 0; i < common; ++i)
        {
            int compare = segments.get(i).compareTo(other.segments.get(i));
            if(compare != 0)
            {
                return compare;
            }
        }
        return Integer.compare(segments.size(), other.segments.size());
    }
}
