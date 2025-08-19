package com.hartwig.hmftools.panelbuilder;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

// Defines the region(s) and base sequence of a probe.
// This may be:
//   - A single reference genome region (the case for most probes)
//       probe = start
//   - 1 or 2 reference genome regions and a custom insert sequence (the case for variant probes).
//       probe = start + insert + end
public record SequenceDefinition(
        @Nullable ChrBaseRegion startRegion,
        // If REVERSE then start region is reverse complemented.
        @Nullable Orientation startOrientation,
        @Nullable String insertSequence,
        @Nullable ChrBaseRegion endRegion,
        // If REVERSE then end region is reverse complemented.
        @Nullable Orientation endOrientation
)
{
    public SequenceDefinition
    {
        boolean valid1 = startRegion != null && startOrientation == null && insertSequence == null && endRegion == null;
        boolean valid2 = insertSequence != null && (startRegion != null || endRegion != null);
        if(!(valid1 || valid2))
        {
            throw new IllegalArgumentException();
        }
        if(startRegion != null && !startRegion.hasValidPositions())
        {
            throw new IllegalArgumentException();
        }
        if(endRegion != null && !endRegion.hasValidPositions())
        {
            throw new IllegalArgumentException();
        }
        if(startRegion == null && startOrientation != null)
        {
            throw new IllegalArgumentException();
        }
        if(endRegion == null && endOrientation != null)
        {
            throw new IllegalArgumentException();
        }
    }

    public static SequenceDefinition exactRegion(final ChrBaseRegion region)
    {
        return new SequenceDefinition(region, null, null, null, null);
    }

    public static SequenceDefinition simpleMutation(final ChrBaseRegion startRegion, final String insertSequence,
            final ChrBaseRegion endRegion)
    {
        Orientation startOrientation = startRegion == null ? null : Orientation.FORWARD;
        Orientation endOrientation = endRegion == null ? null : Orientation.FORWARD;
        return new SequenceDefinition(startRegion, startOrientation, insertSequence, endRegion, endOrientation);
    }

    public static SequenceDefinition structuralVariant(final ChrBaseRegion startRegion, final Orientation startOrientation,
            final String insertSequence, final ChrBaseRegion endRegion, final Orientation endOrientation)
    {
        return new SequenceDefinition(startRegion, startOrientation, insertSequence, endRegion, endOrientation);
    }

    public List<ChrBaseRegion> regions()
    {
        List<ChrBaseRegion> result = new ArrayList<>(2);
        if(startRegion != null)
        {
            result.add(startRegion);
        }
        if(endRegion != null)
        {
            result.add(endRegion);
        }
        return result;
    }

    public int baseLength()
    {
        int length = 0;
        if(startRegion != null)
        {
            length += startRegion.baseLength();
        }
        if(insertSequence != null)
        {
            length += insertSequence.length();
        }
        if(endRegion != null)
        {
            length += endRegion.baseLength();
        }
        return length;
    }

    // Checks if the probe is defined by a single region.
    public boolean isExactRegion()
    {
        return startRegion != null && insertSequence == null && endRegion == null;
    }

    // Gets the single region that the probe is defined by, or throws an exception.
    public ChrBaseRegion exactRegion()
    {
        ChrBaseRegion region = exactRegionOrNull();
        if(region == null)
        {
            throw new IllegalArgumentException("Probe has multiple regions");
        }
        else
        {
            return region;
        }
    }

    @Nullable
    public ChrBaseRegion exactRegionOrNull()
    {
        return isExactRegion() ? startRegion : null;
    }
}
