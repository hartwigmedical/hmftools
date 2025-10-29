package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.SequenceUtils.isDnaSequenceNormal;

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
        String insertSequence,
        @Nullable ChrBaseRegion endRegion,
        // If REVERSE then end region is reverse complemented.
        @Nullable Orientation endOrientation
)
{
    public SequenceDefinition
    {
        // Simple single region.
        boolean single =
                startRegion != null && startOrientation == null && insertSequence.isEmpty() && endRegion == null && endOrientation == null;
        // SNV/INDEL/SV. start (+ insert) + end
        boolean genericVariant =
                startRegion != null && startOrientation != null && endRegion != null && endOrientation != null
                        // Ensure the regions don't join each other; otherwise that should be constructed as a single region.
                        && (!insertSequence.isEmpty() || startRegion.end() + 1 != endRegion.start());
        // SGL SV. start + insert
        boolean sgl1 =
                startRegion != null && startOrientation != null && !insertSequence.isEmpty() && endRegion == null && endOrientation == null;
        // SGL SV. insert + end
        boolean sgl2 =
                startRegion == null && startOrientation == null && !insertSequence.isEmpty() && endRegion != null && endOrientation != null;
        if(!(single || genericVariant || sgl1 || sgl2))
        {
            throw new IllegalArgumentException("Invalid sequence definition");
        }
        if(startRegion != null && !startRegion.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid startRegion");
        }
        if(endRegion != null && !endRegion.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid endRegion");
        }
        if(!insertSequence.isEmpty() && !isDnaSequenceNormal(insertSequence))
        {
            throw new IllegalArgumentException("Invalid insertSequence: " + insertSequence);
        }
    }

    public static SequenceDefinition singleRegion(final ChrBaseRegion region)
    {
        return new SequenceDefinition(region, null, "", null, null);
    }

    public static SequenceDefinition simpleVariant(final ChrBaseRegion startRegion, final String insertSequence,
            final ChrBaseRegion endRegion)
    {
        return new SequenceDefinition(startRegion, Orientation.FORWARD, insertSequence, endRegion, Orientation.FORWARD);
    }

    public static SequenceDefinition structuralVariant(final ChrBaseRegion startRegion, final Orientation startOrientation,
            String insertSequence, final ChrBaseRegion endRegion, final Orientation endOrientation)
    {
        return new SequenceDefinition(startRegion, startOrientation, insertSequence, endRegion, endOrientation);
    }

    public static SequenceDefinition forwardSgl(final ChrBaseRegion startRegion, final String insertSequence)
    {
        return new SequenceDefinition(startRegion, Orientation.FORWARD, insertSequence, null, null);
    }

    public static SequenceDefinition reverseSgl(final String insertSequence, final ChrBaseRegion endRegion)
    {
        return new SequenceDefinition(null, null, insertSequence, endRegion, Orientation.FORWARD);
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
    public boolean isSingleRegion()
    {
        return startRegion != null && insertSequence.isEmpty() && endRegion == null;
    }

    // Gets the single region that the probe is defined by, or throws an exception.
    public ChrBaseRegion singleRegion()
    {
        ChrBaseRegion region = singleRegionOrNull();
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
    public ChrBaseRegion singleRegionOrNull()
    {
        return isSingleRegion() ? startRegion : null;
    }
}
