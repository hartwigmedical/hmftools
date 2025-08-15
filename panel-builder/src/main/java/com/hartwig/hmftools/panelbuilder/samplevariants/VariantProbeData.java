package com.hartwig.hmftools.panelbuilder.samplevariants;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

// Data that defines a probe creates for a variant.
// Exists to decouple the code which determines the probe sequence from the rest of the probe generation code.
public record VariantProbeData(
        // Full probe sequence.
        String sequence,
        // probe = start + insert + end
        @Nullable ChrBaseRegion start,
        @Nullable String insert,
        @Nullable ChrBaseRegion end
)
{
    public VariantProbeData
    {
        if(!(start != null || end != null))
        {
            throw new IllegalArgumentException();
        }
        if(start != null && !start.isValid())
        {
            throw new IllegalArgumentException();
        }
        if(end != null && !end.isValid())
        {
            throw new IllegalArgumentException();
        }
        int startLength = start == null ? 0 : start.baseLength();
        int insertLength = insert == null ? 0 : insert.length();
        int endLength = end == null ? 0 : end.baseLength();
        if(sequence.length() != startLength + insertLength + endLength)
        {
            throw new IllegalArgumentException();
        }
    }

    public List<ChrBaseRegion> regions()
    {
        List<ChrBaseRegion> result = new ArrayList<>();
        if(start != null)
        {
            result.add(start);
        }
        if(end != null)
        {
            result.add(end);
        }
        return result;
    }
}
