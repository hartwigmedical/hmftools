package com.hartwig.hmftools.panelbuilder.samplevariants;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

// Possibilities:
//   - The probe is all ref sequence:
//       `start` is non-null, others are null.
//       probe = start
//   - The probe has an alt sequence or SV:
//       `sequence` is non-null, at least one of `start` or `end` are non-null.
//       probe = start + insert + end
public record VariantProbeData(
        @Nullable String sequence,
        @Nullable ChrBaseRegion start,
        @Nullable String insert,
        @Nullable ChrBaseRegion end
)
{
    public VariantProbeData
    {
        boolean valid1 = sequence == null && start != null && insert == null && end == null;
        boolean valid2 = sequence != null && (start != null || end != null);
        if(!(valid1 || valid2))
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

    public boolean hasAltSequence()
    {
        return insert != null && !insert.isEmpty();
    }
}
