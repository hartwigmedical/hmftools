package com.hartwig.hmftools.tars.liftback.supplementary;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Output DTO from SupplementaryResolver. On success, mergedCigar/mergedStart describe the new primary
// and droppedSupplementaryIndices lists absorbed supps. On failure, rejectReason carries the gate hit.
public record SupplementaryResult(
        boolean merged, String mergedCigar, int mergedStart, List<Integer> droppedSupplementaryIndices,
        List<ChrBaseRegion> introducedIntrons, int chainDepth, SupplementaryRejectReason rejectReason)
{
    private static final SupplementaryResult NO_MERGE_NO_OP = new SupplementaryResult(
            false, null, -1, Collections.emptyList(), Collections.emptyList(), 0,
            SupplementaryRejectReason.NO_TERMINAL_SOFTCLIP);

    public static SupplementaryResult noMerge(final SupplementaryRejectReason reason)
    {
        if(reason == SupplementaryRejectReason.NO_TERMINAL_SOFTCLIP)
        {
            return NO_MERGE_NO_OP;
        }
        return new SupplementaryResult(false, null, -1, Collections.emptyList(), Collections.emptyList(), 0, reason);
    }
}
