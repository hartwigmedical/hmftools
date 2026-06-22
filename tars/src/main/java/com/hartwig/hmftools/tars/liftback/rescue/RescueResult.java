package com.hartwig.hmftools.tars.liftback.rescue;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Output DTO from JunctionRescueResolver. On success, mergedCigar/mergedStart describe the new primary
// and droppedSupplementaryIndices lists absorbed supps. On failure, rejectReason carries the gate hit.
public record RescueResult(
        boolean merged, String mergedCigar, int mergedStart, List<Integer> droppedSupplementaryIndices,
        List<ChrBaseRegion> introducedIntrons, int chainDepth, RescueRejectReason rejectReason)
{
    private static final RescueResult NO_MERGE_NO_OP = new RescueResult(
            false, null, -1, Collections.emptyList(), Collections.emptyList(), 0,
            RescueRejectReason.NO_TERMINAL_SOFTCLIP);

    public static RescueResult noMerge(final RescueRejectReason reason)
    {
        if(reason == RescueRejectReason.NO_TERMINAL_SOFTCLIP)
        {
            return NO_MERGE_NO_OP;
        }
        return new RescueResult(false, null, -1, Collections.emptyList(), Collections.emptyList(), 0, reason);
    }
}
