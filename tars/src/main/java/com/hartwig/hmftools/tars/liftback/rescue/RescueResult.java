package com.hartwig.hmftools.tars.liftback.rescue;

import java.util.Collections;
import java.util.List;

// Output DTO from JunctionRescueResolver. On success, MergedCigar/MergedStart describe the new primary
// and DroppedSupplementaryIndices lists absorbed supps. On failure, RejectReason carries the gate hit.
public class RescueResult
{
    private static final RescueResult NO_MERGE_NO_OP = new RescueResult(
            false, null, -1, Collections.emptyList(), Collections.emptyList(), 0,
            RescueRejectReason.NO_TERMINAL_SOFTCLIP);

    public final boolean Merged;
    public final String MergedCigar;
    public final int MergedStart;
    public final List<Integer> DroppedSupplementaryIndices;
    public final List<ChrIntron> IntroducedIntrons;
    public final int ChainDepth;
    public final RescueRejectReason RejectReason;   // null on success

    public RescueResult(
            final boolean merged, final String mergedCigar, final int mergedStart,
            final List<Integer> droppedSupplementaryIndices, final List<ChrIntron> introducedIntrons,
            final int chainDepth, final RescueRejectReason rejectReason)
    {
        Merged = merged;
        MergedCigar = mergedCigar;
        MergedStart = mergedStart;
        DroppedSupplementaryIndices = droppedSupplementaryIndices;
        IntroducedIntrons = introducedIntrons;
        ChainDepth = chainDepth;
        RejectReason = rejectReason;
    }

    public static RescueResult noMerge(final RescueRejectReason reason)
    {
        if(reason == RescueRejectReason.NO_TERMINAL_SOFTCLIP)
            return NO_MERGE_NO_OP;
        return new RescueResult(false, null, -1, Collections.emptyList(), Collections.emptyList(), 0, reason);
    }

    public boolean merged()
    {
        return Merged;
    }
}
