package com.hartwig.hmftools.redux.splice.rescue;

import java.util.Collections;
import java.util.List;

// Output DTO from JunctionRescueResolver. When merged() is true MergedCigar/MergedStart describe the
// new primary and DroppedSupplementaryIndices lists the supps that were absorbed. When merged() is
// false RejectReason carries the gate that filtered the candidate (NO_TERMINAL_SOFTCLIP if the
// primary never offered a softclip to extend across at all).
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
    public final RescueRejectReason RejectReason;   // null when Merged is true

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
