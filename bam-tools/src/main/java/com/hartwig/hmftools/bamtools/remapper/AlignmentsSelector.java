package com.hartwig.hmftools.bamtools.remapper;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.jetbrains.annotations.NotNull;

public class AlignmentsSelector
{
    @NotNull
    private final AlignmentsList Lefts;
    @NotNull
    private final AlignmentsList Rights;

    public AlignmentsSelector(@NotNull final AlignmentsList lefts, @NotNull final AlignmentsList rights)
    {
        this.Lefts = lefts;
        this.Rights = rights;
    }

    public HlaAlignmentPair bestAlignmentPair(RefGenomeVersion refGenomeVersion)
    {
        final HlaAlignment principalLeft = new HlaAlignment(Lefts.principalAlignment());
        final HlaAlignment principalRight = new HlaAlignment(Rights.principalAlignment());
        final HlaAlignmentPair principalsPair = new HlaAlignmentPair(principalLeft, principalRight);
        if(principalsPair.isConcordantPair() || principalsPair.oneOfPairIsUnmapped())
        {
            return principalsPair;
        }
        // If the principal pairis not an obvious choice then we
        // choose a pair involving one or more supplementary alignments.
        // To make this deterministic we choose the closest pair,
        // but there's no reason this is actually best.
        List<HlaAlignmentPair> allPairs = new ArrayList<>();
        HlaAlignment.hlaAlignments(Lefts.principalAlignment(), refGenomeVersion).forEach(left ->
                HlaAlignment.hlaAlignments(Rights.principalAlignment(), refGenomeVersion)
                        .forEach(right -> allPairs.add(new HlaAlignmentPair(left, right))));
        HlaAlignmentPair bestPair = allPairs.stream()
                .filter(hlaAlignmentPair -> hlaAlignmentPair.interPairDistance() <= 1000).
                sorted().findFirst().orElse(null);
        if(bestPair == null)
        {
            return principalsPair;
        }
        return bestPair;
    }
}
