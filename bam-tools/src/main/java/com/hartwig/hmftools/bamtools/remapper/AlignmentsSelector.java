package com.hartwig.hmftools.bamtools.remapper;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class AlignmentsSelector
{
    private final AlignmentsList mLefts;
    private final AlignmentsList mRights;

    public AlignmentsSelector(final AlignmentsList lefts, final AlignmentsList rights)
    {
        mLefts = lefts;
        mRights = rights;
    }

    public HlaAlignmentPair bestAlignmentPair(RefGenomeVersion refGenomeVersion)
    {
        final HlaAlignment principalLeft = new HlaAlignment(mLefts.principalAlignment());
        final HlaAlignment principalRight = new HlaAlignment(mRights.principalAlignment());
        final HlaAlignmentPair principalsPair = new HlaAlignmentPair(principalLeft, principalRight);
        if(principalsPair.isConcordantPair() || principalsPair.oneOfPairIsUnmapped())
        {
            return principalsPair;
        }
        // If the principal pair is not an obvious choice then we choose a pair involving one or more supplementary alignments.
        // To make this deterministic we choose the closest pair, but there's no reason this is actually best.
        List<HlaAlignmentPair> allPairs = new ArrayList<>();
        HlaAlignment.hlaAlignments(mLefts.principalAlignment(), refGenomeVersion).forEach(left ->
                HlaAlignment.hlaAlignments(mRights.principalAlignment(), refGenomeVersion)
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
