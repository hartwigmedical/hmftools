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

    public HlaAlignmentPair closestAlignmentPair(RefGenomeVersion refGenomeVersion)
    {
        final HlaAlignmentPair principalsPair =
                new HlaAlignmentPair(new HlaAlignment(Lefts.principalAlignment()), new HlaAlignment(Rights.principalAlignment()));
        if(principalsPair.isConcordantPair())
        {
            return principalsPair;
        }
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
