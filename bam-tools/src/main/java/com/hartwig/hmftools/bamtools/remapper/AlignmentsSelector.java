package com.hartwig.hmftools.bamtools.remapper;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;

public class AlignmentsSelector
{
    @NotNull
    private final AlignmentsList lefts;
    @NotNull
    private final AlignmentsList rights;
    private @NotNull
    final RefGenomeVersion refGenomeVersion = RefGenomeVersion.V38; // TODO make parameter

    public AlignmentsSelector(@NotNull final AlignmentsList lefts, @NotNull final AlignmentsList rights)
    {
        this.lefts = lefts;
        this.rights = rights;
    }

    /*
    Alg for best pair
    AlignmentSelector can return Pair<BwaMemAlignment, BwaMemAlignment> consisting of the principal left and right reads
    New class HlaPositions takes a BwaMemAlignment and produces a list of all positions within HLA regions - maybe call these HlaAlignment
    and give them a reference to their originating BwaMemAlignment.
    AlignmentSelector can thus produce a set of paired HlaAlignment objects. Call these objects AlignmentPair
    AlignmentPair has a distance function - the distance between the positions.
    Choose the pair having the least distance.

     */

    public HlaAlignmentPair closestAlignmentPair(RefGenomeVersion refGenomeVersion)
    {
//        return new HlaAlignmentPair(new HlaAlignment(lefts.principalAlignment()), new HlaAlignment(rights.principalAlignment()));

        List<HlaAlignmentPair> allPairs = new ArrayList<>();
        HlaAlignment.hlaAlignments(lefts.principalAlignment(), refGenomeVersion).forEach(left ->
        {
            HlaAlignment.hlaAlignments(rights.principalAlignment(), refGenomeVersion)
                    .forEach(right -> allPairs.add(new HlaAlignmentPair(left, right)));
        });
//        allPairs.forEach(pair -> {
//            System.out.println(" " + pair.interPairDistance());
//        });
        return allPairs.stream().sorted().findFirst().orElseThrow();
    }

    public boolean hasMultipleLeftAlignments()
    {
        return lefts.hasMultipleAlignments();
    }

    public boolean hasMultipleRightAlignments()
    {
        return rights.hasMultipleAlignments();
    }

    public boolean hasMultipleLeftAndRightAlignments()
    {
        return hasMultipleLeftAlignments() && hasMultipleRightAlignments();
    }

    public List<PreferredAlignment> preferredLeftAlignments()
    {
        BwaMemAlignment mate = rights.alignments.get(0);
        return lefts.alignments.stream().map(a -> new PreferredAlignment(a, mate, refGenomeVersion)).collect(Collectors.toList());
    }

    public List<PreferredAlignment> preferredRightAlignments()
    {
        BwaMemAlignment mate = lefts.alignments.get(0);
        return rights.alignments.stream().map(a -> new PreferredAlignment(a, mate, refGenomeVersion)).collect(Collectors.toList());
    }
}
