package com.hartwig.hmftools.bamtools.remapper;

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
