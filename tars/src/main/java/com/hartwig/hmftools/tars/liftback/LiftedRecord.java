package com.hartwig.hmftools.tars.liftback;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.genome.region.Orientation;

import org.jetbrains.annotations.Nullable;

public record LiftedRecord(
        int updatedMapQuality,
        // locus count as decided; deliberately NOT recomputed by later primary-only revisions
        int numLoci,
        // index of the chosen placement in liftedAlignments; 0 is the aligner's own primary, so anything else is a swap
        int primaryIndex,
        List<LiftedAlignment> liftedAlignments)
{
    public static final int NO_PRIMARY = -1;

    public LiftedRecord
    {
        boolean validIndex = liftedAlignments.isEmpty()
                ? primaryIndex == NO_PRIMARY
                : primaryIndex >= 0 && primaryIndex < liftedAlignments.size();

        if(!validIndex)
        {
            throw new IllegalArgumentException(String.format(
                    "primary index(%d) of %d alignment(s)", primaryIndex, liftedAlignments.size()));
        }
    }

    public static LiftedRecord unmapped()
    {
        return new LiftedRecord(0, 0, NO_PRIMARY, List.of());
    }

    // dropped alts included
    public int numXaAlts()
    {
        int altCount = 0;
        for(int i = 1; i < liftedAlignments.size(); ++i)
        {
            if(!liftedAlignments.get(i).hasSupplementaryMerge()) ++altCount;
        }

        return altCount;
    }

    // false when nothing lifted; the placement accessors below throw in that case.
    public boolean hasPlacement()
    {
        return primaryIndex >= 0;
    }

    public LiftedAlignment primaryAlignment()
    {
        return liftedAlignments.get(primaryIndex);
    }

    public String finalChromosome()
    {
        return primaryAlignment().liftedPosition().Chromosome;
    }

    public int finalPos()
    {
        return primaryAlignment().liftedPosition().Position;
    }

    public String finalCigar()
    {
        return primaryAlignment().liftedCigar();
    }

    public boolean negativeStrand()
    {
        return !primaryAlignment().forwardStrand();
    }

    public boolean hasNCigar()
    {
        return primaryAlignment().cigarHasN();
    }

    // +1/-1 for a tx-contig-derived primary; 0 otherwise.
    @Nullable
    public Orientation transcriptStrand()
    {
        return primaryAlignment().transcriptStrand();
    }

    public LiftedRecord withLiftedAlignments(final List<LiftedAlignment> alignments)
    {
        return new LiftedRecord(updatedMapQuality, numLoci, primaryIndex, alignments);
    }

    public LiftedRecord withPrimaryIndex(final int newPrimaryIndex)
    {
        return new LiftedRecord(updatedMapQuality, numLoci, newPrimaryIndex, liftedAlignments);
    }

    public LiftedRecord withMapQuality(final int mapQuality)
    {
        return new LiftedRecord(mapQuality, numLoci, primaryIndex, liftedAlignments);
    }

    public LiftedRecord withPrimaryAlignment(final LiftedAlignment alignment)
    {
        List<LiftedAlignment> revised = new ArrayList<>(liftedAlignments);
        revised.set(primaryIndex, alignment);
        return withLiftedAlignments(revised);
    }

    // alts overlapping the primary's span are excluded: a shared-exon isoform read lifting back onto the primary's
    // coords carries no alternative-position info
    public String xaTag()
    {
        LiftedAlignment primary = primaryAlignment();
        Set<String> altEntries = new LinkedHashSet<>();

        for(int i = 0; i < liftedAlignments.size(); ++i)
        {
            LiftedAlignment alignment = liftedAlignments.get(i);
            if(i == primaryIndex || alignment.dropped() || alignment.overlaps(primary))
            {
                continue;
            }
            altEntries.add(alignment.toXaEntry());
        }

        return altEntries.isEmpty() ? null : String.join("", altEntries);
    }

}
