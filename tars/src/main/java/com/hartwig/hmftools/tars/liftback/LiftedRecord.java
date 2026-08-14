package com.hartwig.hmftools.tars.liftback;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

// One SAMRecord's candidate placements during lift-back: the read's own alignment plus its lifted XA alts, of which
// primaryIndex is the one written back.
public record LiftedRecord(
        int updatedMapQuality,
        // locus count as decided; deliberately NOT recomputed by later primary-only revisions
        int numLoci,
        String notes,
        // index of the chosen placement in liftedAlignments; 0 is the aligner's own primary, so anything else is a swap
        int primaryIndex,
        List<LiftedAlignment> liftedAlignments)
{
    public static final int NO_PRIMARY = -1;

    // Either the record has candidates and one of them is the chosen placement, or it has neither: a mismatch is a
    // lift bug, so fail here rather than downstream.
    public LiftedRecord
    {
        liftedAlignments = List.copyOf(liftedAlignments);
        boolean validIndex = liftedAlignments.isEmpty()
                ? primaryIndex == NO_PRIMARY
                : primaryIndex >= 0 && primaryIndex < liftedAlignments.size();

        if(!validIndex)
        {
            throw new IllegalArgumentException(String.format(
                    "primary index(%d) of %d alignment(s)", primaryIndex, liftedAlignments.size()));
        }
    }

    // No placement: the read arrived unmapped or lift-back gave up. The input record's own unmapped flag says which,
    // so only the reason is carried here.
    public static LiftedRecord unmapped(final String note)
    {
        return new LiftedRecord(0, 0, note, NO_PRIMARY, List.of());
    }

    // XA alts the aligner offered that lifted, dropped ones included: the alignment set is self plus those.
    public int numXaAlts()
    {
        return Math.max(liftedAlignments.size() - 1, 0);
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
        return primaryAlignment().LiftedChromosome;
    }

    public int finalPos()
    {
        return primaryAlignment().LiftedPos;
    }

    public String finalCigar()
    {
        return primaryAlignment().LiftedCigar;
    }

    public boolean negativeStrand()
    {
        return !primaryAlignment().ForwardStrand;
    }

    public boolean hasNCigar()
    {
        return primaryAlignment().cigarHasN();
    }

    // +1/-1 for a tx-contig-derived primary; 0 otherwise.
    public int transcriptStrand()
    {
        return primaryAlignment().TranscriptStrand;
    }

    public LiftedRecord withLiftedAlignments(final List<LiftedAlignment> alignments)
    {
        return new LiftedRecord(updatedMapQuality, numLoci, notes, primaryIndex, alignments);
    }

    public LiftedRecord withPrimaryTranscriptStrand(final int transcriptStrand)
    {
        List<LiftedAlignment> revised = new ArrayList<>(liftedAlignments);
        revised.set(primaryIndex, primaryAlignment().withTranscriptStrand(transcriptStrand));
        return withLiftedAlignments(revised);
    }

    // XA tag: every kept non-primary placement, each distinct locus once, minus alts overlapping the primary's span
    // (a shared-exon isoform read lifting back onto the primary's coords carries no alternative-position info).
    public String xaTag()
    {
        LiftedAlignment primary = primaryAlignment();
        Set<String> altEntries = new LinkedHashSet<>();

        for(int i = 0; i < liftedAlignments.size(); ++i)
        {
            LiftedAlignment alignment = liftedAlignments.get(i);
            if(i == primaryIndex || alignment.Dropped || alignment.overlaps(primary))
            {
                continue;
            }
            altEntries.add(alignment.toXaEntry());
        }

        return altEntries.isEmpty() ? null : String.join("", altEntries);
    }

}
