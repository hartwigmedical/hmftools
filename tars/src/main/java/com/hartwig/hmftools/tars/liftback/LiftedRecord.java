package com.hartwig.hmftools.tars.liftback;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

// Tars' in-flight view of one SAMRecord: 1:N over its candidate placements (own alignment + lifted XA alts), of which
// primaryIndex is the one written back; the placement accessors are valid only once hasPlacement() holds.
public record LiftedRecord(
        int updatedMapQuality,
        // locus count as decided; deliberately NOT recomputed, since a later soft-clip extension can change the alignment set
        int numLoci,
        String notes,
        // index of the chosen placement in liftedAlignments; 0 is the aligner's own primary, so anything else is a swap
        int primaryIndex,
        List<LiftedAlignment> liftedAlignments)
{
    public static final int NO_PRIMARY = -1;

    // Either the record has candidates and one of them is the chosen placement, or it has neither. The placement
    // accessors below assume it, so a mismatch is a lift bug and fails here rather than downstream.
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

    // No placement: either the read arrived unmapped, or lift-back gave up on it. Which one is answered by the input
    // record's own unmapped flag, so only the reason is carried here.
    public static LiftedRecord unmapped(final String note)
    {
        return new LiftedRecord(0, 0, note, NO_PRIMARY, List.of());
    }

    // alts the aligner offered in XA that lifted, dropped ones included: the set is self plus those.
    public int numXaAlts()
    {
        return Math.max(liftedAlignments.size() - 1, 0);
    }

    // the pick moved off the aligner's own primary at element 0
    public boolean swapped()
    {
        return primaryIndex > 0;
    }

    // false when nothing lifted; the accessors below throw in that case.
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

    // Revises the chosen primary's placement (overhang collapse, supplementary merge) with a new MAPQ and an appended note.
    // The revision lands on the alignment, so the placement and the alignment set cannot drift apart.
    public LiftedRecord withRevisedPrimary(
            final int newPos, final String newCigar, final int newUpdatedMapQuality, final String note)
    {
        List<LiftedAlignment> revised = new ArrayList<>(liftedAlignments);
        revised.set(primaryIndex, primaryAlignment().withLiftedCigar(newPos, newCigar));

        return new LiftedRecord(
                newUpdatedMapQuality, numLoci, appendNote(notes, note), primaryIndex, revised);
    }

    // The XA tag for this record: every kept non-primary placement, minus alts whose span overlaps the primary's (a
    // shared-exon isoform read lifting back onto the primary's coords carries no alternative-position info), each
    // distinct locus once. Null when nothing is left to report.
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

    private static String appendNote(final String existing, final String note)
    {
        if(existing == null || existing.isEmpty())
        {
            return note;
        }
        return existing + ";" + note;
    }
}
