package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.firstInPair;

import htsjdk.samtools.SAMRecord;

// One read pair's lifted records, and the mate-field patch they exist to serve: a record's RNEXT/PNEXT/mate-strand/
// MC/TLEN come from its mate's lifted placement, which is only known once that mate has been lifted. Scoped to a
// single name-group, so the pair is two slots rather than a keyed cache.
public class LiftedMatePair
{
    private LiftedRecord mFirstInPair;
    private LiftedRecord mSecondInPair;

    public void recordPrimary(final boolean firstOfPair, final LiftedRecord liftedRecord)
    {
        if(firstOfPair)
        {
            mFirstInPair = liftedRecord;
        }
        else
        {
            mSecondInPair = liftedRecord;
        }
    }

    // null until the mate has been lifted - /1 is decided before /2 exists, so it sees no mate.
    public LiftedRecord mateOf(final boolean firstOfPair)
    {
        return firstOfPair ? mSecondInPair : mFirstInPair;
    }

    // the caller's own side, used to mirror primary coords onto a supplementary whose own lift failed
    public LiftedRecord ownPrimary(final boolean firstOfPair)
    {
        return firstOfPair ? mFirstInPair : mSecondInPair;
    }

    // Patches a record's mate fields (RNEXT/PNEXT/mate-strand/MC) and TLEN from its mate's lifted primary.
    public void patchMateFields(final SAMRecord record)
    {
        if(!record.getReadPairedFlag())
        {
            return;
        }

        LiftedRecord mate = mateOf(firstInPair(record));

        if(mate == null || !mate.hasPlacement())
        {
            applyUnmappedMate(record);
            return;
        }

        applyMappedMate(record, mate);
    }

    private static void applyUnmappedMate(final SAMRecord record)
    {
        if(record.getReadUnmappedFlag())
        {
            record.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
            record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        }
        record.setMateUnmappedFlag(true);
        record.setMateReferenceName(record.getReferenceName());
        record.setMateAlignmentStart(record.getAlignmentStart());
        record.setMateNegativeStrandFlag(false);
        record.setInferredInsertSize(0);
        record.setProperPairFlag(false);
        record.setAttribute(MATE_CIGAR_ATTRIBUTE, null);
    }

    private static void applyMappedMate(final SAMRecord record, final LiftedRecord mate)
    {
        record.setMateUnmappedFlag(false);
        record.setMateReferenceName(mate.finalChromosome());
        record.setMateAlignmentStart(mate.finalPos());
        record.setMateNegativeStrandFlag(mate.negativeStrand());
        record.setAttribute(MATE_CIGAR_ATTRIBUTE, mate.finalCigar());

        if(record.getReadUnmappedFlag())
        {
            record.setReferenceName(mate.finalChromosome());
            record.setAlignmentStart(mate.finalPos());
            record.setInferredInsertSize(0);
            record.setProperPairFlag(false);
            return;
        }

        if(!record.getReferenceName().equals(mate.finalChromosome()))
        {
            record.setInferredInsertSize(0);
            record.setProperPairFlag(false);
            return;
        }

        record.setInferredInsertSize(computeInferredInsertSize(record, mate));
    }

    // TLEN: signed distance between mates' 5' ends (strand-aware), +/-1 so the leftmost-5' mate is positive.
    // Using 5' ends (not alignment-start) fixes sign for same-start pairs and magnitude when softclips extend an end.
    static int computeInferredInsertSize(final SAMRecord record, final LiftedRecord mate)
    {
        int readFivePrime = record.getReadNegativeStrandFlag() ? record.getAlignmentEnd() : record.getAlignmentStart();
        int mateFivePrime = mate.negativeStrand() ? mate.primaryAlignment().alignedEnd() : mate.finalPos();

        int adjustment = mateFivePrime >= readFivePrime ? 1 : -1;
        return mateFivePrime - readFivePrime + adjustment;
    }
}
