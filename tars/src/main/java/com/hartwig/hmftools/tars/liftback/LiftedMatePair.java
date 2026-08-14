package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_INDEX;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAPP_COORDS_DELIM;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.firstInPair;

import htsjdk.samtools.SAMRecord;

// One read pair's lifted records, and the mate-field patch they serve: a record's RNEXT/PNEXT/mate-strand/MC/TLEN come
// from its mate's lifted placement, known only once that mate has been lifted. Scoped to a single name group, so the
// pair is two slots rather than a keyed cache.
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

    // null until the mate has been lifted: /1 is decided before /2 exists, so it sees no mate.
    public LiftedRecord mateOf(final boolean firstOfPair)
    {
        return firstOfPair ? mSecondInPair : mFirstInPair;
    }

    // the caller's own side, used to mirror primary coords onto a supplementary whose own lift failed
    public LiftedRecord ownPrimary(final boolean firstOfPair)
    {
        return firstOfPair ? mFirstInPair : mSecondInPair;
    }

    // TARS already has both final primary decisions for the name group, so write the final pair state directly: park an
    // unmapped read at its mapped mate, or clear both sides if neither maps.
    public void unmapRead(final SAMRecord record)
    {
        LiftedRecord mate = record.getReadPairedFlag() ? mateOf(firstInPair(record)) : null;
        String originalChromosome = record.getReferenceName();
        int originalPosition = record.getAlignmentStart();

        record.setReadUnmappedFlag(true);
        record.setMappingQuality(0);
        record.setProperPairFlag(false);
        record.setInferredInsertSize(0);
        if(!record.isSecondaryOrSupplementary())
        {
            String chromosome = originalChromosome.replaceAll(UNMAPP_COORDS_DELIM, "");
            record.setAttribute(UNMAP_ATTRIBUTE, chromosome + UNMAPP_COORDS_DELIM + originalPosition);
        }

        if(mate != null && mate.hasPlacement())
        {
            record.setReferenceName(mate.finalChromosome());
            record.setAlignmentStart(mate.finalPos());
            record.setMateUnmappedFlag(false);
            record.setMateReferenceName(mate.finalChromosome());
            record.setMateAlignmentStart(mate.finalPos());
            record.setMateNegativeStrandFlag(mate.negativeStrand());
            record.setAttribute(MATE_CIGAR_ATTRIBUTE, mate.finalCigar());
        }
        else
        {
            clearReadCoordinates(record);
            if(record.getReadPairedFlag())
            {
                record.setMateUnmappedFlag(true);
                clearMateCoordinates(record);
                record.setAttribute(MATE_CIGAR_ATTRIBUTE, null);
            }
        }

        record.setCigarString(NO_CIGAR);
        if(record.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
        {
            record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, null);
        }
    }

    // sets RNEXT/PNEXT/mate-strand/MC and TLEN from the mate's lifted primary
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
        record.setMateUnmappedFlag(true);
        record.setProperPairFlag(false);
        record.setInferredInsertSize(0);
        record.setAttribute(MATE_CIGAR_ATTRIBUTE, null);

        if(record.getReadUnmappedFlag())
        {
            clearReadCoordinates(record);
            clearMateCoordinates(record);
        }
        else
        {
            record.setMateReferenceName(record.getReferenceName());
            record.setMateAlignmentStart(record.getAlignmentStart());
        }
    }

    private static void clearReadCoordinates(final SAMRecord record)
    {
        record.setAlignmentStart(NO_POSITION);
        record.setReferenceIndex(NO_CHROMOSOME_INDEX);
        record.setReferenceName(NO_CHROMOSOME_NAME);
    }

    private static void clearMateCoordinates(final SAMRecord record)
    {
        record.setMateAlignmentStart(NO_POSITION);
        record.setMateReferenceIndex(NO_CHROMOSOME_INDEX);
        record.setMateReferenceName(NO_CHROMOSOME_NAME);
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

    // TLEN: signed distance between mates' 5' ends, +/-1 so the leftmost-5' mate is positive. Using 5' ends rather than
    // alignment start fixes the sign for same-start pairs and the magnitude when softclips extend an end.
    static int computeInferredInsertSize(final SAMRecord record, final LiftedRecord mate)
    {
        int readFivePrime = record.getReadNegativeStrandFlag() ? record.getAlignmentEnd() : record.getAlignmentStart();
        int mateFivePrime = mate.negativeStrand() ? mate.primaryAlignment().alignedEnd() : mate.finalPos();

        int adjustment = mateFivePrime >= readFivePrime ? 1 : -1;
        return mateFivePrime - readFivePrime + adjustment;
    }
}
