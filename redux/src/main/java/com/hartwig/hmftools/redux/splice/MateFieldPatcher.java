package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;

import htsjdk.samtools.SAMRecord;

// rewrites a record's mate-related fields (RNEXT / PNEXT / mate-strand / mate-unmapped / TLEN / proper-pair)
// from the partner read's lifted primary info held in the LiftedMateInfoCache.
//
// Applied in pass 2 of SpliceLiftBack after the record's own coordinates have been lifted, so that both
// sides of the pair refer to genomic (not transcript-contig) coordinates and supplementary alignments
// carry the same mate info as their primary.
public class MateFieldPatcher
{
    public static void patchMateFields(final SAMRecord record, final LiftedMateInfoCache liftedMateInfoCache)
    {
        if(!record.getReadPairedFlag())
            return;

        LiftedMateInfo partnerInfo = liftedMateInfoCache.getPartnerMateInfo(
                record.getReadName(), record.getFirstOfPairFlag());

        if(partnerInfo == null)
        {
            // partner primary not seen in pass 1 — clear proper-pair so downstream doesn't trust mate fields
            record.setProperPairFlag(false);
            return;
        }

        if(partnerInfo.Unmapped)
        {
            applyUnmappedPartner(record);
            return;
        }

        applyMappedPartner(record, partnerInfo);
    }

    // SAM convention when the mate is unmapped: place mate at the same position as the read itself,
    // mark mate-unmapped, drop proper-pair, zero out TLEN.
    private static void applyUnmappedPartner(final SAMRecord record)
    {
        record.setMateUnmappedFlag(true);
        record.setMateReferenceName(record.getReferenceName());
        record.setMateAlignmentStart(record.getAlignmentStart());
        record.setMateNegativeStrandFlag(false);
        record.setInferredInsertSize(0);
        record.setProperPairFlag(false);
        // MC is only meaningful for a mapped mate
        record.setAttribute(MATE_CIGAR_ATTRIBUTE, null);
    }

    private static void applyMappedPartner(final SAMRecord record, final LiftedMateInfo partnerInfo)
    {
        record.setMateUnmappedFlag(false);
        record.setMateReferenceName(partnerInfo.Chromosome);
        record.setMateAlignmentStart(partnerInfo.AlignmentStart);
        record.setMateNegativeStrandFlag(partnerInfo.NegativeStrand);
        record.setAttribute(MATE_CIGAR_ATTRIBUTE, partnerInfo.LiftedCigar);

        if(record.getReadUnmappedFlag())
        {
            // SAM convention: an unmapped read with a mapped mate is "placed" with its mate (same RNAME/POS).
            // Without this, the read keeps the pre-lift _tx contig name it inherited from its mate.
            record.setReferenceName(partnerInfo.Chromosome);
            record.setAlignmentStart(partnerInfo.AlignmentStart);
            record.setInferredInsertSize(0);
            record.setProperPairFlag(false);
            return;
        }

        if(!record.getReferenceName().equals(partnerInfo.Chromosome))
        {
            // TLEN is defined only when both reads are on the same reference and both mapped
            record.setInferredInsertSize(0);
            record.setProperPairFlag(false);
            return;
        }

        record.setInferredInsertSize(computeInferredInsertSize(record, partnerInfo));
    }

    // SAM spec TLEN: signed distance from the leftmost mapped base of the pair to the rightmost mapped base.
    // The leftmost read gets a positive value, the rightmost read gets the negation. When both start at the
    // same position, by convention the first-in-pair receives the positive value.
    static int computeInferredInsertSize(final SAMRecord record, final LiftedMateInfo partnerInfo)
    {
        final int readStart = record.getAlignmentStart();
        final int readEnd = record.getAlignmentEnd();
        final int mateStart = partnerInfo.AlignmentStart;
        final int mateEnd = partnerInfo.AlignmentEnd;

        final int leftmost = Math.min(readStart, mateStart);
        final int rightmost = Math.max(readEnd, mateEnd);
        final int span = rightmost - leftmost + 1;

        if(readStart < mateStart)
            return span;
        if(readStart > mateStart)
            return -span;

        // tie on start position — first-in-pair gets the positive value
        return record.getFirstOfPairFlag() ? span : -span;
    }
}
