package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;

import htsjdk.samtools.SAMRecord;

public class MateFieldPatcher
{
    public static void patchMateFields(final SAMRecord record, final LiftedMateInfoCache liftedMateInfoCache)
    {
        if(!record.getReadPairedFlag())
            return;

        LiftedMateInfo partnerInfo = liftedMateInfoCache.getPartnerMateInfo(
                record.getReadName(), record.getFirstOfPairFlag());

        if(partnerInfo == null || partnerInfo.unmapped())
        {
            applyUnmappedPartner(record);
            return;
        }

        applyMappedPartner(record, partnerInfo);
    }

    private static void applyUnmappedPartner(final SAMRecord record)
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

    private static void applyMappedPartner(final SAMRecord record, final LiftedMateInfo partnerInfo)
    {
        record.setMateUnmappedFlag(false);
        record.setMateReferenceName(partnerInfo.chromosome());
        record.setMateAlignmentStart(partnerInfo.alignmentStart());
        record.setMateNegativeStrandFlag(partnerInfo.negativeStrand());
        record.setAttribute(MATE_CIGAR_ATTRIBUTE, partnerInfo.liftedCigar());

        if(record.getReadUnmappedFlag())
        {
            record.setReferenceName(partnerInfo.chromosome());
            record.setAlignmentStart(partnerInfo.alignmentStart());
            record.setInferredInsertSize(0);
            record.setProperPairFlag(false);
            return;
        }

        if(!record.getReferenceName().equals(partnerInfo.chromosome()))
        {
            record.setInferredInsertSize(0);
            record.setProperPairFlag(false);
            return;
        }

        record.setInferredInsertSize(computeInferredInsertSize(record, partnerInfo));
    }

    // SAM spec TLEN: signed pair span. Leftmost read gets +span, rightmost gets -span. On a start-position
    // tie, first-of-pair gets the positive value.
    static int computeInferredInsertSize(final SAMRecord record, final LiftedMateInfo partnerInfo)
    {
        final int readStart = record.getAlignmentStart();
        final int readEnd = record.getAlignmentEnd();
        final int mateStart = partnerInfo.alignmentStart();
        final int mateEnd = partnerInfo.alignmentEnd();

        final int leftmost = Math.min(readStart, mateStart);
        final int rightmost = Math.max(readEnd, mateEnd);
        final int span = rightmost - leftmost + 1;

        if(readStart < mateStart)
            return span;
        if(readStart > mateStart)
            return -span;

        return record.getFirstOfPairFlag() ? span : -span;
    }
}
