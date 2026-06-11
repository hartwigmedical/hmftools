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

    // TLEN: signed distance between mates' 5' ends (strand-aware), +/-1 so the leftmost-5' mate is positive.
    // Using 5' ends (not alignment-start) fixes sign for same-start pairs and magnitude when softclips extend an end.
    static int computeInferredInsertSize(final SAMRecord record, final LiftedMateInfo partnerInfo)
    {
        final int readFivePrime = record.getReadNegativeStrandFlag() ? record.getAlignmentEnd() : record.getAlignmentStart();
        final int mateFivePrime = partnerInfo.negativeStrand() ? partnerInfo.alignmentEnd() : partnerInfo.alignmentStart();

        final int adjustment = mateFivePrime >= readFivePrime ? 1 : -1;
        return mateFivePrime - readFivePrime + adjustment;
    }
}
