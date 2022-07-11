package com.hartwig.hmftools.svprep.reads;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.SAMFlag.MATE_UNMAPPED;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

public class ReadFilters
{
    private final ReadFilterConfig mConfig;

    public ReadFilters(final ReadFilterConfig config)
    {
        mConfig = config;
    }

    public ReadFilterConfig config() { return mConfig; }

    public int checkFilters(final SAMRecord record)
    {
        int filters = 0;

        final Cigar cigar = record.getCigar();
        int alignedBases = cigar.getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

        if(alignedBases < mConfig.MinAlignmentBases)
            filters |= ReadFilterType.MIN_ALIGN_MATCH.flag();

        if(record.getMappingQuality() < mConfig.MinMapQuality)
            filters |= ReadFilterType.MIN_MAP_QUAL.flag();

        int insertAlignmentOverlap = abs(abs(record.getInferredInsertSize()) - alignedBases);

        if(insertAlignmentOverlap < mConfig.MinInsertAlignmentOverlap)
            filters |= ReadFilterType.INSERT_MAP_OVERLAP.flag();

        int maxIndelLength = ReadRecord.maxIndelLength(record.getCigar());

        if(maxIndelLength < mConfig.MinIndelLength)
        {
            int scLeft = cigar.isLeftClipped() ? cigar.getFirstCigarElement().getLength() : 0;
            int scRight = cigar.isRightClipped() ? cigar.getLastCigarElement().getLength() : 0;

            if(scLeft < mConfig.MinSoftClipLength && scRight < mConfig.MinSoftClipLength)
                filters |= ReadFilterType.SOFT_CLIP_LENGTH.flag();

            // base qual in soft clip
            if(scLeft > 0 || scRight > 0)
            {
                final byte[] baseQualities = record.getBaseQualities();
                int scRangeStart = scLeft > scRight ? 0 : baseQualities.length - scRight;
                int scRangeEnd = scLeft > scRight ? scLeft : baseQualities.length;
                double scLength = max(scLeft, scRight);

                int aboveQual = 0;
                for(int i = scRangeStart; i < scRangeEnd; ++i)
                {
                    if(baseQualities[i] >= mConfig.MinSoftClipHighQual)
                        ++aboveQual;
                }

                if(aboveQual / scLength < mConfig.MinSoftClipHighQualPerc)
                    filters |= ReadFilterType.SOFT_CLIP_BASE_QUAL.flag();
            }
        }

        return filters;
    }

    public boolean isCandidateSupportingRead(final SAMRecord record)
    {
        // any read with a supplementary
        if(record.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
            return true;

        // or an fragment length beyond the observed distribution
        if(abs(record.getInferredInsertSize()) > mConfig.fragmentLengthMax()) //  || record.getInferredInsertSize() < mFragmentLengthMin
            return true;

        // an unmapped mate
        if(record.getMateUnmappedFlag())
            return true;

        // interchromosomal
        if(!record.getContig().equals(record.getMateReferenceName()))
            return true;

        // inversion
        if(record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag())
            return true;

        // or with any amount of soft-clipping
        return record.getCigar().isLeftClipped() || record.getCigar().isRightClipped();
    }

}
