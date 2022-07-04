package com.hartwig.hmftools.svprep;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.svprep.SvConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.MAX_FRAGMENT_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_ALIGNMENT_BASES;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_INSERT_ALIGNMENT_OVERLAP;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_INSERT_LENGTH_SUPPORT;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_MAP_QUALITY;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_SOFT_CLIP_HIGH_QUAL_PERC;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_SOFT_CLIP_MIN_BASE_QUAL;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_SUPPORTING_READ_DISTANCE;

import static htsjdk.samtools.CigarOperator.M;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ReadFilterConfig
{
    public final int MinAlignmentBases;
    public final int MinMapQuality;
    public final int MinInsertAlignmentOverlap;
    public final int MinSoftClipLength;
    public final int MinSoftClipHighQual;
    public final double MinSoftClipHighQualPerc;
    public final int MinSupportingReadDistance;
    public final int MinInsertLengthSupport;

    public final int MinJunctionSupport;

    private int mFragmentLengthMin;
    private int mFragmentLengthMax;

    public ReadFilterConfig()
    {
        this(MIN_ALIGNMENT_BASES, MIN_MAP_QUALITY, MIN_INSERT_ALIGNMENT_OVERLAP, MIN_SOFT_CLIP_LENGTH, MIN_SOFT_CLIP_MIN_BASE_QUAL,
                MIN_SOFT_CLIP_HIGH_QUAL_PERC, MIN_SUPPORTING_READ_DISTANCE, MIN_JUNCTION_SUPPORT);
    }

    public ReadFilterConfig(
            final int minAlignmentBases, final int minMapQuality, final int minInsertAlignmentOverlap, final int minSoftClipLength,
            final int minSoftClipHighQual, final double minSoftClipHighQualPerc, final int minSupportingReadDistance,
            final int minJunctionSupport)
    {
        MinAlignmentBases = minAlignmentBases;
        MinMapQuality = minMapQuality;
        MinInsertAlignmentOverlap = minInsertAlignmentOverlap;
        MinSoftClipLength = minSoftClipLength;
        MinSoftClipHighQual = minSoftClipHighQual;
        MinSoftClipHighQualPerc = minSoftClipHighQualPerc;
        MinSupportingReadDistance = minSupportingReadDistance;

        MinInsertLengthSupport = MIN_INSERT_LENGTH_SUPPORT;
        MinJunctionSupport = minJunctionSupport;

        mFragmentLengthMin = DEFAULT_READ_LENGTH;
        mFragmentLengthMax = MAX_FRAGMENT_LENGTH;
    }

    public void setFragmentLengthMin(int minLength, int maxLength)
    {
        mFragmentLengthMin = minLength;
        mFragmentLengthMax = maxLength;
    }

    public int checkFilters(final SAMRecord record)
    {
        int filters = 0;

        final Cigar cigar = record.getCigar();
        int alignedBases = cigar.getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

        if(alignedBases < MinAlignmentBases)
            filters |= ReadFilterType.MIN_ALIGN_MATCH.flag();

        if(record.getMappingQuality() < MinMapQuality)
            filters |= ReadFilterType.MIN_MAP_QUAL.flag();

        int insertAlignmentOverlap = abs(abs(record.getInferredInsertSize()) - alignedBases);

        if(insertAlignmentOverlap < MinInsertAlignmentOverlap)
            filters |= ReadFilterType.INSERT_MAP_OVERLAP.flag();

        int scLeft = cigar.isLeftClipped() ? cigar.getFirstCigarElement().getLength() : 0;
        int scRight = cigar.isRightClipped() ? cigar.getLastCigarElement().getLength() : 0;

        if(scLeft < MinSoftClipLength && scRight < MinSoftClipLength)
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
                if(baseQualities[i] >= MinSoftClipHighQual)
                    ++aboveQual;
            }

            if(aboveQual / scLength < MinSoftClipHighQualPerc)
                filters |= ReadFilterType.SOFT_CLIP_BASE_QUAL.flag();
        }

        return filters;
    }

    public boolean isCandidateSupportingRead(final SAMRecord record)
    {
        int insertLength = record.getCigar().getCigarElements().stream()
                .filter(x -> x.getOperator() == CigarOperator.I).mapToInt(x -> x.getLength()).findFirst().orElse(0);

        if(insertLength >= MinInsertLengthSupport)
            return true;

        /* this is to identify discordant fragments to support known junctions
        if(record.getInferredInsertSize() > mFragmentLengthMax || record.getInferredInsertSize() < mFragmentLengthMin)
            return true;
        */

        return record.getCigar().isLeftClipped() || record.getCigar().isRightClipped();
    }

}
