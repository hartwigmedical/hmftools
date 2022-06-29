package com.hartwig.hmftools.svtools.sv_prep;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.DEFAULT_FRAG_LENGTH_MAX;
import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.DEFAULT_FRAG_LENGTH_MIN;
import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.MIN_ALIGNMENT_BASES;
import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.MIN_INSERT_ALIGNMENT_OVERLAP;
import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.MIN_INSERT_LENGTH_SUPPORT;
import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.MIN_MAP_QUALITY;
import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.MIN_SOFT_CLIP_HIGH_QUAL_PERC;
import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.MIN_SOFT_CLIP_MIN_BASE_QUAL;
import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.MIN_SUPPORTING_READ_DISTANCE;

import static htsjdk.samtools.CigarOperator.M;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ReadFilterConfig
{
    /*
    - with least 1 read with 50 base alignment, MAPQ of 20 and abs(Insert size - M length) > 5 bases (ie test for short fragments with adapter) AND
	- a soft clip length of 30 with >=85% of soft clip bases with qual > 25.
	- at least 1 additional read which has any length soft clipping at the same base OR
	a within 50 bases with exact sequence match in between the soft clip locations.
     */

    public final int MinAlignmentBases;
    public final int MinMapQuality;
    public final int MinInsertAlignmentOverlap;
    public final int MinSoftClipLength;
    public final int MinSoftClipHighQual;
    public final double MinSoftClipHighQualPerc;
    public final int MinSupportingReadDistance;
    public final int MinInsertLengthSupport;

    private int mFragmentLengthMin;
    private int mFragmentLengthMax;

    public ReadFilterConfig()
    {
        this(MIN_ALIGNMENT_BASES, MIN_MAP_QUALITY, MIN_INSERT_ALIGNMENT_OVERLAP, MIN_SOFT_CLIP_LENGTH, MIN_SOFT_CLIP_MIN_BASE_QUAL,
                MIN_SOFT_CLIP_HIGH_QUAL_PERC, MIN_SUPPORTING_READ_DISTANCE);
    }

    public ReadFilterConfig(
            final int minAlignmentBases, final int minMapQuality, final int minInsertAlignmentOverlap, final int minSoftClipLength,
            final int minSoftClipHighQual, final double minSoftClipHighQualPerc, final int minSupportingReadDistance)
    {
        MinAlignmentBases = minAlignmentBases;
        MinMapQuality = minMapQuality;
        MinInsertAlignmentOverlap = minInsertAlignmentOverlap;
        MinSoftClipLength = minSoftClipLength;
        MinSoftClipHighQual = minSoftClipHighQual;
        MinSoftClipHighQualPerc = minSoftClipHighQualPerc;
        MinSupportingReadDistance = minSupportingReadDistance;

        MinInsertLengthSupport = MIN_INSERT_LENGTH_SUPPORT;

        mFragmentLengthMin = DEFAULT_FRAG_LENGTH_MIN;
        mFragmentLengthMax = DEFAULT_FRAG_LENGTH_MAX;
    }

    public int getFragmentLengthMin() { return mFragmentLengthMin; }
    public int getFragmentLengthMax() { return mFragmentLengthMax; }

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

        int insertAlignmentOverlap = abs(record.getInferredInsertSize() - alignedBases);

        if(insertAlignmentOverlap < MinInsertAlignmentOverlap)
            filters |= ReadFilterType.INSERT_MAP_OVERLAP.flag();

        int scLeft = cigar.isLeftClipped() ? cigar.getFirstCigarElement().getLength() : 0;
        int scRight = cigar.isRightClipped() ? cigar.getLastCigarElement().getLength() : 0;

        if(scLeft < MinSoftClipLength && scRight < MinSoftClipLength)
            filters |= ReadFilterType.SOFT_CLIP_LENGTH.flag();

        // base qual in soft clip
        final byte[] baseQualities = record.getBaseQualities();
        int scRangeStart = scLeft > scRight ? 0 : baseQualities.length - scLeft;
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

        return filters;
    }

    public boolean isCandidateSupportingRead(final SAMRecord record)
    {
        int insertLength = record.getCigar().getCigarElements().stream()
                .filter(x -> x.getOperator() == CigarOperator.I).mapToInt(x -> x.getLength()).findFirst().orElse(0);

        if(insertLength >= MinInsertLengthSupport)
            return true;

        if(record.getInferredInsertSize() > mFragmentLengthMax || record.getInferredInsertSize() < mFragmentLengthMin)
            return true;

        return record.getCigar().isLeftClipped() || record.getCigar().isRightClipped();
    }

}
