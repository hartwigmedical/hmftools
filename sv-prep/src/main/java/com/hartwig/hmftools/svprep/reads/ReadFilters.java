package com.hartwig.hmftools.svprep.reads;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.samtools.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.samtools.CigarUtils.leftSoftClipped;
import static com.hartwig.hmftools.common.samtools.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.samtools.CigarUtils.rightSoftClipped;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateUnmapped;
import static com.hartwig.hmftools.common.region.ExcludedRegions.POLY_C_INSERT;
import static com.hartwig.hmftools.common.region.ExcludedRegions.POLY_G_INSERT;
import static com.hartwig.hmftools.common.region.ExcludedRegions.POLY_G_LENGTH;
import static com.hartwig.hmftools.common.sv.LineElements.isMobileLineElement;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.SvConstants.MAX_SOFT_CLIP_LOW_QUAL_COUNT;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_INDEL_SUPPORT_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_LINE_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.REPEAT_BREAK_CHECK_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.REPEAT_BREAK_MIN_MAP_QUAL;
import static com.hartwig.hmftools.svprep.SvConstants.REPEAT_BREAK_MIN_SC_LENGTH;
import static com.hartwig.hmftools.svprep.reads.ReadFilterType.POLY_G_SC;
import static com.hartwig.hmftools.svprep.reads.ReadFilterType.SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.svprep.reads.ReadRecord.getSoftClippedBases;

import static htsjdk.samtools.CigarOperator.M;

import java.util.EnumSet;

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
        // check each filter type that would prevent a read being used to establish a junction
        // NOTE: must check all filters (ie no early) exits since the state of some are subsequently used (eg map-qual & poly-G)
        final EnumSet<ReadFilterType> filters = EnumSet.noneOf(ReadFilterType.class);

        final Cigar cigar = record.getCigar();
        int alignedBases = cigar.getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

        if(alignedBases < mConfig.MinAlignmentBases)
            filters.add(ReadFilterType.MIN_ALIGN_MATCH);

        if(record.getMappingQuality() < mConfig.MinMapQuality)
            filters.add(ReadFilterType.MIN_MAP_QUAL);

        if(!mateUnmapped(record) && !record.getReadUnmappedFlag())
        {
            int insertAlignmentOverlap = abs(abs(record.getInferredInsertSize()) - alignedBases);

            if(insertAlignmentOverlap < mConfig.MinInsertAlignmentOverlap)
                filters.add(ReadFilterType.INSERT_MAP_OVERLAP);
        }

        // check length and quality of soft-clipped bases if not an INDEL
        int scLeft = leftSoftClipLength(record);
        int scRight = rightSoftClipLength(record);

        // a read with an indel junction does not need to meet the min soft-clip length condition, but other SC inserts are checked
        int maxIndelLength = ReadRecord.maxIndelLength(record.getCigar());

        if(maxIndelLength < mConfig.MinIndelLength && scLeft < mConfig.MinSoftClipLength && scRight < mConfig.MinSoftClipLength)
            filters.add(ReadFilterType.SOFT_CLIP_LENGTH);

        // base qual in soft clip
        if(scLeft > 0 || scRight > 0)
        {
            final byte[] baseQualities = record.getBaseQualities();
            int scRangeStart = scLeft > scRight ? 0 : baseQualities.length - scRight;
            int scRangeEnd = scLeft > scRight ? scLeft : baseQualities.length;
            boolean useLeftClip = scLeft > scRight;
            int scLength = useLeftClip ? scLeft : scRight;

            int aboveQual = 0;
            for(int i = scRangeStart; i < scRangeEnd; ++i)
            {
                if(baseQualities[i] >= mConfig.MinSoftClipHighQual)
                    ++aboveQual;
            }

            if(aboveQual / (double)scLength < mConfig.MinSoftClipHighQualPerc)
            {
                filters.add(ReadFilterType.SOFT_CLIP_BASE_QUAL);

                // additional check to exclude a read from any use
                int lowQualCount = scLength - aboveQual;

                if(lowQualCount >= MAX_SOFT_CLIP_LOW_QUAL_COUNT)
                    filters.add(ReadFilterType.SOFT_CLIP_LOW_BASE_QUAL);
            }

            String scBases = getSoftClippedBases(record, useLeftClip);

            // check for poly G/C inserts, then a break in a repetitive section, then poly A/T mobile line insertion
            if(scLength >= POLY_G_LENGTH && (scBases.contains(POLY_G_INSERT) || scBases.contains(POLY_C_INSERT)))
            {
                filters.add(ReadFilterType.POLY_G_SC);
            }
            else if(!filters.contains(SOFT_CLIP_LENGTH))
            {
                if((record.getMappingQuality() < REPEAT_BREAK_MIN_MAP_QUAL || scLength < REPEAT_BREAK_MIN_SC_LENGTH)
                && isRepetitiveSectionBreak(record.getReadBases(), useLeftClip, scLength))
                {
                    filters.add(ReadFilterType.BREAK_IN_REPEAT);
                }
            }
            else if(scLength >= MIN_LINE_SOFT_CLIP_LENGTH)
            {
                // make an exception if the soft-clip sequence meets the LINE criteria
                byte orientation = useLeftClip ? NEG_ORIENT : POS_ORIENT;

                if(isMobileLineElement(orientation, scBases)
                && !isRepetitiveSectionBreak(record.getReadBases(), useLeftClip, scLength))
                {
                    filters.remove(ReadFilterType.SOFT_CLIP_LENGTH);
                }
            }
        }

        return ReadFilterType.toFlags(filters);
    }

    public boolean isCandidateSupportingRead(final SAMRecord record, final int filters)
    {
        // exclude poly-G inserts
        if(ReadFilterType.isSet(filters, POLY_G_SC))
            return false;

        if(isChimericRead(record, mConfig))
            return true;

        // or with any amount of soft or hard clipping or a long INDEL
        return record.getCigar().isLeftClipped() || record.getCigar().isRightClipped()
                || ReadRecord.maxIndelLength(record.getCigar()) >= MIN_INDEL_SUPPORT_LENGTH;
    }

    public static boolean isChimericRead(final SAMRecord record, final ReadFilterConfig config)
    {
        if(record.getReadUnmappedFlag())
            return false;

        // or a fragment length outside the observed distribution
        if(abs(record.getInferredInsertSize()) > config.fragmentLengthMax())
            return true;

        // an unmapped mate
        if(mateUnmapped(record))
            return true;

        if(record.getReadPairedFlag())
        {
            // inter-chromosomal
            if(!record.getReferenceName().equals(record.getMateReferenceName()))
                return true;

            // inversion
            if(record.getReadNegativeStrandFlag() == mateNegativeStrand(record))
                return true;
        }

        return false;
    }

    public static boolean isRepetitiveSectionBreak(final byte[] readBases, boolean leftClipped, int softClippedLength)
    {
        int readLength = readBases.length;

        int startIndex;

        if(leftClipped)
        {
            // 0-9 sc bases, length = 10, checked range is 4-9 and 10-17
            startIndex = max(softClippedLength - REPEAT_BREAK_CHECK_LENGTH, 0);
        }
        else
        {
            // 91-100 sc bases, length = 10, checked range 83-90 and 91-96
            startIndex = max(readLength - softClippedLength - REPEAT_BREAK_CHECK_LENGTH, 0);
        }

        int endIndex = min(startIndex + REPEAT_BREAK_CHECK_LENGTH + REPEAT_BREAK_CHECK_LENGTH, readLength);

        if(endIndex - startIndex < REPEAT_BREAK_CHECK_LENGTH + REPEAT_BREAK_CHECK_LENGTH)
            return false;

        byte firstBase = readBases[startIndex];

        boolean allMatch = true;

        for(int i = startIndex + 1; i < endIndex; ++i)
        {
            if(readBases[i] != firstBase)
            {
                allMatch = false;
                break;
            }
        }

        if(allMatch)
            return true;

        byte secondBase = readBases[startIndex + 1];
        allMatch = true;

        for(int i = startIndex + 2; i < endIndex - 1; i += 2)
        {
            if(readBases[i] != firstBase || readBases[i + 1] != secondBase)
            {
                allMatch = false;
                break;
            }
        }

        if(allMatch)
            return true;

        byte thirdBase = readBases[startIndex + 2];
        allMatch = true;

        for(int i = startIndex + 3; i < endIndex - 2; i += 3)
        {
            if(readBases[i] != firstBase || readBases[i + 1] != secondBase || readBases[i + 2] != thirdBase)
            {
                allMatch = false;
                break;
            }
        }

        return allMatch;
    }
}
