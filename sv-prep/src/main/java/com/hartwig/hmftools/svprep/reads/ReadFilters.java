package com.hartwig.hmftools.svprep.reads;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.svprep.SvConstants.REPEAT_BREAK_MATCH_CHECK_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.REPEAT_BREAK_MIN_MAP_QUAL;
import static com.hartwig.hmftools.svprep.SvConstants.REPEAT_BREAK_MIN_SC_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.REPEAT_BREAK_SC_CHECK_LENGTH;

import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.binaryToEnum;
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

        if(!record.getMateUnmappedFlag() && !record.getReadUnmappedFlag())
        {
            int insertAlignmentOverlap = abs(abs(record.getInferredInsertSize()) - alignedBases);

            if(insertAlignmentOverlap < mConfig.MinInsertAlignmentOverlap)
                filters |= ReadFilterType.INSERT_MAP_OVERLAP.flag();
        }

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
                boolean useLeftClip = scLeft > scRight;
                int scLength = useLeftClip ? scLeft : scRight;

                int aboveQual = 0;
                for(int i = scRangeStart; i < scRangeEnd; ++i)
                {
                    if(baseQualities[i] >= mConfig.MinSoftClipHighQual)
                        ++aboveQual;
                }

                if(aboveQual / scLength < mConfig.MinSoftClipHighQualPerc)
                    filters |= ReadFilterType.SOFT_CLIP_BASE_QUAL.flag();

                if(!ReadFilterType.isSet(filters, ReadFilterType.SOFT_CLIP_LENGTH))
                {
                    if((record.getMappingQuality() < REPEAT_BREAK_MIN_MAP_QUAL || scLength < REPEAT_BREAK_MIN_SC_LENGTH)
                    && isRepetitiveSectionBreak(record.getReadBases(), useLeftClip, scLength))
                    {
                        filters |= ReadFilterType.BREAK_IN_REPEAT.flag();
                    }
                }
            }
        }

        return filters;
    }

    public boolean isCandidateSupportingRead(final SAMRecord record)
    {
        if(isChimericRead(record, mConfig))
            return true;

        // or with any amount of soft-clipping
        return record.getCigar().isLeftClipped() || record.getCigar().isRightClipped();
    }

    public static boolean isChimericRead(final SAMRecord record, final ReadFilterConfig config)
    {
        if(record.getReadUnmappedFlag())
            return false;

        // any read with a supplementary
        if(record.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
            return true;

        // or an fragment length beyond the observed distribution
        if(abs(record.getInferredInsertSize()) > config.fragmentLengthMax()) //  || record.getInferredInsertSize() < mFragmentLengthMin
            return true;

        // an unmapped mate
        if(record.getMateUnmappedFlag())
            return true;

        // interchromosomal
        if(!record.getReferenceName().equals(record.getMateReferenceName()))
            return true;

        // inversion
        if(record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag())
            return true;

        return false;
    }

    public static boolean isRepetitiveSectionBreak(final byte[] readBases, boolean leftClipped, int scLength)
    {
        int readLength = readBases.length;

        int startIndex;

        if(leftClipped)
        {
            // 0-9 sc bases, length = 10, checked range is 4-9 and 10-17
            startIndex = max(scLength - REPEAT_BREAK_SC_CHECK_LENGTH, 0);
        }
        else
        {
            // 91-100 sc bases, length = 10, checked range 83-90 and 91-96
            startIndex = max(readLength - scLength - REPEAT_BREAK_MATCH_CHECK_LENGTH, 0);
        }

        int endIndex = min(startIndex + REPEAT_BREAK_SC_CHECK_LENGTH + REPEAT_BREAK_MATCH_CHECK_LENGTH, readLength);

        if(endIndex - startIndex < REPEAT_BREAK_SC_CHECK_LENGTH + REPEAT_BREAK_MATCH_CHECK_LENGTH)
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
