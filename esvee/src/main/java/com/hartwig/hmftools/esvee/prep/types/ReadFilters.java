package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.maxIndelLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateUnmapped;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.region.ExcludedRegions.POLY_C_INSERT;
import static com.hartwig.hmftools.common.region.ExcludedRegions.POLY_G_INSERT;
import static com.hartwig.hmftools.common.region.ExcludedRegions.POLY_G_LENGTH;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_REQ;
import static com.hartwig.hmftools.common.sv.LineElements.isMobileLineElement;
import static com.hartwig.hmftools.common.utils.Arrays.copyArray;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.calcTrimmedBaseLength;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isDiscordantFragment;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_SUPPORT_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MAX_SOFT_CLIP_LOW_QUAL_COUNT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_LINE_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.REPEAT_BREAK_CHECK_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.REPEAT_BREAK_MIN_MAP_QUAL;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.REPEAT_BREAK_MIN_SC_LENGTH;

import static htsjdk.samtools.CigarOperator.M;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class ReadFilters
{
    private final ReadFilterConfig mConfig;

    public ReadFilters(final ReadFilterConfig config)
    {
        mConfig = config;
    }

    public ReadFilterConfig config() { return mConfig; }

    private static final int MIN_CHECK_SC_BASES = min(LINE_POLY_AT_REQ, POLY_G_LENGTH);

    public int checkFilters(final SAMRecord record)
    {
        // check each filter type that would prevent a read being used to establish a junction
        // NOTE: must check all filters (ie no early) exits since the state of some are subsequently used (eg map-qual & poly-G)
        int filters = 0;

        final Cigar cigar = record.getCigar();
        int alignedBases = cigar.getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

        if(alignedBases < mConfig.MinAlignmentBases)
            filters = ReadFilterType.set(filters, ReadFilterType.MIN_ALIGN_MATCH);

        if(record.getMappingQuality() < mConfig.MinMapQuality)
            filters = ReadFilterType.set(filters, ReadFilterType.MIN_MAP_QUAL);

        if(!mateUnmapped(record) && !record.getReadUnmappedFlag())
        {
            int insertAlignmentOverlap = abs(abs(record.getInferredInsertSize()) - alignedBases);

            if(insertAlignmentOverlap < mConfig.MinInsertAlignmentOverlap)
                filters = ReadFilterType.set(filters, ReadFilterType.INSERT_MAP_OVERLAP);
        }

        // check length and quality of soft-clipped bases if not an INDEL
        int scLeft = leftSoftClipLength(record);
        int scRight = rightSoftClipLength(record);

        // a read with an indel junction does not need to meet the min soft-clip length condition, but other SC inserts are checked
        int maxIndelLength = maxIndelLength(record.getCigar().getCigarElements());

        if(maxIndelLength < mConfig.MinIndelLength && scLeft < mConfig.MinSoftClipLength && scRight < mConfig.MinSoftClipLength)
            filters = ReadFilterType.set(filters, ReadFilterType.SOFT_CLIP_LENGTH);

        // base qual in soft clip
        if(scLeft > 0 || scRight > 0)
        {
            final byte[] baseQualities = record.getBaseQualities();
            boolean useLeftClip = scLeft > scRight;
            int scRangeStart = useLeftClip ? 0 : baseQualities.length - scRight;
            int scRangeEnd = useLeftClip ? scLeft : baseQualities.length;
            int scLength = useLeftClip ? scLeft : scRight;

            int aboveQual = 0;
            StringBuilder scStr = scLength >= MIN_CHECK_SC_BASES ? new StringBuilder() : null;
            for(int i = scRangeStart; i < scRangeEnd; ++i)
            {
                if(baseQualities[i] >= mConfig.MinSoftClipHighQual)
                    ++aboveQual;

                if(scStr != null)
                    scStr.append((char)record.getReadBases()[i]);
            }

            if(aboveQual / (double)scLength < mConfig.MinSoftClipHighQualPerc)
            {
                filters = ReadFilterType.set(filters, ReadFilterType.SOFT_CLIP_BASE_QUAL);

                // additional check to exclude a read from any use
                int lowQualCount = scLength - aboveQual;

                if(lowQualCount >= MAX_SOFT_CLIP_LOW_QUAL_COUNT)
                    filters = ReadFilterType.set(filters, ReadFilterType.SOFT_CLIP_LOW_BASE_QUAL);
            }

            String scBases = scStr != null ? scStr.toString() : "";

            // check for poly G/C inserts, then a break in a repetitive section, then poly A/T mobile line insertion
            if(scLength >= POLY_G_LENGTH && (scBases.contains(POLY_G_INSERT) || scBases.contains(POLY_C_INSERT)))
            {
                filters = ReadFilterType.set(filters, ReadFilterType.POLY_G_SC);
            }
            else if(!ReadFilterType.isSet(filters, ReadFilterType.SOFT_CLIP_LENGTH))
            {
                if((record.getMappingQuality() < REPEAT_BREAK_MIN_MAP_QUAL || scLength < REPEAT_BREAK_MIN_SC_LENGTH)
                && isRepetitiveSectionBreak(record.getReadBases(), useLeftClip, scLength))
                {
                    filters = ReadFilterType.set(filters, ReadFilterType.BREAK_IN_REPEAT);
                }
            }
            else if(scLength >= MIN_LINE_SOFT_CLIP_LENGTH)
            {
                // make an exception if the soft-clip sequence meets the LINE criteria
                Orientation orientation = useLeftClip ? REVERSE : FORWARD;

                if(isMobileLineElement(orientation.asByte(), scBases)
                && !isRepetitiveSectionBreak(record.getReadBases(), useLeftClip, scLength))
                {
                    filters = ReadFilterType.unset(filters, ReadFilterType.SOFT_CLIP_LENGTH);
                }
            }
        }

        return filters;
    }

    public boolean isCandidateSupportingRead(final SAMRecord record, final int filters)
    {
        // exclude poly-G inserts
        if(ReadFilterType.isSet(filters, ReadFilterType.POLY_G_SC))
            return false;

        if(isChimericRead(record, mConfig))
            return true;

        // or with any amount of soft or hard clipping or a long INDEL
        return record.getCigar().isLeftClipped() || record.getCigar().isRightClipped()
                || maxIndelLength(record.getCigar().getCigarElements()) >= MIN_INDEL_SUPPORT_LENGTH;
    }

    public static boolean isChimericRead(final SAMRecord record, final ReadFilterConfig config)
    {
        // only difference from the common discordant fragment method is that reads with unmapped mates are also captured
        // an unmapped mate
        if(mateUnmapped(record))
            return true;

        return isDiscordantFragment(record, config.fragmentLengthMax(), null);
    }

    public static boolean isRepetitiveSectionBreak(final byte[] readBases, boolean leftClipped, int scLength)
    {
        int readLength = readBases.length;

        int startIndex;

        if(leftClipped)
        {
            // 0-9 sc bases, length = 10, checked range is 4-9 and 10-17
            startIndex = max(scLength - REPEAT_BREAK_CHECK_LENGTH, 0);
        }
        else
        {
            // 91-100 sc bases, length = 10, checked range 83-90 and 91-96
            startIndex = max(readLength - scLength - REPEAT_BREAK_CHECK_LENGTH, 0);
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

    public static boolean aboveRepeatTrimmedAlignmentThreshold(final PrepRead read, final int minLength)
    {
        int readIndex = 0;

        int alignedLength = read.cigar().getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

        if(alignedLength == 0)
            return false;

        byte[] alignedBases = new byte[alignedLength];
        int alignedIndex = 0;

        for(CigarElement element : read.cigar().getCigarElements())
        {
            if(element.getOperator() == M)
            {
                copyArray(read.record().getReadBases(), alignedBases, readIndex, readIndex + element.getLength(), alignedIndex);
                alignedIndex += element.getLength();
            }

            if(element.getOperator().consumesReadBases())
                readIndex += element.getLength();
        }

        List<RepeatInfo> repeats = RepeatInfo.findRepeats(alignedBases);

        int repeatTrimmedLength = calcTrimmedBaseLength(0, alignedBases.length - 1, repeats);

        return repeatTrimmedLength >= minLength;
    }
}
