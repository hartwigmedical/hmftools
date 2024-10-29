package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.maxIndelLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateUnmapped;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.region.ExcludedRegions.POLY_C_INSERT;
import static com.hartwig.hmftools.common.region.ExcludedRegions.POLY_G_INSERT;
import static com.hartwig.hmftools.common.region.ExcludedRegions.POLY_G_LENGTH;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_REQ;
import static com.hartwig.hmftools.common.sv.LineElements.isMobileLineElement;
import static com.hartwig.hmftools.common.utils.Arrays.copyArray;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isDiscordantFragment;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_SUPPORT_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MAX_SOFT_CLIP_LOW_QUAL_COUNT;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_LINE_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.REPEAT_BREAK_CHECK_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.REPEAT_BREAK_MIN_MAP_QUAL;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.REPEAT_BREAK_MIN_SC_LENGTH;

import static htsjdk.samtools.CigarOperator.M;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterConfig;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterType;

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

    public boolean ignoreRead(final SAMRecord record)
    {
        // ignore reads with the majority of bases low qual - the assembly routine has this same check for now
        if(filterLowQualRead(record))
            return true;

        // ignore local, non-chimeric, locally-aligned reads
        if(record.getSupplementaryAlignmentFlag() || record.getCigar().getCigarElements().size() > 1)
            return false;

        if(record.getReadUnmappedFlag() || mateUnmapped(record))
            return false;

        if(isChimericRead(record, mConfig))
            return false;

        String mateCigarStr = record.getStringAttribute(MATE_CIGAR_ATTRIBUTE);

        if(mateCigarStr != null && mateCigarStr.equals(format("%dM", record.getReadBases().length)))
            return true;

        return false;
    }

    public static boolean filterLowQualRead(final SAMRecord read)
    {
        // filter any read with 50% + bases classified as low qual
        int baseLength = read.getReadBases().length;
        int qualCountThreshold = baseLength / 2 + 1;
        int lowQualCount = 0;

        for(int i = 0; i < baseLength; ++i)
        {
            if(belowMinQual(read.getBaseQualities()[i]))
            {
                ++lowQualCount;

                if(lowQualCount >= qualCountThreshold)
                    return true;
            }
            else
            {
                // exit early if majority will be high-qual
                int highQualCount = i + 1 - lowQualCount;
                if(highQualCount >= qualCountThreshold)
                    return false;
            }
        }

        return false;
    }

    public boolean isCandidateSupportingRead(final PrepRead read)
    {
        if(isChimericRead(read.record(), mConfig))
            return true;

        // or with any amount of soft or hard clipping or a long INDEL
        return read.isLeftClipped() || read.isRightClipped() || read.maxIndelLength() >= MIN_INDEL_SUPPORT_LENGTH;
    }

    public static boolean isChimericRead(final SAMRecord record, final ReadFilterConfig config)
    {
        // only difference from the common discordant fragment method is that reads with unmapped mates are also captured
        // an unmapped mate
        if(mateUnmapped(record))
            return true;

        return isDiscordantFragment(record, config.fragmentLengthMax(), null);
    }

    public void checkFilters(final PrepRead read)
    {
        // check each filter type that would prevent a read being used to establish a junction
        // NOTE: must check all filters (ie no early) exits since the state of some are subsequently used (eg map-qual & poly-G)
        int alignedBases = read.alignedBaseLength();

        SAMRecord record = read.record();

        if(alignedBases < mConfig.MinAlignmentBases)
            read.addFilter(ReadFilterType.MIN_ALIGN_MATCH);

        if(read.record().getMappingQuality() < mConfig.MinMapQuality)
            read.addFilter(ReadFilterType.MIN_MAP_QUAL);

        if(!mateUnmapped(record) && !record.getReadUnmappedFlag())
        {
            int insertAlignmentOverlap = abs(abs(record.getInferredInsertSize()) - alignedBases);

            if(insertAlignmentOverlap < mConfig.MinInsertAlignmentOverlap)
                read.addFilter(ReadFilterType.INSERT_MAP_OVERLAP);
        }

        // check length and quality of soft-clipped bases if not an INDEL
        int scLeft = read.leftClipLength();
        int scRight = read.rightClipLength();

        // a read with an indel junction does not need to meet the min soft-clip length condition, but other SC inserts are checked
        if(read.maxIndelLength() < mConfig.MinIndelLength)
        {
            if(scLeft < mConfig.MinSoftClipLength && scRight < mConfig.MinSoftClipLength)
                read.addFilter(ReadFilterType.SOFT_CLIP_LENGTH);
        }

        // base qual in soft clip
        if(scLeft > 0 || scRight > 0)
        {
            final byte[] baseQualities = record.getBaseQualities();
            boolean useLeftClip = scLeft > scRight;
            int scRangeStart = useLeftClip ? 0 : baseQualities.length - scRight;
            int scRangeEnd = useLeftClip ? scLeft : baseQualities.length;
            int scLength = useLeftClip ? scLeft : scRight;

            int aboveQual = 0;
            byte[] scBaseArray = scLength >= MIN_CHECK_SC_BASES ? new byte[scLength] : null;

            int scIndex = 0;
            for(int i = scRangeStart; i < scRangeEnd; ++i)
            {
                if(baseQualities[i] >= mConfig.MinSoftClipHighQual)
                    ++aboveQual;

                if(scBaseArray != null)
                    scBaseArray[scIndex++] = record.getReadBases()[i];
            }

            if(aboveQual / (double)scLength < mConfig.MinSoftClipHighQualPerc)
            {
                read.addFilter(ReadFilterType.SOFT_CLIP_BASE_QUAL);

                // additional check to exclude a read from any use
                int lowQualCount = scLength - aboveQual;

                if(lowQualCount >= MAX_SOFT_CLIP_LOW_QUAL_COUNT)
                    read.addFilter(ReadFilterType.SOFT_CLIP_LOW_BASE_QUAL);
            }

            String scBases = scBaseArray != null ? new String(scBaseArray) : "";

            // check for poly G/C inserts, then a break in a repetitive section, then poly A/T mobile line insertion
            if(scLength >= POLY_G_LENGTH && (scBases.contains(POLY_G_INSERT) || scBases.contains(POLY_C_INSERT)))
            {
                read.addFilter(ReadFilterType.POLY_G_SC);
            }
            else if(!read.hasFilter(ReadFilterType.SOFT_CLIP_LENGTH))
            {
                if((record.getMappingQuality() < REPEAT_BREAK_MIN_MAP_QUAL || scLength < REPEAT_BREAK_MIN_SC_LENGTH)
                && isRepetitiveSectionBreak(record.getReadBases(), useLeftClip, scLength))
                {
                    read.addFilter(ReadFilterType.BREAK_IN_REPEAT);
                }
            }
            else if(scLength >= MIN_LINE_SOFT_CLIP_LENGTH)
            {
                // make an exception if the soft-clip sequence meets the LINE criteria
                Orientation orientation = useLeftClip ? REVERSE : FORWARD;

                if(isMobileLineElement(orientation.asByte(), scBases)
                && !isRepetitiveSectionBreak(record.getReadBases(), useLeftClip, scLength))
                {
                    read.removefilter(ReadFilterType.SOFT_CLIP_LENGTH);
                }
            }
        }
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

    public static boolean aboveRepeatTrimmedAlignmentThreshold(final PrepRead read, final int minLength, boolean applyPowerAdjust)
    {
        // at least one junction supporting read with [AS - Len + TrimmedLen] * (AS/SUM(M))^2 > 50
        int readIndex = 0;

        int alignedLength = read.alignedBaseLength();

        if(alignedLength < minLength)
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

        int repeatReduction = 0;

        if(repeats != null)
        {
            for(RepeatInfo repeatInfo : repeats)
            {
                // each repeat only counts twice at most and other bases will detract the aligned score
                repeatReduction += max(repeatInfo.Count - 2, 0) * repeatInfo.Bases.length();
            }
        }

        int alignmentScore = 0;

        if(read.record().hasAttribute(ALIGNMENT_SCORE_ATTRIBUTE))
        {
            alignmentScore = read.record().getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE).intValue();
        }
        else
        {
            alignmentScore = alignedLength;

            if(read.record().hasAttribute(NUM_MUTATONS_ATTRIBUTE))
                alignmentScore -= read.record().getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue();
        }

        double adjustedAlignScore = (alignmentScore - repeatReduction);

        if(applyPowerAdjust)
            adjustedAlignScore *= pow(alignmentScore / (double)alignedLength, 2);

        return adjustedAlignScore >= minLength;
    }
}
