package com.hartwig.hmftools.sage.candidate;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.NO_POSITION_INFO;
import static com.hartwig.hmftools.common.bam.CigarUtils.getPositionFromReadIndex;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.SageConstants.MIN_INSERT_ALIGNMENT_OVERLAP;
import static com.hartwig.hmftools.sage.SageConstants.REGION_BLOCK_SIZE;
import static com.hartwig.hmftools.sage.SageConstants.SC_INSERT_REF_TEST_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.SC_INSERT_MIN_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.SC_READ_EVENTS_FACTOR;
import static com.hartwig.hmftools.sage.common.Microhomology.findLeftHomologyShift;
import static com.hartwig.hmftools.sage.common.NumberEvents.rawNM;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.isImproperPair;

import java.util.List;
import java.util.Set;
import java.util.function.Function;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.bam.CigarHandler;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.common.NumberEvents;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class RefContextConsumer
{
    private final SageConfig mConfig;
    private final ChrBaseRegion mBounds;
    private final RefSequence mRefSequence;
    private final RefContextCache mRefContextCache;
    private final VariantReadContextBuilder mReadContextBuilder;
    private final Set<Integer> mHotspotPositions;

    private final List<RegionBlock> mRegionBlocks;
    private int mCurrentRegionBlockIndex;
    private RegionBlock mCurrentRegionBlock;
    private RegionBlock mNextRegionBlock;

    private int mReadCount;

    public RefContextConsumer(
            final SageConfig config, final ChrBaseRegion regionBounds, final RefSequence refSequence, final RefContextCache refContextCache,
            final List<SimpleVariant> regionHotspots)
    {
        mBounds = regionBounds;
        mRefSequence = refSequence;
        mRefContextCache = refContextCache;
        mReadContextBuilder = new VariantReadContextBuilder(config.ReadContextFlankLength);
        mConfig = config;

        mHotspotPositions = Sets.newHashSet();
        regionHotspots.forEach(x -> mHotspotPositions.add(x.position()));

        mRegionBlocks = RegionBlock.buildRegionBlocks(
                REGION_BLOCK_SIZE, regionBounds, mRefContextCache.panelSelector(), regionHotspots,
                mConfig.MaxReadDepthPanel, mConfig.MaxReadDepth);

        mCurrentRegionBlockIndex = 0;
        mCurrentRegionBlock = mRegionBlocks.get(0);
        mNextRegionBlock = null;

        mReadCount = 0;
    }

    public int getReadCount() { return mReadCount; }

    private void updateCurrentRegionBlocks(int readStart)
    {
        if(mCurrentRegionBlock.containsPosition(readStart))
            return;

        while(readStart > mCurrentRegionBlock.end() && mCurrentRegionBlockIndex < mRegionBlocks.size())
        {
            ++mCurrentRegionBlockIndex;
            mCurrentRegionBlock = mRegionBlocks.get(mCurrentRegionBlockIndex);
        }

        if(mCurrentRegionBlockIndex < mRegionBlocks.size() - 1)
            mNextRegionBlock = mRegionBlocks.get(mCurrentRegionBlockIndex + 1);
    }

    public void processRead(final SAMRecord record)
    {
        int readStart = record.getAlignmentStart();
        int readEnd = record.getAlignmentEnd();

        if(!positionsOverlap(readStart, readEnd, mBounds.start(), mBounds.end()))
            return;

        updateCurrentRegionBlocks(readStart);

        ++mReadCount;

        if(reachedDepthLimit(readStart) && reachedDepthLimit(readEnd))
            return;

        ReadInfo readInfo = buildReadInfo(record);

        int scEvents = (int)NumberEvents.calcSoftClipAdjustment(readInfo.SoftClipLength);
        int adjustedMapQual = calcAdjustedMapQualLessEventsPenalty(
                record, readInfo.NumberOfEvents, applyMapQualEventPenalty(readStart, readEnd));

        readInfo.ReadExceedsQuality = adjustedMapQual > 0;

        boolean readCoversHotspot = readCoversHotspot(readStart, readEnd);

        // if the read is below the required threshold and does not cover a hotspot position, then stop processing it
        if(!readInfo.ReadExceedsQuality && !readCoversHotspot)
            return;

        updateRegionBlockDepth(readStart, readEnd);

        int scAdjustedMapQual = (int)round(adjustedMapQual - scEvents * mConfig.Quality.ReadMapQualEventsPenalty);
        readInfo.ReadExceedsScAdjustedQuality = scAdjustedMapQual > 0;
        boolean ignoreScAdapter = scEvents > 0 && ignoreSoftClipAdapter(record, readInfo.AlignedLength);

        final CigarHandler handler = new CigarHandler()
        {
            @Override
            public void handleAlignment(
                    final SAMRecord record, final CigarElement element, final int readIndex, final int refPosition)
            {
                if(record.isSecondaryOrSupplementary())
                    return;

                if(!readInfo.ReadExceedsScAdjustedQuality && !readCoversHotspot)
                    return;

                processAlignment(record, readIndex, refPosition, element.getLength(), readInfo);
            }

            @Override
            public void handleInsert(final SAMRecord record, final CigarElement element, final int readIndex, final int refPosition)
            {
                if(record.isSecondaryOrSupplementary())
                    return;

                processInsert(element, record, readIndex, refPosition, readInfo);
            }

            @Override
            public void handleDelete(final SAMRecord record, final CigarElement element, final int readIndex, final int refPosition)
            {
                if(record.isSecondaryOrSupplementary())
                    return;

                processDel(element, record, readIndex, refPosition, readInfo);
            }

            @Override
            public void handleLeftSoftClip(final SAMRecord record, final CigarElement element)
            {
                if(ignoreScAdapter)
                    return;

                processSoftClip(record, element.getLength(), 0, mRefSequence, readInfo, true);
            }

            @Override
            public void handleRightSoftClip(final SAMRecord record, final CigarElement element, int readIndex, int refPosition)
            {
                if(ignoreScAdapter)
                    return;

                processSoftClip(record, element.getLength(), readIndex, mRefSequence, readInfo, false);
            }
        };

        CigarHandler.traverseCigar(record, handler);

        if(readInfo.AltReads != null)
        {
            for(AltRead altRead : readInfo.AltReads)
            {
                altRead.updateRefContext(mReadContextBuilder, mRefSequence);
            }
        }
    }

    private class ReadInfo
    {
        public int NumberOfEvents;
        public int AlignedLength;
        public int SoftClipLength;
        public boolean ReadExceedsQuality;
        public boolean ReadExceedsScAdjustedQuality;

        public List<AltRead> AltReads;

        public ReadInfo()
        {
            AlignedLength = 0;
            SoftClipLength = 0;
            NumberOfEvents = 0;
            ReadExceedsQuality = false;
            ReadExceedsScAdjustedQuality = false;
            AltReads = null;
        }

        public void addAltRead(final AltRead altRead)
        {
            if(AltReads == null)
                AltReads = Lists.newArrayList();

            AltReads.add(altRead);
        }
    }

    private ReadInfo buildReadInfo(final SAMRecord record)
    {
        ReadInfo readInfo = new ReadInfo();

        int additionalIndels = 0;

        for(CigarElement cigarElement : record.getCigar())
        {
            switch(cigarElement.getOperator())
            {
                case M:
                    readInfo.AlignedLength += cigarElement.getLength();
                    break;

                case S:
                    readInfo.SoftClipLength += cigarElement.getLength();
                    break;

                case D:
                case I:
                    additionalIndels += cigarElement.getLength() - 1;
                    break;
            }
        }

        if(!record.getSupplementaryAlignmentFlag())
        {
            readInfo.NumberOfEvents = max(rawNM(record, mRefSequence) - additionalIndels, 0);
        }

        return readInfo;
    }

    private static boolean ignoreSoftClipAdapter(final SAMRecord record, int alignedLength)
    {
        // ignore soft-clips from short fragments indicating adapter bases are present
        int fragmentLength = abs(record.getInferredInsertSize());

        if(fragmentLength >= record.getReadBases().length + MIN_INSERT_ALIGNMENT_OVERLAP)
            return false;

        int insertAlignmentOverlap = abs(fragmentLength - alignedLength);
        return insertAlignmentOverlap < MIN_INSERT_ALIGNMENT_OVERLAP;
    }

    private boolean reachedDepthLimit(int position)
    {
        if(mCurrentRegionBlock.containsPosition(position))
            return mCurrentRegionBlock.depthLimitReached();

        if(mNextRegionBlock != null && mNextRegionBlock.containsPosition(position))
            return mNextRegionBlock.depthLimitReached();

        return false;
    }

    private void updateRegionBlockDepth(int readStart, int readEnd)
    {
        if(mCurrentRegionBlock.containsPosition(readStart))
            mCurrentRegionBlock.incrementDepth();

        if(mNextRegionBlock != null && mNextRegionBlock.containsPosition(readEnd))
            mNextRegionBlock.incrementDepth();
    }

    private boolean readCoversHotspot(int readStart, int readEnd)
    {
        if(mCurrentRegionBlock.coversHotspot(readStart, readEnd))
            return true;

        return mNextRegionBlock != null && mNextRegionBlock.coversHotspot(readStart, readEnd);
    }

    private boolean applyMapQualEventPenalty(int readStart, int readEnd)
    {
        // if the read overlaps any block where the event penalty should not be applied, then don't apply it
        if(mCurrentRegionBlock.containsPosition(readStart) && !mCurrentRegionBlock.applyEventPenalty())
            return false;

        if(mNextRegionBlock != null && mNextRegionBlock.containsPosition(readEnd) && !mNextRegionBlock.applyEventPenalty())
            return false;

        return true;
    }

    private int calcAdjustedMapQualLessEventsPenalty(final SAMRecord record, int numberOfEvents, boolean applyEventPenalty)
    {
        if(!applyEventPenalty)
            return record.getMappingQuality();;

        int eventPenalty = (int)round((numberOfEvents - 1) * mConfig.Quality.ReadMapQualEventsPenalty);

        int improperPenalty = isImproperPair(record) || record.getSupplementaryAlignmentFlag() ?
                mConfig.Quality.ImproperPairPenalty : 0;

        return record.getMappingQuality() - mConfig.Quality.FixedMapQualPenalty - eventPenalty - improperPenalty;
    }

    private boolean isHotspotPosition(int position)
    {
        return mCurrentRegionBlock.isHotspot(position) || mNextRegionBlock != null && mNextRegionBlock.isHotspot(position);
    }

    private void processInsert(
            final CigarElement element, final SAMRecord record, int readIndex, int refPosition, final ReadInfo readInfo)
    {
        if(!mBounds.containsPosition(refPosition))
            return;

        boolean exceedsQuality = element.getLength() <= SC_READ_EVENTS_FACTOR ?
                readInfo.ReadExceedsScAdjustedQuality : readInfo.ReadExceedsQuality;

        if(!exceedsQuality && !isHotspotPosition(refPosition))
            return;

        int refIndex = mRefSequence.index(refPosition);
        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        String ref = new String(mRefSequence.Bases, refIndex, 1);
        String alt = new String(record.getReadBases(), readIndex, element.getLength() + 1);

        RefContext refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);

        readInfo.addAltRead(new AltRead(refContext, ref, alt, readInfo.NumberOfEvents, sufficientMapQuality, record, readIndex));
    }

    private void processDel(
            final CigarElement element, final SAMRecord record, int readIndex, int refPosition, final ReadInfo readInfo)
    {
        if(!mBounds.containsPosition(refPosition))
            return;

        boolean exceedsQuality = element.getLength() <= SC_READ_EVENTS_FACTOR ?
                readInfo.ReadExceedsScAdjustedQuality : readInfo.ReadExceedsQuality;

        if(!exceedsQuality && !isHotspotPosition(refPosition))
            return;

        int refIndex = mRefSequence.index(refPosition);
        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        final String ref = new String(mRefSequence.Bases, refIndex, element.getLength() + 1);
        final String alt = new String(record.getReadBases(), readIndex, 1);

        final RefContext refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);

        if(refContext != null)
        {
            readInfo.addAltRead(new AltRead(refContext, ref, alt, readInfo.NumberOfEvents, sufficientMapQuality, record, readIndex));
        }
    }

    private void processAlignment(
            final SAMRecord record, int readBasesStartIndex, int refPositionStart, int alignmentLength, final ReadInfo readInfo)
    {
        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        int refIndex = mRefSequence.index(refPositionStart);

        for(int i = 0; i < alignmentLength; i++)
        {
            int refPosition = refPositionStart + i;
            int readBaseIndex = readBasesStartIndex + i;
            int refBaseIndex = refIndex + i;

            if(refPosition < mBounds.start())
                continue;

            if(refPosition > mBounds.end())
                break;

            if(!readInfo.ReadExceedsQuality && !isHotspotPosition(refPosition))
                continue;

            byte refByte = mRefSequence.Bases[refBaseIndex];
            byte readByte = record.getReadBases()[readBaseIndex];

            if(readByte != refByte)
            {
                final RefContext refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);
                if(refContext == null)
                    continue;

                String alt = String.valueOf((char) readByte);
                String ref = String.valueOf((char) refByte);

                readInfo.addAltRead(new AltRead(refContext, ref, alt, readInfo.NumberOfEvents, sufficientMapQuality, record, readBaseIndex));

                int mnvMaxLength = mnvLength(readBaseIndex, refBaseIndex, record.getReadBases(), mRefSequence.Bases);

                int nextReadIndex = i;
                for(int mnvLength = 2; mnvLength <= mnvMaxLength; mnvLength++)
                {
                    ++nextReadIndex;

                    // MNVs cannot extend past the end of this Cigar element
                    if(nextReadIndex >= alignmentLength)
                        break;

                    String mnvRef = new String(mRefSequence.Bases, refBaseIndex, mnvLength);
                    String mnvAlt = new String(record.getReadBases(), readBaseIndex, mnvLength);

                    // Only check last base because some subsets may not be valid,
                    // ie CA > TA is not a valid subset of CAC > TAT
                    if(mnvRef.charAt(mnvLength - 1) != mnvAlt.charAt(mnvLength - 1))
                    {
                        readInfo.addAltRead(new AltRead(
                                refContext, mnvRef, mnvAlt, NumberEvents.calcWithMnvRaw(readInfo.NumberOfEvents, mnvRef, mnvAlt),
                                sufficientMapQuality, record, readBaseIndex));
                    }
                }
            }
        }
    }

    private void processSoftClip(
            final SAMRecord record, int scLength, int scReadIndex, final RefSequence refSequence, final ReadInfo readInfo, boolean onLeft)
    {
        if(!readInfo.ReadExceedsQuality)
            return;

        if(scLength < SC_INSERT_REF_TEST_LENGTH + 1)
            return;

        AltRead altRead = processSoftClip(
                record.getAlignmentStart(), record.getAlignmentEnd(), record.getReadString(), scLength, scReadIndex, refSequence, onLeft);

        if(altRead == null)
            return;

        int refPosition;
        int readIndex;

        if(onLeft)
        {
            refPosition = record.getAlignmentStart() - 1;

            // set to start at the ref/alt base prior to the insert
            readIndex = scLength - altRead.Alt.length();
        }
        else
        {
            refPosition = record.getAlignmentEnd();
            readIndex = record.getReadBases().length - scLength - 1;
        }

        if(!mBounds.containsPosition(refPosition))
            return;

        if(!withinReadContext(readIndex, record))
            return;

        SimpleVariant variant = new SimpleVariant(record.getContig(), refPosition, altRead.Ref, altRead.Alt);

        if(variant.isInsert() && !onLeft)
        {
            int leftHomologyShift = findLeftHomologyShift(variant, mRefSequence, record.getReadBases(), readIndex);

            if(leftHomologyShift > 0)
            {
                int newReadIndex = readIndex - leftHomologyShift;

                // recompute the new reference position, taking into consideration indels in the read
                int[] posInfo = getPositionFromReadIndex(
                        record.getAlignmentStart(), record.getCigar().getCigarElements(), newReadIndex, true, true);

                if(posInfo == NO_POSITION_INFO) // don't revert to original variant if cannot align to the new index
                    return;

                refPosition = posInfo[0];
                readIndex = newReadIndex + posInfo[1];

                String newAltBases, newRefBases;

                if(variant.isInsert())
                {
                    newRefBases = mRefSequence.positionBases(refPosition, refPosition);
                    newAltBases = newRefBases + new String(Arrays.subsetArray(
                            record.getReadBases(), readIndex + 1, readIndex + variant.altLength() - 1));
                }
                else
                {
                    newRefBases = mRefSequence.positionBases(refPosition, refPosition + variant.refLength() - 1);
                    newAltBases = newRefBases.substring(0, 1);
                }

                variant = new SimpleVariant(record.getContig(), refPosition, newRefBases, newAltBases);
            }
        }

        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        RefContext refContext = mRefContextCache.getOrCreateRefContext(variant.Chromosome, variant.Position);

        AltRead altReadFull = new AltRead(refContext, variant.Ref, variant.Alt, readInfo.NumberOfEvents, sufficientMapQuality, record, readIndex);
        readInfo.addAltRead(altReadFull);
    }

    public static AltRead processSoftClip(
            int readStart, int readEnd, final String readBases, int scLength, int scReadIndex, final RefSequence refSequence, boolean onLeft)
    {
        // longer insertions or duplications may be aligned as a soft clipping instead of as a cigar insert
        // to ensure these insertions are captured, searches for candidates in soft-clip by taking the first 12 bases of the ref
        // at the location of the soft clip and then searching in the soft-clip out from the read bounds 5+ bases
        if(onLeft)
        {
            int prevRefPos = readStart - 1;
            int refIndexOffset = prevRefPos - refSequence.Start;

            int refIndexStart = refIndexOffset - SC_INSERT_REF_TEST_LENGTH + 1;
            int refIndexEnd = refIndexStart + SC_INSERT_REF_TEST_LENGTH;

            if(refIndexStart < 0 || refIndexEnd > refSequence.Bases.length) // can occur with SCs at the start of a chromosome
                return null;

            String requiredRefBases = new String(refSequence.Bases, refIndexStart, SC_INSERT_REF_TEST_LENGTH);

            String scBases = readBases.substring(0, scLength);
            int maxMaxIndex = scLength - SC_INSERT_REF_TEST_LENGTH - SC_INSERT_MIN_LENGTH;

            int scMatchIndex = -1;
            int scNextMatchIndex = scBases.indexOf(requiredRefBases, 0);

            while(scNextMatchIndex >= 0)
            {
                // must match at least 5 bases from the end of the soft-clip
                // eg soft-clip = 20, 12 bases of ref then 5 of inserted, so must match 20 - 17 = 3 or earlier
                if(scNextMatchIndex > maxMaxIndex)
                    break;

                scMatchIndex = scNextMatchIndex;
                scNextMatchIndex = scBases.indexOf(requiredRefBases, scNextMatchIndex + 1);
            }

            if(scMatchIndex <= 0)
                return null;

            int impliedVarIndex = scMatchIndex + SC_INSERT_REF_TEST_LENGTH - 1;
            int altLength = scLength - impliedVarIndex - 1;

            String ref = readBases.substring(impliedVarIndex, impliedVarIndex + 1);
            String alt = readBases.substring(impliedVarIndex, impliedVarIndex + altLength + 1);

            return new AltRead(null, ref, alt, 0, false, null, -1);
        }
        else
        {
            int nextRefPos = readEnd + 1;
            int refIndexOffset = nextRefPos - refSequence.Start;

            int refIndexStart = refIndexOffset;
            int refIndexEnd = refIndexStart + SC_INSERT_REF_TEST_LENGTH;

            if(refIndexStart < 0 || refIndexEnd > refSequence.Bases.length) // can occur with SCs at the start of a chromosome
                return null;

            String requiredRefBases = new String(refSequence.Bases, refIndexStart, SC_INSERT_REF_TEST_LENGTH);

            String scBases = readBases.substring(scReadIndex);
            int scMatchIndex = scBases.indexOf(requiredRefBases);

            if(scMatchIndex <= 0 || scMatchIndex < SC_INSERT_MIN_LENGTH)
                return null;

            int impliedVarIndex = scReadIndex - 1; // last aligned/ref base of the read
            int altLength = scMatchIndex;

            String ref = readBases.substring(impliedVarIndex, impliedVarIndex + 1);
            String alt = readBases.substring(impliedVarIndex, impliedVarIndex + altLength + 1);

            return new AltRead(null, ref, alt, 0, false, null, -1);
        }
    }

    private boolean withinReadContext(int readIndex, final SAMRecord record)
    {
        return readIndex >= mConfig.ReadContextFlankLength && readIndex < record.getReadLength() - mConfig.ReadContextFlankLength;
    }

    @VisibleForTesting
    static int mnvLength(int readIndex, int refIndex, byte[] readBases, byte[] refBases)
    {
        final Function<Integer, Boolean> isDifferent =
                i -> refIndex + i < refBases.length && readIndex + i < readBases.length && refBases[refIndex + i] != readBases[readIndex
                        + i];

        if(isDifferent.apply(2))
            return 3;

        return isDifferent.apply((1)) ? 2 : 1;
    }
}
