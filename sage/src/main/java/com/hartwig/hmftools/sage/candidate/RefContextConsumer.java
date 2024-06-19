package com.hartwig.hmftools.sage.candidate;

import static java.lang.Math.abs;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.bam.CigarUtils.getPositionFromReadIndex;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.SageConstants.MIN_INSERT_ALIGNMENT_OVERLAP;
import static com.hartwig.hmftools.sage.SageConstants.SC_INSERT_REF_TEST_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.SC_INSERT_MIN_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.SC_READ_EVENTS_FACTOR;
import static com.hartwig.hmftools.sage.common.Microhomology.findLeftHomologyShift;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.isImproperPair;

import static htsjdk.samtools.CigarOperator.M;

import java.util.List;
import java.util.Set;
import java.util.function.Function;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.hla.HlaCommon;
import com.hartwig.hmftools.common.bam.CigarHandler;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.common.NumberEvents;
import com.hartwig.hmftools.sage.select.ReadPanelStatus;

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

    private int mReadCount;

    public RefContextConsumer(
            final SageConfig config, final ChrBaseRegion bounds, final RefSequence refSequence, final RefContextCache refContextCache,
            final List<SimpleVariant> regionHotspots)
    {
        mBounds = bounds;
        mRefSequence = refSequence;
        mRefContextCache = refContextCache;
        mReadContextBuilder = new VariantReadContextBuilder(config.ReadContextFlankLength);
        mConfig = config;

        mHotspotPositions = Sets.newHashSet();
        regionHotspots.forEach(x -> mHotspotPositions.add(x.position()));

        mReadCount = 0;
    }

    public int getReadCount() { return mReadCount; }

    public void processRead(final SAMRecord record)
    {
        int readStart = record.getAlignmentStart();
        int readEnd = record.getAlignmentEnd();

        if(!positionsOverlap(readStart, readEnd, mBounds.start(), mBounds.end()))
            return;

        ++mReadCount;

        // check if the read falls within or overlaps a panel region, since this impacts the depth limits
        ReadPanelStatus panelStatus = mRefContextCache.panelSelector().panelStatus(readStart, readEnd);

        if(reachedDepthLimit(readStart, panelStatus) || reachedDepthLimit(readEnd, panelStatus))
            return;

        int numberOfEvents = !record.getSupplementaryAlignmentFlag() ? NumberEvents.calc(record, mRefSequence) : 0;
        int scEvents = (int)NumberEvents.calcSoftClipAdjustment(record);
        int adjustedMapQual = calcAdjustedMapQualLessEventsPenalty(record, numberOfEvents);
        boolean readExceedsQuality = adjustedMapQual > 0;

        final Boolean readCoversHotspot = !readExceedsQuality ?
                mHotspotPositions.stream().anyMatch(x -> positionWithin(x, readStart, readEnd)) : null;

        // if the read is below the required threshold and does not cover a hotspot position, then stop processing it
        if(!readExceedsQuality && !readCoversHotspot)
            return;

        int scAdjustedMapQual = (int)round(adjustedMapQual - scEvents * mConfig.Quality.ReadEventsPenalty);
        boolean readExceedsScAdjustedQuality = scAdjustedMapQual > 0;
        boolean ignoreScAdapter = scEvents > 0 && ignoreSoftClipAdapter(record);

        final List<AltRead> altReads = Lists.newArrayList();

        final CigarHandler handler = new CigarHandler()
        {
            @Override
            public void handleAlignment(
                    final SAMRecord record, final CigarElement element, final int readIndex, final int refPosition)
            {
                if(record.isSecondaryOrSupplementary())
                    return;

                if(!readExceedsScAdjustedQuality)
                {
                    if((readCoversHotspot != null && !readCoversHotspot)
                    || mHotspotPositions.stream().noneMatch(x -> positionWithin(x, readStart, readEnd)))
                    {
                        return;
                    }
                }

                altReads.addAll(processAlignment(
                        record, readIndex, refPosition, element.getLength(), panelStatus, numberOfEvents,
                        readExceedsScAdjustedQuality));
            }

            @Override
            public void handleInsert(final SAMRecord record, final CigarElement element, final int readIndex, final int refPosition)
            {
                if(record.isSecondaryOrSupplementary())
                    return;

                AltRead altRead = processInsert(
                        element, record, readIndex, refPosition, panelStatus, numberOfEvents, readExceedsQuality,
                        readExceedsScAdjustedQuality);

                if(altRead != null)
                    altReads.add(altRead);
            }

            @Override
            public void handleDelete(final SAMRecord record, final CigarElement element, final int readIndex, final int refPosition)
            {
                if(record.isSecondaryOrSupplementary())
                    return;

                AltRead altRead = processDel(
                        element, record, readIndex, refPosition, panelStatus, numberOfEvents, readExceedsQuality,
                        readExceedsScAdjustedQuality);

                if(altRead != null)
                    altReads.add(altRead);
            }

            @Override
            public void handleLeftSoftClip(final SAMRecord record, final CigarElement element)
            {
                if(ignoreScAdapter)
                    return;

                AltRead altRead = processSoftClip(
                        record, element.getLength(), 0, panelStatus, mRefSequence, readExceedsQuality, numberOfEvents, true);

                if(altRead != null)
                    altReads.add(altRead);
            }

            @Override
            public void handleRightSoftClip(final SAMRecord record, final CigarElement element, int readIndex, int refPosition)
            {
                if(ignoreScAdapter)
                    return;

                AltRead altRead = processSoftClip(
                        record, element.getLength(), readIndex, panelStatus, mRefSequence, readExceedsQuality, numberOfEvents, false);

                if(altRead != null)
                    altReads.add(altRead);
            }
        };

        CigarHandler.traverseCigar(record, handler);

        for(AltRead altRead : altReads)
        {
            altRead.updateRefContext();

            if(altRead.SufficientMapQuality)
                mRefContextCache.incrementDepth(altRead.position());
        }
    }

    public static boolean ignoreSoftClipAdapter(final SAMRecord record)
    {
        // ignore soft-clips from short fragments indicating adapter bases are present
        int fragmentLength = abs(record.getInferredInsertSize());

        if(fragmentLength >= record.getReadBases().length + MIN_INSERT_ALIGNMENT_OVERLAP)
            return false;

        int alignedBases = record.getCigar().getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();
        int insertAlignmentOverlap = abs(fragmentLength - alignedBases);
        return insertAlignmentOverlap < MIN_INSERT_ALIGNMENT_OVERLAP;
    }

    private boolean reachedDepthLimit(int position, final ReadPanelStatus panelStatus)
    {
        Boolean exceedsLimit = mRefContextCache.exceedsDepthLimit(position);

        if(exceedsLimit != null)
            return exceedsLimit;

        // set depth limit on the first time this position is processed
        int depthLimit = depthLimit(panelStatus, position);
        mRefContextCache.registerDepthLimit(position, depthLimit);
        return false;
    }

    private int depthLimit(final ReadPanelStatus panelStatus, int position)
    {
        if(mConfig.IncludeMT && MitochondrialChromosome.contains(mBounds.Chromosome))
            return mConfig.MaxReadDepthPanel;

        if(panelStatus == ReadPanelStatus.WITHIN_PANEL)
            return mConfig.MaxReadDepthPanel;
        else if(panelStatus == ReadPanelStatus.OUTSIDE_PANEL)
            return mConfig.MaxReadDepth;

        return mRefContextCache.panelSelector().panelStatus(position) == ReadPanelStatus.WITHIN_PANEL ?
                mConfig.MaxReadDepthPanel : mConfig.MaxReadDepth;
    }

    private int calcAdjustedMapQualLessEventsPenalty(final SAMRecord record, int numberOfEvents)
    {
        // ignore HLA genes
        if(HlaCommon.overlaps(record.getContig(), record.getStart(), record.getEnd()))
            return record.getMappingQuality();

        int eventPenalty = (int)round((numberOfEvents - 1) * mConfig.Quality.ReadEventsPenalty);

        int improperPenalty = isImproperPair(record) || record.getSupplementaryAlignmentFlag() ?
                mConfig.Quality.ImproperPairPenalty : 0;

        return record.getMappingQuality() - mConfig.Quality.FixedPenalty - eventPenalty - improperPenalty;
    }

    private boolean isHotspotPosition(int position) { return mHotspotPositions.contains(position); }

    private AltRead processInsert(
            final CigarElement element, final SAMRecord record, int readIndex, int refPosition, final ReadPanelStatus panelStatus,
            int numberOfEvents, boolean readExceedsQuality, boolean readExceedsScAdjustedQuality)
    {
        if(!mBounds.containsPosition(refPosition))
            return null;

        boolean exceedsQuality = element.getLength() <= SC_READ_EVENTS_FACTOR ? readExceedsScAdjustedQuality : readExceedsQuality;

        if(!exceedsQuality && !isHotspotPosition(refPosition))
            return null;

        if(reachedDepthLimit(refPosition, panelStatus))
            return null;

        int refIndex = mRefSequence.index(refPosition);
        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        String ref = new String(mRefSequence.Bases, refIndex, 1);
        String alt = new String(record.getReadBases(), readIndex, element.getLength() + 1);

        RefContext refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);

        int baseQuality = baseQuality(readIndex, record, alt.length());

        SimpleVariant variant = new SimpleVariant(record.getContig(), refPosition, ref, alt);

        VariantReadContext readContext = mReadContextBuilder.createContext(variant, record, readIndex, mRefSequence);

        if(readContext == null)
            return null;

        return new AltRead(refContext, ref, alt, baseQuality, numberOfEvents, sufficientMapQuality, readContext);
    }

    private AltRead processDel(
            final CigarElement element, final SAMRecord record, int readIndex, int refPosition, final ReadPanelStatus panelStatus,
            int numberOfEvents, boolean readExceedsQuality, boolean readExceedsScAdjustedQuality)
    {
        if(!mBounds.containsPosition(refPosition))
            return null;

        boolean exceedsQuality = element.getLength() <= SC_READ_EVENTS_FACTOR ? readExceedsScAdjustedQuality : readExceedsQuality;

        if(!exceedsQuality && !isHotspotPosition(refPosition))
            return null;

        if(reachedDepthLimit(refPosition, panelStatus))
            return null;

        int refIndex = mRefSequence.index(refPosition);
        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        final String ref = new String(mRefSequence.Bases, refIndex, element.getLength() + 1);
        final String alt = new String(record.getReadBases(), readIndex, 1);

        final RefContext refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);
        if(refContext != null)
        {
            int baseQuality = baseQuality(readIndex, record, 2);

            SimpleVariant variant = new SimpleVariant(record.getContig(), refPosition, ref, alt);

            VariantReadContext readContext = mReadContextBuilder.createContext(variant, record, readIndex, mRefSequence);

            if(readContext == null)
                return null;

            return new AltRead(refContext, ref, alt, baseQuality, numberOfEvents, sufficientMapQuality, readContext);
        }

        return null;
    }

    private List<AltRead> processAlignment(
            final SAMRecord record, int readBasesStartIndex, int refPositionStart, int alignmentLength,
            final ReadPanelStatus panelStatus, int numberOfEvents, boolean readExceedsQuality)
    {
        List<AltRead> result = Lists.newArrayList();
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

            if(!readExceedsQuality && !isHotspotPosition(refPosition))
                continue;

            if(reachedDepthLimit(refPosition, panelStatus))
                continue;

            final byte refByte = mRefSequence.Bases[refBaseIndex];
            final String ref = String.valueOf((char) refByte);
            final byte readByte = record.getReadBases()[readBaseIndex];

            if(readByte != refByte)
            {
                final RefContext refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);
                if(refContext == null)
                    continue;

                int baseQuality = record.getBaseQualities()[readBaseIndex];
                final String alt = String.valueOf((char) readByte);

                SimpleVariant variant = new SimpleVariant(record.getContig(), refPosition, ref, alt);
                VariantReadContext readContext = mReadContextBuilder.createContext(variant, record, readBaseIndex, mRefSequence);

                if(readContext != null)
                    result.add(new AltRead(refContext, ref, alt, baseQuality, numberOfEvents, sufficientMapQuality, readContext));

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
                        SimpleVariant mnv = new SimpleVariant(record.getContig(), refPosition, mnvRef, mnvAlt);

                        VariantReadContext mnvReadContext = mReadContextBuilder.createContext(mnv, record, readBaseIndex, mRefSequence);

                        if(readContext != null)
                        {
                            result.add(new AltRead(
                                    refContext, mnvRef, mnvAlt, baseQuality, NumberEvents.calcWithMnvRaw(numberOfEvents, mnvRef, mnvAlt),
                                    sufficientMapQuality, mnvReadContext));
                        }
                    }
                }
            }
            else
            {
                if(sufficientMapQuality)
                    mRefContextCache.incrementDepth(refPosition);
            }
        }

        return result;
    }

    private AltRead processSoftClip(
            final SAMRecord record, int scLength, int scReadIndex, final ReadPanelStatus panelStatus, final RefSequence refSequence,
            boolean readExceedsQuality, int numberOfEvents, boolean onLeft)
    {
        if(!readExceedsQuality)
            return null;

        if(scLength < SC_INSERT_REF_TEST_LENGTH + 1)
            return null;

        AltRead altRead = processSoftClip(
                record.getAlignmentStart(), record.getAlignmentEnd(), record.getReadString(), scLength, scReadIndex, refSequence, onLeft);

        if(altRead == null)
            return null;

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
            return null;

        if(!withinReadContext(readIndex, record))
            return null;

        if(reachedDepthLimit(refPosition, panelStatus))
            return null;

        RefContext refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);

        int baseQuality = baseQuality(readIndex, record, altRead.Alt.length());

        SimpleVariant variant = new SimpleVariant(record.getContig(), refPosition, altRead.Ref, altRead.Alt);

        if(variant.isIndel())
        {
            int leftHomologyShift = findLeftHomologyShift(variant, mRefSequence, record.getReadBases(), readIndex);

            if(leftHomologyShift > 0)
            {
                readIndex -= leftHomologyShift;

                // recompute the new reference position, taking into consideration indels in the read
                refPosition = getPositionFromReadIndex(record.getAlignmentStart(), record.getCigar().getCigarElements(), readIndex);

                if(refPosition != NO_POSITION)
                {
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
        }

        VariantReadContext readContext = mReadContextBuilder.createContext(variant, record, readIndex, mRefSequence);

        if(readContext == null)
            return null;

        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        AltRead altReadFull = new AltRead(refContext, altRead.Ref, altRead.Alt, baseQuality, numberOfEvents, sufficientMapQuality, readContext);
        return altReadFull;
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
            int scMatchIndex = scBases.lastIndexOf(requiredRefBases);

            if(scMatchIndex <= 0)
                return null;

            // must match at least 5 bases from the end of the soft-clip
            // eg soft-clip = 20, 12 bases of ref then 5 of inserted, so must match 20 - 17 = 3 or earlier
            int maxMaxIndex = scLength - SC_INSERT_REF_TEST_LENGTH - SC_INSERT_MIN_LENGTH;
            if(scMatchIndex > maxMaxIndex)
                return null;

            int impliedVarIndex = scMatchIndex + SC_INSERT_REF_TEST_LENGTH - 1;
            int altLength = scLength - impliedVarIndex - 1;

            String ref = readBases.substring(impliedVarIndex, impliedVarIndex + 1);
            String alt = readBases.substring(impliedVarIndex, impliedVarIndex + altLength + 1);

            return new AltRead(null, ref, alt, 0, 0, false, null);
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

            return new AltRead(null, ref, alt, 0, 0, false, null);
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

    private int baseQuality(int readIndex, final SAMRecord record, int length)
    {
        int maxIndex = Math.min(readIndex + length, record.getBaseQualities().length) - 1;
        int quality = Integer.MAX_VALUE;
        for(int i = readIndex; i <= maxIndex; i++)
        {
            quality = Math.min(quality, record.getBaseQualities()[i]);
        }
        return quality;
    }
}
