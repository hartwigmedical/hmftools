package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.SC_INSERT_MIN_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.SC_INSERT_MIN_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.SC_READ_EVENTS_FACTOR;

import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.function.Consumer;
import java.util.function.Function;

import com.beust.jcommander.internal.Sets;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hla.HlaCommon;
import com.hartwig.hmftools.common.samtools.CigarHandler;
import com.hartwig.hmftools.common.samtools.CigarTraversal;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextFactory;
import com.hartwig.hmftools.sage.read.NumberEvents;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class RefContextConsumer implements Consumer<SAMRecord>
{
    private final SageConfig mConfig;
    private final ChrBaseRegion mBounds;
    private final RefSequence mRefGenome;
    private final RefContextCache mRefContextCache;
    private final ReadContextFactory mReadContextFactory;
    private final Set<Integer> mHotspotPositions;

    private int mReadCount;

    public RefContextConsumer(
            final SageConfig config, final ChrBaseRegion bounds, final RefSequence refGenome, final RefContextCache refContextCache,
            final List<VariantHotspot> regionHotspots)
    {
        mBounds = bounds;
        mRefGenome = refGenome;
        mRefContextCache = refContextCache;
        mReadContextFactory = new ReadContextFactory(config.ReadContextFlankSize);
        mConfig = config;

        mHotspotPositions = Sets.newHashSet();
        regionHotspots.forEach(x -> mHotspotPositions.add(x.position()));

        mReadCount = 0;
    }

    public int getReadCount() { return mReadCount; }

    @Override
    public void accept(final SAMRecord record)
    {
        if(!inBounds(record))
            return;

        ++mReadCount;

        if(reachedDepthLimit(record))
            return;

        int numberOfEvents = !record.getSupplementaryAlignmentFlag() ? NumberEvents.calc(record, mRefGenome) : 0;
        int scEvents = (int)NumberEvents.calcSoftClipAdjustment(record);
        int adjustedMapQual = calcAdjustedMapQualLessEventsPenalty(record, numberOfEvents);
        boolean readExceedsQuality = adjustedMapQual > 0;

        if(!readExceedsQuality)
        {
            // if the read is below the required threshold and does not cover a hotspot position, then stop processing it
            if(mHotspotPositions.stream().noneMatch(x -> positionWithin(x, record.getStart(), record.getEnd())))
                return;
        }

        int scAdjustedMapQual = adjustedMapQual - scEvents * mConfig.Quality.MapQualityReadEventsPenalty;
        boolean readExceedsScAdjustedQuality = scAdjustedMapQual > 0;

        final List<AltRead> altReads = Lists.newArrayList();
        final IndexedBases refBases = mRefGenome.alignment();

        final CigarHandler handler = new CigarHandler()
        {
            @Override
            public void handleAlignment(final SAMRecord record, final CigarElement element, final int readIndex, final int refPosition)
            {
                if(record.isSecondaryOrSupplementary())
                    return;

                if(!readExceedsScAdjustedQuality)
                {
                    if(mHotspotPositions.stream().noneMatch(x -> positionWithin(x, record.getStart(), record.getEnd())))
                        return;
                }

                altReads.addAll(processAlignment(
                        record, readIndex, refPosition, element.getLength(), refBases, numberOfEvents, readExceedsScAdjustedQuality));
            }

            @Override
            public void handleInsert(final SAMRecord record, final CigarElement element, final int readIndex, final int refPosition)
            {
                if(record.isSecondaryOrSupplementary())
                    return;

                AltRead altRead = processInsert(
                        element, record, readIndex, refPosition, refBases, numberOfEvents, readExceedsQuality, readExceedsScAdjustedQuality);

                if(altRead != null)
                    altReads.add(altRead);
            }

            @Override
            public void handleDelete(final SAMRecord record, final CigarElement element, final int readIndex, final int refPosition)
            {
                if(record.isSecondaryOrSupplementary())
                    return;

                AltRead altRead = processDel(
                        element, record, readIndex, refPosition, refBases, numberOfEvents, readExceedsQuality, readExceedsScAdjustedQuality);

                if(altRead != null)
                    altReads.add(altRead);
            }

            @Override
            public void handleLeftSoftClip(final SAMRecord record, final CigarElement element)
            {
                if(!mConfig.FindSoftClipInserts)
                    return;

                AltRead altRead = processSoftClip(
                        record, element.getLength(), 0, refBases, readExceedsQuality, numberOfEvents, true);

                if(altRead != null)
                    altReads.add(altRead);
            }

            @Override
            public void handleRightSoftClip(final SAMRecord record, final CigarElement element, int readIndex, int refPosition)
            {
                if(!mConfig.FindSoftClipInserts)
                    return;

                AltRead altRead = processSoftClip(
                        record, element.getLength(), readIndex, refBases, readExceedsQuality, numberOfEvents, false);

                if(altRead != null)
                    altReads.add(altRead);
            }
        };

        CigarTraversal.traverseCigar(record, handler);

        checkCoreExtension(altReads);

        altReads.forEach(AltRead::updateRefContext);
    }

    private int calcAdjustedMapQualLessEventsPenalty(final SAMRecord record, int numberOfEvents)
    {
        // ignore HLA genes
        if(HlaCommon.overlaps(record.getContig(), record.getStart(), record.getEnd()))
            return record.getMappingQuality();

        int eventsPenalty = (numberOfEvents - 1) * mConfig.Quality.MapQualityReadEventsPenalty;

        int improperPenalty = !record.getReadPairedFlag() || !record.getProperPairFlag() || record.getSupplementaryAlignmentFlag() ?
                mConfig.Quality.MapQualityImproperPairPenalty : 0;

        return record.getMappingQuality() - mConfig.Quality.MapQualityFixedPenalty - eventsPenalty - improperPenalty;
    }

    private boolean isHotspotPosition(int position) { return mHotspotPositions.contains(position); }

    private AltRead processInsert(
            final CigarElement element, final SAMRecord record, int readIndex, int refPosition,
            final IndexedBases refBases, int numberOfEvents, boolean readExceedsQuality, boolean readExceedsScAdjustedQuality)
    {
        if(refPosition > mBounds.end() || refPosition < mBounds.start())
            return null;

        boolean exceedsQuality = element.getLength() <= SC_READ_EVENTS_FACTOR ? readExceedsScAdjustedQuality : readExceedsQuality;

        if(!exceedsQuality && !isHotspotPosition(refPosition))
            return null;

        int refIndex = refBases.index(refPosition);
        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        final String ref = new String(refBases.Bases, refIndex, 1);
        final String alt = new String(record.getReadBases(), readIndex, element.getLength() + 1);
        boolean findReadContext = withinReadContext(readIndex, record);

        final RefContext refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);
        if(!reachedDepthLimit(refContext))
        {
            final int baseQuality = baseQuality(readIndex, record, alt.length());
            final ReadContext readContext =
                    findReadContext ? mReadContextFactory.createInsertContext(alt, refPosition, readIndex, record, refBases) : null;
            return new AltRead(refContext, ref, alt, baseQuality, numberOfEvents, sufficientMapQuality, readContext);
        }

        return null;
    }

    private AltRead processDel(
            final CigarElement element, final SAMRecord record, int readIndex, int refPosition,
            final IndexedBases refBases, int numberOfEvents, boolean readExceedsQuality, boolean readExceedsScAdjustedQuality)
    {
        if(refPosition > mBounds.end() || refPosition < mBounds.start())
            return null;

        boolean exceedsQuality = element.getLength() <= SC_READ_EVENTS_FACTOR ? readExceedsScAdjustedQuality : readExceedsQuality;

        if(!exceedsQuality && !isHotspotPosition(refPosition))
            return null;

        int refIndex = refBases.index(refPosition);
        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        final String ref = new String(refBases.Bases, refIndex, element.getLength() + 1);
        final String alt = new String(record.getReadBases(), readIndex, 1);
        boolean findReadContext = withinReadContext(readIndex, record);

        final RefContext refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);
        if(refContext != null && !reachedDepthLimit(refContext))
        {
            final int baseQuality = baseQuality(readIndex, record, 2);
            final ReadContext readContext =
                    findReadContext ? mReadContextFactory.createDelContext(ref, refPosition, readIndex, record, refBases) : null;
            return new AltRead(refContext, ref, alt, baseQuality, numberOfEvents, sufficientMapQuality, readContext);
        }

        return null;
    }

    private List<AltRead> processAlignment(
            final SAMRecord record, int readBasesStartIndex, int refPositionStart,
            int alignmentLength, final IndexedBases refBases, int numberOfEvents, boolean readExceedsQuality)
    {
        final List<AltRead> result = Lists.newArrayList();
        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        int refIndex = refBases.index(refPositionStart);

        for(int i = 0; i < alignmentLength; i++)
        {
            int refPosition = refPositionStart + i;
            int readBaseIndex = readBasesStartIndex + i;
            int refBaseIndex = refIndex + i;

            if(!mBounds.containsPosition(refPosition))
                continue;

            if(!readExceedsQuality && !isHotspotPosition(refPosition))
                continue;

            final byte refByte = refBases.Bases[refBaseIndex];
            final String ref = String.valueOf((char) refByte);
            final byte readByte = record.getReadBases()[readBaseIndex];
            boolean isWithinReadContext = withinReadContext(readBaseIndex, record);

            final RefContext refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);
            if(refContext != null && !reachedDepthLimit(refContext))
            {
                if(readByte != refByte)
                {
                    int baseQuality = record.getBaseQualities()[readBaseIndex];
                    final String alt = String.valueOf((char) readByte);
                    final ReadContext readContext = isWithinReadContext ?
                            mReadContextFactory.createSNVContext(refPosition, readBaseIndex, record, refBases) : null;

                    result.add(new AltRead(refContext, ref, alt, baseQuality, numberOfEvents, sufficientMapQuality, readContext));

                    if(mConfig.MnvEnabled)
                    {
                        int mnvMaxLength = mnvLength(readBaseIndex, refBaseIndex, record.getReadBases(), refBases.Bases);
                        for(int mnvLength = 2; mnvLength <= mnvMaxLength; mnvLength++)
                        {
                            final String mnvRef = new String(refBases.Bases, refBaseIndex, mnvLength);
                            final String mnvAlt = new String(record.getReadBases(), readBaseIndex, mnvLength);

                            // Only check last base because some subsets may not be valid,
                            // ie CA > TA is not a valid subset of CAC > TAT
                            if(mnvRef.charAt(mnvLength - 1) != mnvAlt.charAt(mnvLength - 1))
                            {
                                final ReadContext mnvReadContext = isWithinReadContext ? mReadContextFactory.createMNVContext(refPosition,
                                        readBaseIndex,
                                        mnvLength,
                                        record,
                                        refBases) : null;

                                result.add(new AltRead(refContext,
                                        mnvRef,
                                        mnvAlt,
                                        baseQuality,
                                        NumberEvents.calcWithMnvRaw(numberOfEvents, mnvRef, mnvAlt),
                                        sufficientMapQuality,
                                        mnvReadContext));
                            }
                        }
                    }
                }
                else
                {
                    refContext.refRead(sufficientMapQuality);
                }
            }
        }

        return result;
    }

    private AltRead processSoftClip(
            final SAMRecord record, int scLength, int scReadIndex, final IndexedBases refBases, boolean readExceedsQuality,
            int numberOfEvents, boolean onLeft)
    {
        if(!readExceedsQuality)
            return null;

        if(scLength < SC_INSERT_MIN_FLANK_LENGTH + 1)
            return null;

        AltRead altRead = processSoftClip(
                record.getAlignmentStart(), record.getAlignmentEnd(), record.getReadString(), scLength, scReadIndex, refBases, onLeft);

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

        final RefContext refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);

        if(reachedDepthLimit(refContext))
            return null;

        final int baseQuality = baseQuality(readIndex, record, altRead.Alt.length());

        final ReadContext readContext = mReadContextFactory.createInsertContext(altRead.Alt, refPosition, readIndex, record, refBases);

        SG_LOGGER.trace("soft-clipped insert({}:{} {}>{}) read(index={} {}) softClip(len={} index={} on {})",
                record.getContig(), refPosition, altRead.Ref, altRead.Alt, readIndex, record.getReadName(),
                scLength, scReadIndex, onLeft ? "left" : "right");

        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        AltRead altReadFull = new AltRead(refContext, altRead.Ref, altRead.Alt, baseQuality, numberOfEvents, sufficientMapQuality, readContext);
        return altReadFull;
    }

    public static AltRead processSoftClip(
            int readStart, int readEnd, final String readBases, int scLength, int scReadIndex, final IndexedBases refBases, boolean onLeft)
    {
        // look for an insert of X bases starting at the soft-clip followed by at least 10 bases of matching ref bases in the soft-clipping
        if(onLeft)
        {
            int prevRefPos = readStart - 1;
            int refIndexOffset = prevRefPos - refBases.Position;

            int refIndexStart = refBases.Index + refIndexOffset - SC_INSERT_MIN_FLANK_LENGTH + 1;
            int refIndexEnd = refIndexStart + SC_INSERT_MIN_FLANK_LENGTH;

            if(refIndexStart < 0 || refIndexEnd > refBases.Bases.length) // can occur with SCs at the start of a chromosome
                return null;

            String requiredRefBases = new String(refBases.Bases, refIndexStart, SC_INSERT_MIN_FLANK_LENGTH);

            String scBases = readBases.substring(0, scLength);
            int scMatchIndex = scBases.lastIndexOf(requiredRefBases);

            if(scMatchIndex <= 0)
                return null;

            if(scMatchIndex < SC_INSERT_MIN_LENGTH)
                return null;

            int scIndexMatchEnd = scMatchIndex + SC_INSERT_MIN_FLANK_LENGTH - 1;
            int altLength = scLength - scIndexMatchEnd - 1;

            try
            {
                String ref = readBases.substring(scIndexMatchEnd, scIndexMatchEnd + 1);
                String alt = readBases.substring(scIndexMatchEnd, scIndexMatchEnd + altLength + 1);

                return new AltRead(null, ref, alt, 0, 0, false, null);
            }
            catch(Exception e)
            {
                return null;
            }
        }
        else
        {
            int nextRefPos = readEnd + 1;
            int refIndexOffset = nextRefPos - refBases.Position;

            int refIndexStart = refBases.Index + refIndexOffset;
            int refIndexEnd = refIndexStart + SC_INSERT_MIN_FLANK_LENGTH;

            if(refIndexStart < 0 || refIndexEnd > refBases.Bases.length) // can occur with SCs at the start of a chromosome
                return null;

            String requiredRefBases = new String(refBases.Bases, refIndexStart, SC_INSERT_MIN_FLANK_LENGTH);

            String scBases = readBases.substring(scReadIndex);
            int scMatchIndex = scBases.indexOf(requiredRefBases);

            if(scMatchIndex <= 0)
                return null;

            try
            {
                String ref = readBases.substring(scReadIndex - 1, scReadIndex);
                String alt = readBases.substring(scReadIndex - 1, scReadIndex + scMatchIndex);

                return new AltRead(null, ref, alt, 0, 0, false, null);
            }
            catch(Exception e)
            {
                return null;
            }
        }
    }

    private boolean withinReadContext(int readIndex, final SAMRecord record)
    {
        return readIndex >= mConfig.ReadContextFlankSize && readIndex < record.getReadLength() - mConfig.ReadContextFlankSize;
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

    private boolean inBounds(final SAMRecord record)
    {
        return positionsOverlap(record.getStart(), record.getEnd(), mBounds.start(), mBounds.end());
    }

    private boolean reachedDepthLimit(final RefContext refContext)
    {
        return refContext.exceedsDepthLimit(mConfig.MaxReadDepth, mConfig.MaxReadDepthPanel);
    }

    private boolean reachedDepthLimit(final SAMRecord record)
    {
        RefContext startRefContext = mRefContextCache.getOrCreateRefContext(mBounds.Chromosome, record.getAlignmentStart());

        if(reachedDepthLimit(startRefContext))
            return true;

        RefContext endRefContext = mRefContextCache.getOrCreateRefContext(mBounds.Chromosome, record.getAlignmentEnd());

        return reachedDepthLimit(endRefContext);
    }

    private void checkCoreExtension(final List<AltRead> altReads)
    {
        if(altReads.size() < 2)
            return;

        // if an SNV core overlaps an indel core, then extend the cores of both
        for(int i = 0; i < altReads.size(); i++)
        {
            final AltRead snv = altReads.get(i);
            if(!snv.isIndel() && snv.containsReadContext())
            {
                for(int j = altReads.size() - 1; j > i; j--)
                {
                    final AltRead nextIndel = altReads.get(j);
                    if(nextIndel != null && nextIndel.isIndel() && nextIndel.containsReadContext())
                    {
                        if(nextIndel.leftCoreIndex() - nextIndel.length() <= snv.rightCoreIndex())
                        {
                            snv.extend(nextIndel);
                            nextIndel.extend(snv);
                        }
                    }
                }
            }
        }

        for(int i = altReads.size() - 1; i >= 0; i--)
        {
            final AltRead snv = altReads.get(i);

            if(!snv.isIndel() && snv.containsReadContext())
            {
                for(int j = 0; j < i; j++)
                {
                    final AltRead previousIndel = altReads.get(j);
                    if(previousIndel != null && previousIndel.isIndel() && previousIndel.containsReadContext())
                    {
                        if(previousIndel.rightCoreIndex() + previousIndel.length() >= snv.leftCoreIndex())
                        {
                            previousIndel.extend(snv);
                            snv.extend(previousIndel);
                        }
                    }
                }

            }
        }
    }
}
