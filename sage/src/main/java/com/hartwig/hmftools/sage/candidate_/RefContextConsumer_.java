package com.hartwig.hmftools.sage.candidate_;

import static java.lang.Math.abs;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.SageConstants.MIN_INSERT_ALIGNMENT_OVERLAP;
import static com.hartwig.hmftools.sage.SageConstants.SC_INSERT_MIN_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.SC_INSERT_MIN_SC_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.SC_READ_EVENTS_FACTOR;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.isImproperPair;

import static htsjdk.samtools.CigarOperator.M;

import java.util.List;
import java.util.Set;
import java.util.function.Function;

import com.beust.jcommander.internal.Sets;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.hla.HlaCommon;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.CigarHandler;
import com.hartwig.hmftools.common.samtools.CigarTraversal;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.read.NumberEvents;
import com.hartwig.hmftools.sage.select.ReadPanelStatus;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

/**
 * Processes reads and generates AltReads that are stored in a RefContextCache.
 */
public class RefContextConsumer_
{
    private final SageConfig mConfig;
    private final ChrBaseRegion mBounds;
    private final RefSequence mRefGenome;
    private final RefContextCache_ mRefContextCache;
    private final ReadContextFactory_ mReadContextFactory;
    private final Set<Integer> mHotspotPositions;

    private int mReadCount;

    public RefContextConsumer_(
            final SageConfig config, final ChrBaseRegion bounds, final RefSequence refGenome, final RefContextCache_ refContextCache,
            final List<VariantHotspot> regionHotspots)
    {
        mBounds = bounds;
        mRefGenome = refGenome;
        mRefContextCache = refContextCache;
        mReadContextFactory = new ReadContextFactory_(config.ReadContextFlankSize);
        mConfig = config;

        mHotspotPositions = Sets.newHashSet();
        regionHotspots.forEach(x -> mHotspotPositions.add(x.position()));

        mReadCount = 0;
    }

    /**
     * Has the depth limit in the ref context cache has been exceeded. If depth limit has not been set then set it to the panel status at
     * that position and return false.
     */
    private boolean reachedDepthLimit(int position, final ReadPanelStatus panelStatus)
    {
        Boolean exceedsLimit = mRefContextCache.exceedsDepthLimit(position);

        if(exceedsLimit != null)
        {
            return exceedsLimit;
        }

        // set depth limit on the first time this position is processed
        int depthLimit = depthLimit(panelStatus, position);
        mRefContextCache.registerDepthLimit(position, depthLimit);
        return false;
    }

    /**
     * Returns the max depth limit based on the panel status at position.
     */
    private int depthLimit(final ReadPanelStatus panelStatus, int position)
    {
        if(mConfig.IncludeMT && MitochondrialChromosome.contains(mBounds.Chromosome))
        {
            return mConfig.MaxReadDepthPanel;
        }

        if(panelStatus == ReadPanelStatus.WITHIN_PANEL)
        {
            return mConfig.MaxReadDepthPanel;
        }
        else if(panelStatus == ReadPanelStatus.OUTSIDE_PANEL)
        {
            return mConfig.MaxReadDepth;
        }

        return mRefContextCache.panelSelector().panelStatus(position) == ReadPanelStatus.WITHIN_PANEL ?
                mConfig.MaxReadDepthPanel : mConfig.MaxReadDepth;
    }

    /**
     * Returns the MapQ penalised by events (downshifted by 1) and if the record is part of an improper pair, we also subtracts a fixed
     * penalty.
     * <p>
     * NOTE: If it overlaps an HLA gene then we do not penalise.
     */
    private int calcAdjustedMapQualLessEventsPenalty(final SAMRecord record, int numberOfEvents)
    {
        // ignore HLA genes
        if(HlaCommon.overlaps(record.getContig(), record.getStart(), record.getEnd()))
        {
            return record.getMappingQuality();
        }

        int eventPenalty = (int) round((numberOfEvents - 1) * mConfig.Quality.ReadEventsPenalty);

        int improperPenalty = isImproperPair(record) || record.getSupplementaryAlignmentFlag() ?
                mConfig.Quality.ImproperPairPenalty : 0;

        return record.getMappingQuality() - mConfig.Quality.FixedPenalty - eventPenalty - improperPenalty;
    }

    /**
     * Processes a read generating AltReads in the RefContextCache.
     */
    public void processRead(final SAMRecord record)
    {
        int readStart = record.getAlignmentStart();
        int readEnd = record.getAlignmentEnd();

        // Filter out if the read does no overlap the bounds of the current region.
        if(!positionsOverlap(readStart, readEnd, mBounds.start(), mBounds.end()))
            return;

        ++mReadCount;

        // check if the read falls within or overlaps a panel region, since this impacts the depth limits
        ReadPanelStatus panelStatus = mRefContextCache.panelSelector().panelStatus(readStart, readEnd);

        // Filter out if we reached the depth limit at the start or end of read bases on whether we are in a panel region.
        if(reachedDepthLimit(readStart, panelStatus) || reachedDepthLimit(readEnd, panelStatus))
            return;

        // Compute number of events. Note that we set this to zero for supp reads. Use this to compute adjusted MapQ, and set a flag
        // whether the read exceeds the quality, i.e. it has a positive adjusted MapQ.
        int numberOfEvents = !record.getSupplementaryAlignmentFlag() ? NumberEvents.calc(record, mRefGenome) : 0;
        int scEvents = (int)NumberEvents.calcSoftClipAdjustment(record);
        int adjustedMapQual = calcAdjustedMapQualLessEventsPenalty(record, numberOfEvents);
        boolean readExceedsQuality = adjustedMapQual > 0;

        // If the read doesn't exceed quality, does it cover a hotspot?
        final Boolean readCoversHotspot = !readExceedsQuality ?
                mHotspotPositions.stream().anyMatch(x -> positionWithin(x, readStart, readEnd)) : null;

        // If the read is below the required quality threshold and does not cover a hotspot position, then stop processing it.
        if(!readExceedsQuality && !readCoversHotspot)
            return;

        // Calculate softlcip adjusted MapQ, and whether we exceed the required quality theshold.
        int scAdjustedMapQual = (int)round(adjustedMapQual - scEvents * mConfig.Quality.ReadEventsPenalty);
        boolean readExceedsScAdjustedQuality = scAdjustedMapQual > 0;

        // Should we ignore the softclip adaptor?
        boolean ignoreScAdapter = scEvents > 0 && ignoreSoftClipAdapter(record);

        // Process the cigar and collect AltReads.
        // These AltReads will be sorted by position.
        final List<AltRead_> altReads = Lists.newArrayList();
        final IndexedBases_ refBases = mRefGenome.alignment();

        final CigarHandler handler = new CigarHandler()
        {
            @Override
            public void handleAlignment(
                    final SAMRecord record, final CigarElement element, boolean beforeIndel, final int readIndex, final int refPosition)
            {
                // Filter out secondary or supp reads.
                if(record.isSecondaryOrSupplementary())
                    return;

                // Filter if we do not exceed softclip adjusted quality and we are not in a hotspot.
                if(!readExceedsScAdjustedQuality)
                {
                    if((readCoversHotspot != null && !readCoversHotspot)
                            || mHotspotPositions.stream().noneMatch(x -> positionWithin(x, readStart, readEnd)))
                    {
                        return;
                    }
                }

                // Collect all AltReads.
                altReads.addAll(processAlignment(
                        record, readIndex, refPosition, element.getLength(), panelStatus, refBases, numberOfEvents,
                        readExceedsScAdjustedQuality));
            }

            @Override
            public void handleInsert(final SAMRecord record, final CigarElement element, final int readIndex, final int refPosition)
            {
                if(record.isSecondaryOrSupplementary())
                    return;

                AltRead_ altRead = processInsert(
                        element, record, readIndex, refPosition, panelStatus, refBases, numberOfEvents, readExceedsQuality,
                        readExceedsScAdjustedQuality);

                if(altRead != null)
                    altReads.add(altRead);
            }

            @Override
            public void handleDelete(final SAMRecord record, final CigarElement element, final int readIndex, final int refPosition)
            {
                if(record.isSecondaryOrSupplementary())
                    return;

                AltRead_ altRead = processDel(
                        element, record, readIndex, refPosition, panelStatus, refBases, numberOfEvents, readExceedsQuality,
                        readExceedsScAdjustedQuality);

                if(altRead != null)
                    altReads.add(altRead);
            }

            @Override
            public void handleLeftSoftClip(final SAMRecord record, final CigarElement element)
            {
                // Filter if we are ignore the SC adaptor.
                // Perhaps this is the only way a soft clip can disguise an insert.
                if(ignoreScAdapter)
                    return;

                AltRead_ altRead = processSoftClip(
                        record, element.getLength(), 0, panelStatus, refBases, readExceedsQuality, numberOfEvents, true);

                if(altRead != null)
                    altReads.add(altRead);
            }

            @Override
            public void handleRightSoftClip(final SAMRecord record, final CigarElement element, int readIndex, int refPosition)
            {
                // Filter if we are ignore the SC adaptor.
                // Perhaps this is the only way a soft clip can disguise an insert.
                if(ignoreScAdapter)
                    return;

                AltRead_ altRead = processSoftClip(
                        record, element.getLength(), readIndex, panelStatus, refBases, readExceedsQuality, numberOfEvents, false);

                if(altRead != null)
                    altReads.add(altRead);
            }
        };

        CigarTraversal.traverseCigar(record, handler);

        // Extend AltReads of overlapping SNV and Indel cores.
        checkCoreExtension(altReads);

        // Update the associated ref contexts for each alt read generated.
        // If it has sufficient MapQ then we increment the depth at its position.
        for(AltRead_ altRead : altReads)
        {
            altRead.updateRefContext();

            if(altRead.SufficientMapQuality)
                mRefContextCache.incrementDepth(altRead.position());
        }
    }

    /**
     * Returns whether it is possible for the soft clip can disguise an insert.
     */
    public static boolean ignoreSoftClipAdapter(final SAMRecord record)
    {
        // Don't really understand how this does what it does.
        int fragmentLength = abs(record.getInferredInsertSize());
        if(fragmentLength >= record.getReadBases().length + MIN_INSERT_ALIGNMENT_OVERLAP)
            return false;

        int alignedBases = record.getCigar().getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();
        int insertAlignmentOverlap = abs(fragmentLength - alignedBases);
        return insertAlignmentOverlap < MIN_INSERT_ALIGNMENT_OVERLAP;
    }

    /**
     * Extend AltReads of overlapping SNV and Indel cores.
     */
    private void checkCoreExtension(final List<AltRead_> altReads)
    {
        // Note: That all AltReads come from the same read, and these AltReads are sorted by position.

        // We have to at least two alt reads.
        if(altReads.size() < 2)
        {
            return;
        }

        // If an SNV core overlaps an indel core, then extend the cores of both. Why?
        for(int i = 0; i < altReads.size(); i++)
        {
            final AltRead_ snv = altReads.get(i);
            if(!snv.isIndel() && snv.containsReadContext())
            {
                for(int j = altReads.size() - 1; j > i; j--)
                {
                    final AltRead_ nextIndel = altReads.get(j);
                    if(nextIndel != null && nextIndel.isIndel() && nextIndel.containsReadContext())
                    {
                        // Remember that the AltReads are sorted by position.
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
            final AltRead_ snv = altReads.get(i);

            if(!snv.isIndel() && snv.containsReadContext())
            {
                for(int j = 0; j < i; j++)
                {
                    final AltRead_ previousIndel = altReads.get(j);
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

    /**
     * Process this M region and return AltReads based on different bases.
     */
    private List<AltRead_> processAlignment(
            final SAMRecord record, int readBasesStartIndex, int refPositionStart, int alignmentLength,
            final ReadPanelStatus panelStatus, final IndexedBases_ refBases, int numberOfEvents, boolean readExceedsQuality)
    {
        List<AltRead_> result = Lists.newArrayList();

        // Is it of sufficient quality?
        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        int refIndex = refBases.index(refPositionStart);

        // Go through each base in the alignment.
        for(int i = 0; i < alignmentLength; i++)
        {
            int refPosition = refPositionStart + i;
            int readBaseIndex = readBasesStartIndex + i;
            int refBaseIndex = refIndex + i;

            // Filter out positions that are not within the bounds of the region.
            if(refPosition < mBounds.start())
                continue;

            if(refPosition > mBounds.end())
                break;

            // Filter out position if it isn't in a hotspot and the read does not exceed quality.
            if(!readExceedsQuality && !isHotspotPosition(refPosition))
                continue;

            // Filter out if we have reached the depth limit.
            if(reachedDepthLimit(refPosition, panelStatus))
                continue;

            final byte refByte = refBases.Bases[refBaseIndex];
            final String ref = String.valueOf((char) refByte);
            final byte readByte = record.getReadBases()[readBaseIndex];
            boolean isNotWithinFlanks = notWithinFlanks(readBaseIndex, record);

            // We have a mismatch.
            if(readByte != refByte)
            {
                // Get ref context for this GenomePosition.
                final RefContext_ refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);

                // Has the position been evicted from the ref context cache?
                if(refContext == null)
                    continue;

                int baseQuality = record.getBaseQualities()[readBaseIndex];
                final String alt = String.valueOf((char) readByte);

                // Create an SNV read context if the base is not too close to the end of the read.
                final ReadContext_ readContext = isNotWithinFlanks ?
                        mReadContextFactory.createSNVContext(refPosition, readBaseIndex, record, refBases) : null;

                // Add a new AltRead to the results.
                result.add(new AltRead_(refContext, ref, alt, baseQuality, numberOfEvents, sufficientMapQuality, readContext));

                int mnvMaxLength = mnvLength(readBaseIndex, refBaseIndex, record.getReadBases(), refBases.Bases);

                // If MVN, i.e. mnvMaxLength > 1, then loop through each base in the mnv.
                int nextReadIndex = i;
                for(int mnvLength = 2; mnvLength <= mnvMaxLength; mnvLength++)
                {
                    ++nextReadIndex;

                    // MNVs cannot extend past the end of this Cigar element
                    if(nextReadIndex >= alignmentLength)
                        break;

                    // Get the ref and alt bases at this MVN truncated to mnvLength.
                    final String mnvRef = new String(refBases.Bases, refBaseIndex, mnvLength);
                    final String mnvAlt = new String(record.getReadBases(), readBaseIndex, mnvLength);

                    // Check whether the alt and ref bases of the MVN are different.
                    // Only check last base because some subsets may not be valid,
                    // ie CA > TA is not a valid subset of CAC > TAT
                    if(mnvRef.charAt(mnvLength - 1) != mnvAlt.charAt(mnvLength - 1))
                    {
                        // Create and add a new Alt read to the result if the base is not too close to the end of the read.
                        // This is based of a MVN context.
                        ReadContext_ mnvReadContext = isNotWithinFlanks ? mReadContextFactory.createMNVContext(refPosition,
                                readBaseIndex,
                                mnvLength,
                                record.getReadBases(),
                                refBases) : null;

                        result.add(new AltRead_(refContext,
                                mnvRef,
                                mnvAlt,
                                baseQuality,
                                NumberEvents.calcWithMnvRaw(numberOfEvents, mnvRef, mnvAlt),
                                sufficientMapQuality,
                                mnvReadContext));
                    }
                }
            }
            else
            {
                // If MapQ is sufficient then increment the depth at the refPosition.
                if(sufficientMapQuality)
                    mRefContextCache.incrementDepth(refPosition);
            }
        }

        return result;
    }

    /**
     * Returns the length of the MVN (max MVN length is three) assuming that the read base is different than refBase.
     */
    @VisibleForTesting
    public static int mnvLength(int readIndex, int refIndex, byte[] readBases, byte[] refBases)
    {
        // Returns true when the refBase and readBase is valid at offset i, and they are different values, otherwise returns false.
        final Function<Integer, Boolean> isDifferent = (final Integer i) ->
        {
            if(refIndex + i >= refBases.length)
            {
                return false;
            }

            if(readIndex + i >= readBases.length)
            {
                return false;
            }

            return refBases[refIndex + i] != readBases[readIndex + i];
        };

        // Two bases after a different, then MVN has length 3. I guess this is the max MVN length.
        if(isDifferent.apply(2))
        {
            return 3;
        }

        return isDifferent.apply(1) ? 2 : 1;
    }

    /**
     * Returns the AltRead associated to the Del ReadContext. There are various filters, if one of these is triggered then returns null.
     */
    @Nullable
    private AltRead_ processDel(
            final CigarElement element, final SAMRecord record, int readIndex, int refPosition, final ReadPanelStatus panelStatus,
            final IndexedBases_ refBases, int numberOfEvents, boolean readExceedsQuality, boolean readExceedsScAdjustedQuality)
    {
        // Are we out of bounds?
        if(!mBounds.containsPosition(refPosition))
        {
            return null;
        }

        // Do we exceed quality based on whether it is a short del or a long del?
        boolean exceedsQuality = element.getLength() <= SC_READ_EVENTS_FACTOR ? readExceedsScAdjustedQuality : readExceedsQuality;

        // Return null if we are not in a hotpot and the quality is not exceeded.
        if(!exceedsQuality && !isHotspotPosition(refPosition))
        {
            return null;
        }

        // Return null if we have reached the depth limit.
        if(reachedDepthLimit(refPosition, panelStatus))
        {
            return null;
        }

        int refIndex = refBases.index(refPosition);

        // Do we have sufficient MapQ?
        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        // Get the ref covering the del plus one base before, and the alt base before the del.
        final String ref = new String(refBases.Bases, refIndex, element.getLength() + 1);
        final String alt = new String(record.getReadBases(), readIndex, 1);
        boolean isNotWithinFlanks = notWithinFlanks(readIndex, record);

        // Get the ref context from the ref context cache at this GenomePosition.
        final RefContext_ refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);
        if(refContext != null)
        {
            // min base quality from readIndex to readIndex + 2.
            final int minBaseQuality = minBaseQuality(readIndex, record, 2);

            // Create ReadContext (DelContext) if the readIndex is not too close to the ends of the read.
            final ReadContext_ readContext = isNotWithinFlanks ?
                    mReadContextFactory.createDelContext(ref, refPosition, readIndex, record.getReadBases(), refBases) : null;

            // Return the AltRead wrapping the ReadContext.
            return new AltRead_(refContext, ref, alt, minBaseQuality, numberOfEvents, sufficientMapQuality, readContext);
        }

        return null;
    }

    /**
     * Returns the AltRead associated to the Insert ReadContext. There are various filters, if one of these is triggered then returns null.
     */
    private AltRead_ processInsert(
            final CigarElement element, final SAMRecord record, int readIndex, int refPosition, final ReadPanelStatus panelStatus,
            final IndexedBases_ refBases, int numberOfEvents, boolean readExceedsQuality, boolean readExceedsScAdjustedQuality)
    {
        // Are we out of bounds?
        if(!mBounds.containsPosition(refPosition))
        {
            return null;
        }

        // Do we exceed quality based on whether it is a short ins or a long ins?
        boolean exceedsQuality = element.getLength() <= SC_READ_EVENTS_FACTOR ? readExceedsScAdjustedQuality : readExceedsQuality;

        // Return null if we are not in a hotpot and the quality is not exceeded.
        if(!exceedsQuality && !isHotspotPosition(refPosition))
        {
            return null;
        }

        // Return null if the depth limit has been reached.
        if(reachedDepthLimit(refPosition, panelStatus))
        {
            return null;
        }

        int refIndex = refBases.index(refPosition);

        // Is there sufficient MapQ?
        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        // One ref base before the ins.
        final String ref = new String(refBases.Bases, refIndex, 1);

        // One read base before the ins + the ins.
        final String alt = new String(record.getReadBases(), readIndex, element.getLength() + 1);

        boolean isNotWithinFlanks = notWithinFlanks(readIndex, record);

        // Get ref context from the ref context cache based on this GenomePosition.
        RefContext_ refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);

        int minBaseQuality = minBaseQuality(readIndex, record, alt.length());

        // Create ReadContext (InsertContext) if the readIndex is not too close to the ends of the read.
        ReadContext_ readContext = isNotWithinFlanks ?
                mReadContextFactory.createInsertContext(alt, refPosition, readIndex, record.getReadBases(), refBases) : null;

        // Return the AltRead wrapping the ReadContext.
        return new AltRead_(refContext, ref, alt, minBaseQuality, numberOfEvents, sufficientMapQuality, readContext);
    }

    /**
     * Returns the AltRead associated to an Insert ReadContext, if the soft clip disguises an insert, otherwise returns null.
     * Moreover, there are various filters, if one of these is triggered then returns null.
     */
    private AltRead_ processSoftClip(
            final SAMRecord record, int scLength, int scReadIndex, final ReadPanelStatus panelStatus, final IndexedBases_ refBases,
            boolean readExceedsQuality, int numberOfEvents, boolean onLeft)
    {
        // If read does not exceed quality then return null.
        if(!readExceedsQuality)
            return null;

        // If soft-clip is too short then we return null.
        if(scLength < SC_INSERT_MIN_SC_LENGTH + 1)
            return null;

        // Is this soft clip disguising an insert? If not, then just return null.
        AltRead_ altRead = isSoftClipInsert(
                record.getAlignmentStart(), record.getAlignmentEnd(), record.getReadString(), scLength, scReadIndex, refBases, onLeft);
        if(altRead == null)
            return null;

        int refPosition;
        int readIndex;

        // Not exactly sure why these values are set as they are?
        if(onLeft)
        {
            // Last index of the left soft clip.
            refPosition = record.getAlignmentStart() - 1;

            // Set to start at the ref/alt base prior to the insert.
            readIndex = scLength - altRead.Alt.length();
        }
        else
        {
            // One index before right soft clip.
            refPosition = record.getAlignmentEnd();

            // Index in the read of the base before the right soft clip.
            readIndex = record.getReadBases().length - scLength - 1;
        }

        // Filter out if refPosition is not in bounds.
        if(!mBounds.containsPosition(refPosition))
            return null;

        // Filter out if the readIndex is too close to the start/end of the read.
        if(!notWithinFlanks(readIndex, record))
            return null;

        // Filter out if depth limit is reached.
        if(reachedDepthLimit(refPosition, panelStatus))
            return null;

        // Get RefContext at GenomePosition from cache.
        RefContext_ refContext = mRefContextCache.getOrCreateRefContext(record.getContig(), refPosition);

        // Get min base quality of record at [readIndex, readIndex + altRead.Alt.length()]
        int minBaseQuality = minBaseQuality(readIndex, record, altRead.Alt.length());

        // Create read context, and wrap the old AltRead_.
        ReadContext_ readContext = mReadContextFactory.createInsertContext(altRead.Alt, refPosition, readIndex, record.getReadBases(), refBases);

        boolean sufficientMapQuality = record.getMappingQuality() >= mConfig.MinMapQuality;

        AltRead_ altReadFull =
                new AltRead_(refContext, altRead.Ref, altRead.Alt, minBaseQuality, numberOfEvents, sufficientMapQuality, readContext);
        return altReadFull;
    }

    /**
     * Checks to see if a soft clip is disguising an insert.
     *
     * @return An AltRead_ associated to the insert if the soft clip is disguising an insert, null otherwise.
     */
    public static AltRead_ isSoftClipInsert(
            int readStart, int readEnd, final String readBases, int scLength, int scReadIndex, final IndexedBases_ refBases, boolean onLeft)
    {
        if(onLeft)
        {
            // Check whether SC_INSERT_MIN_SC_LENGTH of ref bases before the start of the read alignment match in the soft clip.
            // Note that this won't match the end of the soft clip, because if it did then it would have been aligned.
            // We do not consider it a match if it matches too close to the start of the soft clip. Why?
            // If this does happen we can consider the gap between the matched ref bases and the start of the read alignment as a insert.
            int refIndexAtReadStart = refBases.Index + readStart - refBases.Position;

            // We look back SC_INSERT_MIN_SC_LENGTH bases from the start of the read, and then to the ref base of the first base of the read.
            int refIndexStart = refIndexAtReadStart - SC_INSERT_MIN_SC_LENGTH;

            // Does the soft clip go over the start/end of a chromosome.
            if(refIndexStart < 0 || refIndexAtReadStart > refBases.Bases.length)
                return null;

            // The SC_INSERT_MIN_SC_LENGTH ref bases just before the start of the read.
            String requiredRefBases = new String(refBases.Bases, refIndexStart, SC_INSERT_MIN_SC_LENGTH);

            // Get the left soft clipped bases of the read.
            String scBases = readBases.substring(0, scLength);

            // Last index of requiredRefBases as a substring of scBases.
            int scMatchIndexStart = scBases.lastIndexOf(requiredRefBases);
            if(scMatchIndexStart < SC_INSERT_MIN_LENGTH)
                return null;

            // The last index of requiredRefBases within scBases.
            int scIndexMatchEnd = scMatchIndexStart + SC_INSERT_MIN_SC_LENGTH - 1;

            // How many bases in the gap/insert?
            int gapSize = scLength - scIndexMatchEnd - 1;

            try
            {
                // Ref and alt for the insert.
                String ref = readBases.substring(scIndexMatchEnd, scIndexMatchEnd + 1);
                String alt = readBases.substring(scIndexMatchEnd, scIndexMatchEnd + gapSize + 1);

                // Null ReadContext. Why?
                return new AltRead_(null, ref, alt, 0, 0, false, null);
            }
            catch(Exception e)
            {
                return null;
            }
        }
        else
        {
            // Similar but now the soft clip on the right side.
            int refIndexAtReadEnd = refBases.Index + readEnd - refBases.Position;
            int refIndexStart = refIndexAtReadEnd + 1;
            int refIndexEnd = refIndexStart + SC_INSERT_MIN_SC_LENGTH - 1;

            if(refIndexStart < 0 || refIndexEnd >= refBases.Bases.length)
                return null;

            String requiredRefBases = new String(refBases.Bases, refIndexStart, SC_INSERT_MIN_SC_LENGTH);

            String scBases = readBases.substring(scReadIndex);
            int scMatchIndex = scBases.indexOf(requiredRefBases);

            // We don't care in this case if the gap is too large like in the left case?
            if(scMatchIndex <= 0)
                return null;

            try
            {
                String ref = readBases.substring(scReadIndex - 1, scReadIndex);
                // Ref + gap/inserted bases.
                String alt = readBases.substring(scReadIndex - 1, scReadIndex + scMatchIndex);

                // Null ReadContext. Why?
                return new AltRead_(null, ref, alt, 0, 0, false, null);
            }
            catch(Exception e)
            {
                return null;
            }
        }
    }

    /**
     * Returns whether the readIndex is not near the ends (i.e. not in the potential flanks) of the read.
     */
    private boolean notWithinFlanks(int readIndex, final SAMRecord record)
    {
        // Read index cannot be in the flank, has to be within the
        return readIndex >= mConfig.ReadContextFlankSize && readIndex < record.getReadLength() - mConfig.ReadContextFlankSize;
    }

    /**
     * Returns min BaseQ from readIndex until readIndex + length (or earlier if we run over the end of the read).
     */
    private int minBaseQuality(int readIndex, final SAMRecord record, int length)
    {
        int maxIndex = Math.min(readIndex + length, record.getBaseQualities().length) - 1;
        int quality = Integer.MAX_VALUE;
        for(int i = readIndex; i <= maxIndex; i++)
        {
            quality = Math.min(quality, record.getBaseQualities()[i]);
        }
        return quality;
    }

    public int getReadCount()
    {
        return mReadCount;
    }

    private boolean isHotspotPosition(int position)
    {
        return mHotspotPositions.contains(position); }
}
