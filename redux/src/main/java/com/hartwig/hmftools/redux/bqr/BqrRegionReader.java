package com.hartwig.hmftools.redux.bqr;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_N_BYTE;
import static com.hartwig.hmftools.common.redux.BqrReadType.extractReadType;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ULTIMA;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConstants.BQR_MIN_MAP_QUAL;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.bam.CigarHandler;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.perf.PerformanceCounter;
import com.hartwig.hmftools.common.redux.BqrKey;
import com.hartwig.hmftools.common.redux.BqrReadType;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.sequencing.SequencingType;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class BqrRegionReader implements CigarHandler
{
    private final RefGenomeInterface mRefGenome;
    private final BaseQualityResults mResults;

    private final List<ChrBaseRegion> mAllBqrRegions;

    private boolean mHasActiveRegion;
    private ChrBaseRegion mPartitionOverlapRegion; // overlap between a BQR region and the current partition
    private RefSequence mCurrentRefSequence;

    private final Set<BqrKeyCounter> mQualityCounts; // summarised counts with position removed
    private final Map<BqrKey, Integer> mKeyCountsMap;

    private int mNextPurgePosition;
    private int mPurgeIndex;
    private int mMaxIndex;

    private BaseQualityData[] mBaseQualityData; // base-qual data by position for this region

    private final PerformanceCounter mPerfCounter;

    private long mTotalReadsUsed;
    private int mTotalAltsFiltered;

    private BqrReadType mCurrentReadType;
    private final SequencingType mSequencingType;

    public static final int REF_BASE_REGION_SIZE = 100000;
    private static final int BASE_DATA_PURGE_POS_BUFFER = 10;
    private static final int BASE_DATA_PURGE_POS_CHECK = 1000;

    public BqrRegionReader(
            final SequencingType sequencingType, final RefGenomeInterface refGenome, final BaseQualityResults results,
            final List<ChrBaseRegion> allRegions)
    {
        mRefGenome = refGenome;
        mResults = results;
        mAllBqrRegions = allRegions;

        mPartitionOverlapRegion = null;
        mCurrentRefSequence = null;
        mHasActiveRegion = false;

        mCurrentReadType = BqrReadType.NONE;
        mSequencingType = sequencingType;

        mBaseQualityData = null;
        mQualityCounts = Sets.newHashSet();
        mKeyCountsMap = Maps.newHashMap();
        mPurgeIndex = 0;
        mNextPurgePosition = 0;
        mMaxIndex = 0;

        mPerfCounter = new PerformanceCounter("BaseQualBuild");

        mTotalReadsUsed = 0;
        mTotalAltsFiltered = 0;
    }

    public boolean isActive() { return mHasActiveRegion; }
    public long totalReadsUsed() { return mTotalReadsUsed; }
    public int totalAltsFiltered() { return mTotalAltsFiltered; }

    public void initialise(final ChrBaseRegion region)
    {
        if(mRefGenome == null)
            return;

        List<ChrBaseRegion> overlappedBqrRegions = mAllBqrRegions.stream().filter(x -> x.overlaps(region)).collect(Collectors.toList());

        if(overlappedBqrRegions.isEmpty())
        {
            mHasActiveRegion = false;
            mPartitionOverlapRegion = null;
            mCurrentRefSequence = null;
            return;
        }

        ChrBaseRegion firstRegion = overlappedBqrRegions.get(0);
        ChrBaseRegion lastRegion = overlappedBqrRegions.get(overlappedBqrRegions.size() - 1);
        mPartitionOverlapRegion = new ChrBaseRegion(firstRegion.Chromosome, firstRegion.start(), lastRegion.end());
        mHasActiveRegion = true;
    }

    private void checkLocalRegion(int position)
    {
        // work in blocks of 100K (could consider even less for panels)
        if(mCurrentRefSequence != null && mCurrentRefSequence.positionWithinBounds(position))
            return;

        if(position < mPartitionOverlapRegion.start())
            return;

        if(mHasActiveRegion)
            processCompletedLocalRegion();

        if(position > mPartitionOverlapRegion.end())
        {
            // de-activate if the reads are now outside the required BQR region(s)
            mHasActiveRegion = false;
            return;
        }

        // extract the next block of ref bases
        int refBaseStart = position;
        int refBaseEnd = min(position + REF_BASE_REGION_SIZE, mPartitionOverlapRegion.end());

        if(refBaseEnd <= refBaseStart)
            return;

        mCurrentRefSequence = new RefSequence(mPartitionOverlapRegion.Chromosome, refBaseStart, refBaseEnd, mRefGenome);

        if(!mCurrentRefSequence.IsValid)
            return;

        mNextPurgePosition = refBaseStart + BASE_DATA_PURGE_POS_CHECK;

        if(mBaseQualityData == null || mBaseQualityData.length != mCurrentRefSequence.length())
        {
            mBaseQualityData = new BaseQualityData[mCurrentRefSequence.length()];
        }

        RD_LOGGER.trace("bqr starting  new sub-region({})", mCurrentRefSequence);

        mPerfCounter.reset();
        mPerfCounter.start();
    }

    private void processCompletedLocalRegion()
    {
        if(mBaseQualityData == null)
            return;

        RD_LOGGER.trace("bqr completed sub-region({})", mCurrentRefSequence);

        buildQualityCounts();

        mResults.addBaseQualityRegionCounter(this);
        mResults.addPerfCounter(mPerfCounter);

        for(int i = 0; i < mBaseQualityData.length; ++i)
        {
            mBaseQualityData[i] = null;
        }

        mKeyCountsMap.clear();
        mQualityCounts.clear();
        mNextPurgePosition = 0;
        mPurgeIndex = 0;
        mMaxIndex = 0;

        mPerfCounter.stop();

        /*
        if(mConfig.PerfWarnTime > 0 && mPerfCounter.getLastTime() > mConfig.PerfWarnTime)
        {
            SG_LOGGER.warn("BQR region({}) time({}) reads({})",
                    mRegion, String.format("%.1f", mPerfCounter.getLastTime()), mReadCounter);
        }
        */

    }

    public Collection<BqrKeyCounter> getQualityCounts() { return mQualityCounts; }

    public void onRegionComplete()
    {
        processCompletedLocalRegion();
    }

    @VisibleForTesting
    public void buildQualityCounts()
    {
        for(int i = mPurgeIndex; i <= mMaxIndex; ++i)
        {
            BaseQualityData baseQualityData = mBaseQualityData[i];

            if(baseQualityData != null)
                convertToKeyCounts(baseQualityData);
        }

        for(Map.Entry<BqrKey, Integer> entry : mKeyCountsMap.entrySet())
        {
            BqrKeyCounter counter = new BqrKeyCounter(entry.getKey());
            counter.increment(entry.getValue());
            mQualityCounts.add(counter);
        }
    }

    private void convertToKeyCounts(final BaseQualityData bqData)
    {
        if(bqData == null)
            return;

        if(bqData.hasIndel())
            return;

        Map<BqrKey,Integer> keyCounts = bqData.formKeyCounts();

        mTotalAltsFiltered += bqData.filteredAltCount();

        for(Map.Entry<BqrKey,Integer> entry : keyCounts.entrySet())
        {
            if(!entry.getKey().isValid())
                continue;

            Integer count = mKeyCountsMap.get(entry.getKey());
            mKeyCountsMap.put(entry.getKey(), count != null ? count + entry.getValue() : entry.getValue());
        }
    }

    public void processRecord(final SAMRecord record)
    {
        checkLocalRegion(record.getAlignmentStart());

        if(!mHasActiveRegion || mCurrentRefSequence == null || !mCurrentRefSequence.IsValid)
            return;

        if(record.getMappingQuality() < BQR_MIN_MAP_QUAL)
            return;

        setShortFragmentBoundaries(record);

        mCurrentReadType = extractReadType(record, mSequencingType);

        CigarHandler.traverseCigar(record, this);

        if(record.getAlignmentStart() >= mNextPurgePosition)
        {
            purgeBaseDataList(record.getAlignmentStart());
        }
    }

    private static final int SHORT_FRAG_BOUNDARY_NONE = -1;
    private int mMinReadStartPosition = SHORT_FRAG_BOUNDARY_NONE;
    private int mMaxReadEndPosition = SHORT_FRAG_BOUNDARY_NONE;

    private void setShortFragmentBoundaries(final SAMRecord record)
    {
        // ensure we don’t read past insert size on 3’ end of short fragments
        mMinReadStartPosition = SHORT_FRAG_BOUNDARY_NONE;
        mMaxReadEndPosition = SHORT_FRAG_BOUNDARY_NONE;

        // CHECK: was using Sage config getReadLength()
        if(record.getInferredInsertSize() == 0 || abs(record.getInferredInsertSize()) >= record.getReadBases().length)
            return;

        if(record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag())
            return;

        if(record.getReadNegativeStrandFlag())
        {
            // set boundary at start
            mMinReadStartPosition = record.getMateAlignmentStart();
        }
        else
        {
            // boundary set at end
            String mateCigar = record.getStringAttribute(MATE_CIGAR_ATTRIBUTE);

            if(mateCigar != null)
            {
                mMaxReadEndPosition = getMateAlignmentEnd(record.getMateAlignmentStart(), mateCigar);
            }
            else
            {
                mMaxReadEndPosition = record.getAlignmentStart() + abs(abs(record.getInferredInsertSize())) - 1;
            }
        }
    }

    @Override
    public void handleInsert(final SAMRecord record, final CigarElement e, final int readIndex, final int refPos)
    {
        // note: ref position here is the last base of the previous aligned element - likewise for deletes
        markIndelPosition(refPos);
    }

    @Override
    public void handleDelete(final SAMRecord record, final CigarElement e, final int readIndex, final int refPos)
    {
        markIndelPosition(refPos);
    }

    private void markIndelPosition(int position)
    {
        if(!mCurrentRefSequence.positionWithinBounds(position))
            return;

        byte ref = mCurrentRefSequence.base(position);
        byte[] trinucleotideContext = mCurrentRefSequence.trinucleotideContext(position);
        BaseQualityData bqData = getOrCreateBaseQualData(position, ref, trinucleotideContext);
        bqData.setHasIndel();
    }

    // debug only
    private static final boolean LOG_READ_INFO = false;
    private static final List<String> LOG_TNCS = List.of("AAA");
    private static final List<Byte> LOG_QUAL = List.of((byte)11);
    private static final List<Byte> LOG_POSITIONS = List.of();
    private static final List<BqrReadType> LOG_READ_TYPES = List.of();

    @Override
    public void handleAlignment(final SAMRecord record, final CigarElement cigarElement, final int startReadIndex, final int refPos)
    {
        byte[] trinucleotideContext = new byte[3];
        boolean readPosStrand = !record.getReadNegativeStrandFlag();
        boolean readUsed = false;

        for(int i = 0; i < cigarElement.getLength(); i++)
        {
            int position = refPos + i;

            if(mCurrentRefSequence.afterEnd(position))
                return;

            if(mCurrentRefSequence.beforeStart(position))
                continue;

            if(mMinReadStartPosition > 0 && position < mMinReadStartPosition)
                continue;

            if(mMaxReadEndPosition > 0 && position > mMaxReadEndPosition)
                break;

            int readIndex = startReadIndex + i;

            byte ref = mCurrentRefSequence.base(position);
            byte alt = record.getReadBases()[readIndex];

            if(alt == DNA_N_BYTE)
                continue;

            byte quality = record.getBaseQualities()[readIndex];

            if(mSequencingType == ULTIMA && quality != ULTIMA_MAX_QUAL)
                continue;

            mCurrentRefSequence.populateTrinucleotideContext(position, trinucleotideContext);
            readUsed = true;

            if(LOG_READ_INFO)
            {
                String tncStr = new String(trinucleotideContext);

                if((LOG_QUAL.isEmpty() || LOG_QUAL.contains(quality))
                || (LOG_POSITIONS.isEmpty() || LOG_POSITIONS.contains(position))
                || (LOG_READ_TYPES.isEmpty() || LOG_READ_TYPES.contains(mCurrentReadType))
                || (LOG_TNCS.isEmpty() || LOG_TNCS.contains(tncStr)))
                {
                    RD_LOGGER.debug("BQR: read({}) position({}:{}) context({}) alt({}) qual({}) readType({}) posStrand({})",
                            record.getReadName(), mPartitionOverlapRegion.Chromosome, position, tncStr, alt, quality, mCurrentReadType, readPosStrand);
                }
            }

            BaseQualityData baseQualityData = getOrCreateBaseQualData(position, ref, trinucleotideContext);
            baseQualityData.processReadBase(mCurrentReadType, alt, quality, readPosStrand);
        }

        if(readUsed)
            ++mTotalReadsUsed;
    }

    private void purgeBaseDataList(int currentReadStartPos)
    {
        for(; mPurgeIndex <= mMaxIndex; ++mPurgeIndex)
        {
            int position = mCurrentRefSequence.position(mPurgeIndex);

            if(position <= 0)
            {
                RD_LOGGER.error("invalid calc position({}) purgeIndex({} -> {}) readPos({}) refSeq({})",
                        position, mPurgeIndex, mMaxIndex, currentReadStartPos, mCurrentRefSequence);
                System.exit(1);
            }

            if(position >= currentReadStartPos - BASE_DATA_PURGE_POS_BUFFER)
                break;

            BaseQualityData bqData = mBaseQualityData[mPurgeIndex];

            if(bqData == null)
                continue;

            convertToKeyCounts(bqData);
            mBaseQualityData[mPurgeIndex] = null;
        }

        mNextPurgePosition = currentReadStartPos + BASE_DATA_PURGE_POS_CHECK;
    }

    @VisibleForTesting
    public BaseQualityData getOrCreateBaseQualData(int position, final byte ref, final byte[] trinucleotideContext)
    {
        int posIndex = mCurrentRefSequence.index(position);
        BaseQualityData baseQualityData = mBaseQualityData[posIndex];

        if(baseQualityData == null)
        {
            baseQualityData = new BaseQualityData(trinucleotideContext);
            mBaseQualityData[posIndex] = baseQualityData;
            mMaxIndex = posIndex;
        }

        return baseQualityData;
    }

    @VisibleForTesting
    public void initialise(final ChrBaseRegion region, final RefGenomeInterface refGenome)
    {
        mPartitionOverlapRegion = region;
        mCurrentRefSequence = new RefSequence(region.Chromosome, region.start(), region.end(), refGenome);

        mBaseQualityData = new BaseQualityData[mCurrentRefSequence.length()];

        mKeyCountsMap.clear();
        mQualityCounts.clear();
        mPurgeIndex = 0;
        mMaxIndex = 0;
    }
}
