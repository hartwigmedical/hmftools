package com.hartwig.hmftools.redux.bqr;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractUmiType;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_N_BYTE;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ULTIMA;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractConsensusType;
import static com.hartwig.hmftools.redux.common.Constants.BQR_MIN_MAP_QUAL;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.bam.CigarHandler;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.perf.PerformanceCounter;
import com.hartwig.hmftools.common.qual.BqrKey;
import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.sequencing.SequencingType;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class BqrRegionReader implements CigarHandler
{
    private final RefGenomeInterface mRefGenome;
    private final BaseQualityResults mResults;

    private final List<ChrBaseRegion> mAllBqrRegions;
    private ChrBaseRegion mCurrentRegion;
    private RefSequence mCurrentRefSequence;

    private final Set<BqrKeyCounter> mQualityCounts; // summarised counts with position removed
    private final Map<BqrKey, Integer> mKeyCountsMap;
    private int mPurgeIndex;
    private int mMaxIndex;

    private BaseQualityDataCollection[] mBaseQualityData; // base-qual data by position for this region

    private final PerformanceCounter mPerfCounter;
    private int mReadCounter;

    private BqrReadType mCurrentReadType;
    private final SequencingType mSequencingType;

    private static final int BASE_DATA_POS_BUFFER = 100;
    public static final int REF_BASE_REGION_SIZE = 100000;

    public BqrRegionReader(
            final SequencingType sequencingType, final RefGenomeInterface refGenome, final BaseQualityResults results,
            final List<ChrBaseRegion> allRegions)
    {
        mRefGenome = refGenome;
        mResults = results;
        mAllBqrRegions = allRegions;
        mCurrentRegion = null;
        mCurrentRefSequence = null;

        mCurrentReadType = BqrReadType.NONE;
        mSequencingType = sequencingType;

        mBaseQualityData = null;
        mQualityCounts = Sets.newHashSet();
        mKeyCountsMap = Maps.newHashMap();
        mPurgeIndex = 0;
        mMaxIndex = 0;

        mPerfCounter = new PerformanceCounter("BaseQualBuild");
        mReadCounter = 0;
    }

    public boolean isActive() { return mCurrentRegion != null; }

    public void initialise(final ChrBaseRegion region)
    {
        if(mRefGenome == null)
            return;

        mCurrentRegion = mAllBqrRegions.stream().filter(x -> x.overlaps(region)).findFirst().orElse(null);

        if(mCurrentRegion == null)
        {
            mCurrentRefSequence = null;
            return;
        }

        mCurrentRefSequence = new RefSequence(
                mCurrentRegion.Chromosome, mCurrentRegion.start(), mCurrentRegion.start() + REF_BASE_REGION_SIZE, mRefGenome);

        if(mBaseQualityData == null || mBaseQualityData.length != region.baseLength())
        {
            int regionPositionCount = region.baseLength();
            mBaseQualityData = new BaseQualityDataCollection[regionPositionCount];
        }
        else
        {
            for(int i = 0; i < mBaseQualityData.length; ++i)
            {
                mBaseQualityData[i] = null;
            }
        }

        mKeyCountsMap.clear();
        mQualityCounts.clear();
        mReadCounter = 0;
        mPurgeIndex = 0;
        mMaxIndex = 0;

        mPerfCounter.reset();
    }

    public Collection<BqrKeyCounter> getQualityCounts() { return mQualityCounts; }

    public void onRegionComplete()
    {
        if(mBaseQualityData == null)
            return;

        buildQualityCounts();

        mKeyCountsMap.clear();

        mPerfCounter.stop();

        /*
        if(mConfig.PerfWarnTime > 0 && mPerfCounter.getLastTime() > mConfig.PerfWarnTime)
        {
            SG_LOGGER.warn("BQR region({}) time({}) reads({})",
                    mRegion, String.format("%.1f", mPerfCounter.getLastTime()), mReadCounter);
        }
        */

        mResults.addBaseQualityRegionCounter(this);
        mResults.addPerfCounter(mPerfCounter);
    }

    @VisibleForTesting
    protected void buildQualityCounts()
    {
        for(int i = mPurgeIndex; i <= mMaxIndex; ++i)
        {
            BaseQualityDataCollection bqDataCollection = mBaseQualityData[i];

            if(bqDataCollection == null)
                continue;

            buildSummaryData(bqDataCollection);
        }

        for(Map.Entry<BqrKey, Integer> entry : mKeyCountsMap.entrySet())
        {
            BqrKeyCounter counter = new BqrKeyCounter(entry.getKey());
            counter.increment(entry.getValue());
            mQualityCounts.add(counter);
        }
    }

    private void buildSummaryData(final BaseQualityDataCollection bqDataCollection)
    {
        if(bqDataCollection != null)
        {
            bqDataCollection.DataMap.values().forEach(x -> buildSummaryData(x));
        }
    }

    private void buildSummaryData(final BaseQualityData bqData)
    {
        if(bqData == null)
            return;

        if(bqData.hasIndel())
            return;

        Map<BqrKey, Integer> keyCounts = bqData.formKeyCounts();

        for(Map.Entry<BqrKey, Integer> entry : keyCounts.entrySet())
        {
            Integer count = mKeyCountsMap.get(entry.getKey());
            mKeyCountsMap.put(entry.getKey(), count != null ? count + entry.getValue() : entry.getValue());
        }
    }

    public void processRecord(final SAMRecord record)
    {
        if(record.getAlignmentStart() < mCurrentRegion.start())
            return;

        if(record.getAlignmentStart() > mCurrentRegion.end())
        {
            mCurrentRegion = null;
            mCurrentRefSequence = null;
            return;
        }

        // check need to update the current ref sequence
        if(record.getAlignmentStart() >= mCurrentRefSequence.End - 3) // a buffer for TNC
        {
            mCurrentRefSequence = new RefSequence(
                    mCurrentRegion.Chromosome, mCurrentRefSequence.End, mCurrentRefSequence.End + REF_BASE_REGION_SIZE, mRefGenome);
        }

        if(record.getMappingQuality() < BQR_MIN_MAP_QUAL)
            return;

        ++mReadCounter;

        setShortFragmentBoundaries(record);

        mCurrentReadType = extractReadType(record, mSequencingType);

        CigarHandler.traverseCigar(record, this);

        if(mReadCounter > 0 && (mReadCounter % 1000) == 0)
        {
            purgeBaseDataList(record.getAlignmentStart());
        }
    }

    public static BqrReadType extractReadType(final SAMRecord record, final SequencingType sequencingType)
    {
        if(sequencingType == SequencingType.ILLUMINA)
            return BqrReadType.fromUmiType(extractUmiType(record));
        else
            return BqrReadType.fromUltimaType(extractConsensusType(record));
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
        if(!mCurrentRegion.containsPosition(position))
            return;

        byte ref = mCurrentRefSequence.base(position);
        byte[] trinucleotideContext = mCurrentRefSequence.trinucleotideContext(position);
        BaseQualityData bqData = getOrCreateBaseQualData(position, ref, trinucleotideContext, mCurrentReadType);
        bqData.setHasIndel();
    }

    @Override
    public void handleAlignment(final SAMRecord record, final CigarElement cigarElement, final int startReadIndex, final int refPos)
    {
        for(int i = 0; i < cigarElement.getLength(); i++)
        {
            int position = refPos + i;

            if(position > mCurrentRegion.end())
                return;

            if(position < mCurrentRegion.start())
                continue;

            if(mMinReadStartPosition > 0 && position < mMinReadStartPosition)
                continue;

            if(mMaxReadEndPosition > 0 && position > mMaxReadEndPosition)
                break;

            int readIndex = startReadIndex + i;

            byte ref = mCurrentRefSequence.base(position);
            byte alt = record.getReadBases()[readIndex];
            byte quality = record.getBaseQualities()[readIndex];

            if(mSequencingType == ULTIMA && quality != ULTIMA_MAX_QUAL)
                continue;

            byte[] trinucleotideContext = mCurrentRefSequence.trinucleotideContext(position);

            if(alt == DNA_N_BYTE || !isValid(trinucleotideContext))
                continue;

            BaseQualityData baseQualityData = getOrCreateBaseQualData(position, ref, trinucleotideContext, mCurrentReadType);
            baseQualityData.processReadBase(alt, quality);
        }
    }

    private void purgeBaseDataList(int currentReadStartPos)
    {
        for(; mPurgeIndex <= mMaxIndex; ++mPurgeIndex)
        {
            BaseQualityDataCollection bqData = mBaseQualityData[mPurgeIndex];

            if(bqData == null)
                continue;

            if(bqData.Position >= currentReadStartPos - BASE_DATA_POS_BUFFER)
                break;

            buildSummaryData(bqData);
            mBaseQualityData[mPurgeIndex] = null;
        }
    }

    protected BaseQualityData getOrCreateBaseQualData(
            int position, final byte ref, final byte[] trinucleotideContext, final BqrReadType readType)
    {
        int posIndex = position - mCurrentRegion.start();
        BaseQualityDataCollection bqDataCollection = mBaseQualityData[posIndex];

        if(bqDataCollection == null)
        {
            bqDataCollection = new BaseQualityDataCollection(position);
            mBaseQualityData[posIndex] = bqDataCollection;
            mMaxIndex = posIndex;
        }

        BaseQualityData bqData = bqDataCollection.DataMap.get(readType);

        if(bqData == null)
        {
            bqData = new BaseQualityData(ref, trinucleotideContext, readType);
            bqDataCollection.DataMap.put(readType, bqData);
        }

        return bqData;
    }

    private static boolean isValid(final byte[] trinucleotideContext)
    {
        for(byte b : trinucleotideContext)
        {
            if(!Nucleotides.isValidDnaBase((char) b))
                return false;
        }

        return trinucleotideContext.length == 3;
    }

    private class BaseQualityDataCollection
    {
        public final int Position;
        public final Map<BqrReadType, BaseQualityData> DataMap;

        public BaseQualityDataCollection(final int position)
        {
            Position = position;
            DataMap = Maps.newHashMap();
        }
    }

    private class RefSequence
    {
        public final int Start;
        public final int End;
        public final byte[] Bases;

        private static final int BUFFER = 5;

        public RefSequence(final String chromosome, final int posStart, final int posEnd, final RefGenomeInterface refGenome)
        {
            final int sequenceEnd = refGenome.getChromosomeLength(chromosome);

            Start = max(1, posStart - BUFFER);
            End = min(sequenceEnd, posEnd) + BUFFER;
            Bases = refGenome.getBases(chromosome, Start, End);
        }

        public int index(int position)
        {
            return position - Start;
        }

        public int length()
        {
            return Bases.length;
        }

        public byte base(int position)
        {
            return Bases[index(position)];
        }

        public byte[] trinucleotideContext(int position)
        {
            int index = index(position);
            return new byte[] { Bases[index - 1], Bases[index], Bases[index + 1] };
        }

        public String toString()
        {
            return format("%d-%d len(%d)", Start, End, length());
        }

    }
}
