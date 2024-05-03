package com.hartwig.hmftools.sage.bqr;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.extractUmiType;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ULTIMA;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_MAX_QUAL;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractConsensusType;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.bqr.BqrConfig.useReadType;

import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.CompletionException;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.bam.CigarHandler;
import com.hartwig.hmftools.common.qual.BqrKey;
import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.RefSequence;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class BqrRegionReader implements CigarHandler
{
    private final SageConfig mConfig;
    private final ReferenceSequenceFile mRefGenome;
    private final SamReader mBamReader;
    private final BaseQualityResults mResults;
    private final BqrRecordWriter mRecordWriter;
    private final boolean mWritePositionData;
    private final boolean mWriteReadData;

    private ChrBaseRegion mRegion;
    private RefSequence mRefSequence;
    private final Set<Integer> mKnownVariants;

    private final Set<BqrKeyCounter> mQualityCounts; // summarised counts with position removed
    private final Map<BqrKey,Integer> mKeyCountsMap;
    private int mPurgeIndex;
    private int mMaxIndex;

    private BaseQualityDataCollection[] mBaseQualityData; // base-qual data by position for this region

    private final PerformanceCounter mPerfCounter;
    private int mReadCounter;

    private BqrReadType mCurrentReadType;
    private final boolean mUseReadType;
    private final SequencingType mSequencingType;

    private static final CigarElement SINGLE = new CigarElement(1, CigarOperator.M);
    private static final byte N = (byte) 'N';
    private static final byte M = (byte) 'M';
    private static final int BASE_DATA_POS_BUFFER = 100;

    public BqrRegionReader(
            final SageConfig config, final SamReader bamReader, final ReferenceSequenceFile refGenome, final BaseQualityResults results,
            final BqrRecordWriter recordWriter)
    {
        mConfig = config;
        mBamReader = bamReader;

        mRefGenome = refGenome;
        mResults = results;
        mRecordWriter = recordWriter;
        mWriteReadData = mConfig.BQR.WriteReads;
        mWritePositionData = mConfig.BQR.WritePositions;

        mUseReadType = useReadType(mConfig);
        mCurrentReadType = BqrReadType.NONE;
        mSequencingType = mConfig.Sequencing.Type;

        mBaseQualityData = null;
        mQualityCounts = Sets.newHashSet();
        mKnownVariants = Sets.newHashSet();
        mKeyCountsMap = Maps.newHashMap();
        mPurgeIndex = 0;
        mMaxIndex = 0;

        mPerfCounter = new PerformanceCounter("BaseQualBuild");
        mReadCounter = 0;
    }

    public void initialise(final ChrBaseRegion region, final Set<Integer> knownVariants)
    {
        mRegion = region;
        mKnownVariants.clear();
        mKnownVariants.addAll(knownVariants);

        if(mRefGenome != null)
        {
            mRefSequence = new RefSequence(mRegion, mRefGenome);
        }
        else
        {
            mRefSequence = null;
        }

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

    public void run()
    {
        // SG_LOGGER.trace("processing BQR region {}", mRegion);

        mPerfCounter.start();

        readBam();

        buildQualityCounts();

        mKeyCountsMap.clear();

        mPerfCounter.stop();

        if(mConfig.PerfWarnTime > 0 && mPerfCounter.getLastTime() > mConfig.PerfWarnTime)
        {
            SG_LOGGER.warn("BQR region({}) time({}) reads({})",
                    mRegion, String.format("%.1f", mPerfCounter.getLastTime()), mReadCounter);
        }

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

        for(Map.Entry<BqrKey,Integer> entry : mKeyCountsMap.entrySet())
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
            if(mWritePositionData)
            {
                for(Map.Entry<BqrReadType,BaseQualityData> entry : bqDataCollection.DataMap.entrySet())
                {
                    BqrReadType readType = entry.getKey();

                    if(readType != BqrReadType.DUAL)
                        continue;

                    BaseQualityData bqData = entry.getValue();

                    if(bqData.hasIndel())
                        continue;

                    for(AltQualityCount altQualityCount : bqData.altQualityCounts())
                    {
                        mRecordWriter.writePositionData(
                                mRegion.Chromosome, bqDataCollection.Position, bqData.Ref, altQualityCount.Alt,
                                bqData.TrinucleotideContext, altQualityCount.Quality, readType, altQualityCount.Count);
                    }
                }
            }

            bqDataCollection.DataMap.values().forEach(x -> buildSummaryData(x));
        }
    }

    private void buildSummaryData(final BaseQualityData bqData)
    {
        if(bqData == null)
            return;

        if(bqData.hasIndel())
            return;

        Map<BqrKey,Integer> keyCounts = bqData.formKeyCounts();

        for(Map.Entry<BqrKey,Integer> entry : keyCounts.entrySet())
        {
            Integer count = mKeyCountsMap.get(entry.getKey());
            mKeyCountsMap.put(entry.getKey(), count != null ? count + entry.getValue() : entry.getValue());
        }
    }

    private void readBam()
    {
        if(mBamReader == null)
            return;

        BamSlicer slicer = new BamSlicer(mConfig.BQR.MinMapQuality);

        try
        {
            slicer.slice(mBamReader, mRegion, this::processRecord);
        }
        catch(Exception e)
        {
            throw new CompletionException(e);
        }
    }

    public void processRecord(final SAMRecord record)
    {
        ++mReadCounter;

        setShortFragmentBoundaries(record);

        if(mUseReadType)
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

        if(record.getInferredInsertSize() == 0 || abs(record.getInferredInsertSize()) >= mConfig.getReadLength())
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
        if(!mRegion.containsPosition(position))
            return;

        byte ref = mRefSequence.base(position);
        byte[] trinucleotideContext = mRefSequence.trinucleotideContext(position);
        BaseQualityData bqData = getOrCreateBaseQualData(position, ref, trinucleotideContext, mCurrentReadType);
        bqData.setHasIndel();
    }

    @Override
    public void handleAlignment(final SAMRecord record, final CigarElement cigarElement, final int startReadIndex, final int refPos)
    {
        for(int i = 0; i < cigarElement.getLength(); i++)
        {
            int position = refPos + i;

            if(position > mRegion.end())
                return;

            if(position < mRegion.start())
                continue;

            if(mMinReadStartPosition > 0 && position < mMinReadStartPosition)
                continue;

            if(mMaxReadEndPosition > 0 && position > mMaxReadEndPosition)
                break;

            if(mKnownVariants.contains(position))
                continue;

            int readIndex = startReadIndex + i;

            byte ref = mRefSequence.base(position);
            byte alt = record.getReadBases()[readIndex];
            byte quality = record.getBaseQualities()[readIndex];

            if(mSequencingType == ULTIMA && quality != ULTIMA_MAX_QUAL)
                continue;

            byte[] trinucleotideContext = mRefSequence.trinucleotideContext(position);

            if(alt == N || !isValid(trinucleotideContext))
                continue;

            BaseQualityData baseQualityData = getOrCreateBaseQualData(position, ref, trinucleotideContext, mCurrentReadType);
            baseQualityData.processReadBase(alt, quality);

            if(mWriteReadData && ref != alt)
                mRecordWriter.writeRecordData(record, position, readIndex, ref, alt, trinucleotideContext, quality, mCurrentReadType);
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
        int posIndex = position - mRegion.start();
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
            if(!Nucleotides.isValidDnaBase((char)b))
                return false;
        }

        return trinucleotideContext.length == 3;
    }

    private class BaseQualityDataCollection
    {
        public final int Position;
        public final Map<BqrReadType,BaseQualityData> DataMap;

        public BaseQualityDataCollection(final int position)
        {
            Position = position;
            DataMap = Maps.newHashMap();
        }
    }
}
