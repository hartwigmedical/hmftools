package com.hartwig.hmftools.sage.bqr;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.samtools.CigarUtils.getEndPosition;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_MIN_MAP_QUALITY;

import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.CompletionException;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.samtools.CigarHandler;
import com.hartwig.hmftools.common.samtools.CigarTraversal;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.RefSequence;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class BqrRegionReader implements CigarHandler
{
    private final SamReader mBamReader;
    private final ReferenceSequenceFile mRefGenome;
    private final SageConfig mConfig;
    private final BaseQualityResults mResults;

    private ChrBaseRegion mRegion;
    private IndexedBases mIndexedBases;

    private final Set<BqrKeyCounter> mQualityCounts; // summarised counts with position removed
    private final Map<BqrKey,Integer> mKeyCountsMap;
    private int mPurgeIndex;
    private int mMaxIndex;

    private BaseQualityData[] mBaseQualityData; // base-qual data by position for this region

    private final PerformanceCounter mPerfCounter;
    private int mReadCounter;

    private static final CigarElement SINGLE = new CigarElement(1, CigarOperator.M);
    private static final byte N = (byte) 'N';
    private static final byte M = (byte) 'M';
    private static final int BASE_DATA_POS_BUFFER = 100;

    public BqrRegionReader(
            final SageConfig config, final SamReader bamReader, final ReferenceSequenceFile refGenome, final BaseQualityResults results)
    {
        mConfig = config;
        mBamReader = bamReader;

        mRefGenome = refGenome;
        mResults = results;

        mBaseQualityData = null;
        mQualityCounts = Sets.newHashSet();
        mKeyCountsMap = Maps.newHashMap();
        mPurgeIndex = 0;
        mMaxIndex = 0;

        mPerfCounter = new PerformanceCounter("BaseQualBuild");
        mReadCounter = 0;
    }

    public void initialise(final ChrBaseRegion region)
    {
        mRegion = region;

        if(mRefGenome != null)
        {
            final RefSequence refSequence = new RefSequence(mRegion, mRefGenome);
            mIndexedBases = refSequence.alignment();
        }
        else
        {
            mIndexedBases = null;
        }

        if(mBaseQualityData == null || mBaseQualityData.length != region.baseLength())
        {
            int regionPositionCount = region.baseLength();
            mBaseQualityData = new BaseQualityData[regionPositionCount];
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

        for(int i = mPurgeIndex; i <= mMaxIndex; ++i)
        {
            mapBaseQualityData(mBaseQualityData[i]);
        }

        for(Map.Entry<BqrKey,Integer> entry : mKeyCountsMap.entrySet())
        {
            BqrKeyCounter counter = new BqrKeyCounter(entry.getKey());
            counter.increment(entry.getValue());
            mQualityCounts.add(counter);
        }

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

    private void mapBaseQualityData(final BaseQualityData bqData)
    {
        if(bqData == null)
            return;

        if(bqData.hasIndel())
            return;

        Map<BqrKey,Integer> keyCounts = bqData.formKeyCounts(
                mConfig.QualityRecalibration.MaxAltCount, mConfig.QualityRecalibration.MaxAltPerc);

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

        BamSlicer slicer = new BamSlicer(DEFAULT_MIN_MAP_QUALITY);

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


        CigarTraversal.traverseCigar(record, this);

        if(mReadCounter > 0 && (mReadCounter % 1000) == 0)
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
                mMaxReadEndPosition = getEndPosition(record.getMateAlignmentStart(), mateCigar, false, false);
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
        // need to add one because indel is actually AFTER this by convention
        int indelPos = refPos + 1;
        handleAlignment(record, SINGLE, false, readIndex, refPos);
        markIndelPosition(indelPos);
    }

    @Override
    public void handleDelete(final SAMRecord record, final CigarElement e, final int readIndex, final int refPos)
    {
        int indelPos = refPos + 1;
        handleAlignment(record, SINGLE, false, readIndex, refPos);
        markIndelPosition(indelPos);
    }

    private void markIndelPosition(int position)
    {
        if(!mRegion.containsPosition(position))
            return;

        byte ref = mIndexedBases.base(position);
        byte[] trinucleotideContext = mIndexedBases.trinucleotideContext(position);
        BaseQualityData bqData = getOrCreateBaseQualData(position, ref, trinucleotideContext);
        bqData.setHasIndel();
    }

    @Override
    public void handleAlignment(final SAMRecord record, final CigarElement cigarElement, boolean beforeIndel, final int startReadIndex, final int refPos)
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

            int readIndex = startReadIndex + i;

            byte ref = mIndexedBases.base(position);
            byte alt = record.getReadBases()[readIndex];
            byte quality = record.getBaseQualities()[readIndex];
            byte[] trinucleotideContext = mIndexedBases.trinucleotideContext(position);

            if(alt == N || !isValid(trinucleotideContext))
                continue;

            BaseQualityData baseQualityData = getOrCreateBaseQualData(position, ref, trinucleotideContext);
            baseQualityData.processReadBase(alt, quality);
        }
    }

    private void purgeBaseDataList(int currentReadStartPos)
    {
        for(; mPurgeIndex <= mMaxIndex; ++mPurgeIndex)
        {
            BaseQualityData bqData = mBaseQualityData[mPurgeIndex];

            if(bqData == null)
                continue;

            if(bqData.Position >= currentReadStartPos - BASE_DATA_POS_BUFFER)
                break;

            mapBaseQualityData(bqData);
            mBaseQualityData[mPurgeIndex] = null;
        }
    }

    protected BaseQualityData getOrCreateBaseQualData(int position, final byte ref, final byte[] trinucleotideContext)
    {
        int posIndex = position - mRegion.start();
        BaseQualityData baseQualityData = mBaseQualityData[posIndex];

        if(baseQualityData == null)
        {
            baseQualityData = new BaseQualityData(position, ref, trinucleotideContext);
            mBaseQualityData[posIndex] = baseQualityData;
            mMaxIndex = posIndex;
        }

        return baseQualityData;
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
}
