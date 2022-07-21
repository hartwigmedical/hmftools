package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.CompletionException;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.samtools.CigarHandler;
import com.hartwig.hmftools.common.samtools.CigarTraversal;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.RefSequence;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class BaseQualityRegionCounter implements CigarHandler
{
    private final SamReader mBamReader;
    private final ReferenceSequenceFile mRefGenome;
    private final SageConfig mConfig;
    private final BaseQualityResults mResults;

    private ChrBaseRegion mRegion;
    private IndexedBases mIndexedBases;

    private final Set<QualityCounter> mQualityCounts; // summarised counts with position removed
    private final Map<BaseQualityKey,Integer> mKeyCountsMap;
    private int mPurgeIndex;
    private int mMaxIndex;

    private BaseQualityData[] mBaseQualityData; // base-qual data by position for this region

    private final PerformanceCounter mPerfCounter;
    private int mReadCounter;

    private static final CigarElement SINGLE = new CigarElement(1, CigarOperator.M);
    private static final byte N = (byte) 'N';
    private static final int BASE_DATA_POS_BUFFER = 100;

    public BaseQualityRegionCounter(
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

            mKeyCountsMap.clear();
            mQualityCounts.clear();
            mReadCounter = 0;
            mPurgeIndex = 0;
            mMaxIndex = 0;
            mPerfCounter.reset();
        }
    }

    public Collection<QualityCounter> getQualityCounts() { return mQualityCounts; }

    public void run()
    {
        SG_LOGGER.trace("processing BQR region {}", mRegion);

        mPerfCounter.start();

        readBam();

        for(int i = mPurgeIndex; i <= mMaxIndex; ++i)
        {
            mapBaseQualityData(mBaseQualityData[i]);
        }

        for(Map.Entry<BaseQualityKey,Integer> entry : mKeyCountsMap.entrySet())
        {
            QualityCounter counter = new QualityCounter(entry.getKey());
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

        Map<BaseQualityKey,Integer> keyCounts = bqData.formKeyCounts(
                mConfig.QualityRecalibration.MaxAltCount, mConfig.QualityRecalibration.MaxAltPerc);

        for(Map.Entry<BaseQualityKey,Integer> entry : keyCounts.entrySet())
        {
            Integer count = mKeyCountsMap.get(entry.getKey());
            mKeyCountsMap.put(entry.getKey(), count != null ? count + entry.getValue() : entry.getValue());
        }
    }

    private void readBam()
    {
        if(mBamReader == null)
            return;

        BamSlicer slicer = new BamSlicer(mConfig.MinMapQuality);

        try
        {
            slicer.slice(mBamReader, Lists.newArrayList(mRegion), this::processRecord);
        }
        catch(Exception e)
        {
            throw new CompletionException(e);
        }
    }

    public void processRecord(@NotNull final SAMRecord record)
    {
        ++mReadCounter;
        CigarTraversal.traverseCigar(record, this);

        if(mReadCounter > 0 && (mReadCounter % 1000) == 0)
            purgeBaseDataList(record.getAlignmentStart());
    }

    @Override
    public void handleInsert(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int readIndex, final int refPos)
    {
        // need to add one because indel is actually AFTER this by convention
        int indelPos = refPos + 1;
        handleAlignment(record, SINGLE, readIndex, refPos);
        markIndelPosition(indelPos);
    }

    @Override
    public void handleDelete(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int readIndex, final int refPos)
    {
        int indelPos = refPos + 1;
        handleAlignment(record, SINGLE, readIndex, refPos);
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
    public void handleAlignment(final SAMRecord record, final CigarElement cigarElement, final int startReadIndex, final int refPos)
    {
        for(int i = 0; i < cigarElement.getLength(); i++)
        {
            int position = refPos + i;

            if(position > mRegion.end())
                return;

            if(position < mRegion.start())
                continue;

            int readIndex = startReadIndex + i;

            byte ref = mIndexedBases.base(position);
            byte alt = record.getReadBases()[readIndex];
            byte quality = record.getBaseQualities()[readIndex];
            byte[] trinucleotideContext = mIndexedBases.trinucleotideContext(position);

            if(alt == N || !isValid(trinucleotideContext))
                continue;

            BaseQualityData baseQualityData = getOrCreateBaseQualData(position, ref, trinucleotideContext);
            baseQualityData.processRead(alt, quality);
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
            if(b == N)
                return false;
        }

        return trinucleotideContext.length == 3;
    }
}
