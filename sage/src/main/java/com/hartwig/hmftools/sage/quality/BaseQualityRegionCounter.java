package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.CompletionException;
import java.util.stream.Collectors;

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
    private final ChrBaseRegion mRegion;
    private final ReferenceSequenceFile mRefGenome;
    private final IndexedBases mIndexedBases;
    private final SageConfig mConfig;
    private final BaseQualityResults mResults;

    private final Set<Integer> mIndelPositions = Sets.newHashSet();

    private final Set<QualityCounter> mQualityCounts; // summarised counts with position removed

    private BaseQualityData[] mBaseQualityData; // base-qual data by position for this region

    private static final CigarElement SINGLE = new CigarElement(1, CigarOperator.M);
    private static final byte N = (byte) 'N';

    private final PerformanceCounter mPerfCounter;

    public BaseQualityRegionCounter(
            final SageConfig config, final SamReader bamReader, final ReferenceSequenceFile refGenome, final ChrBaseRegion region,
            final BaseQualityResults results)
    {
        mConfig = config;
        mBamReader = bamReader;

        mRegion = region;
        mRefGenome = refGenome;
        mResults = results;

        if(mRefGenome != null)
        {
            final RefSequence refSequence = new RefSequence(mRegion, mRefGenome);
            mIndexedBases = refSequence.alignment();
        }
        else
        {
            mIndexedBases = null;
        }

        int regionPositionCount = region.baseLength();
        mBaseQualityData = new BaseQualityData[regionPositionCount];

        mQualityCounts = Sets.newHashSet();

        mPerfCounter = new PerformanceCounter("BaseQualBuild");
    }

    public Collection<QualityCounter> getQualityCounts() { return mQualityCounts; }

    public void run()
    {
        SG_LOGGER.trace("processing BQR region {}", mRegion);

        mPerfCounter.start();

        readBam();

        Map<BaseQualityKey,Integer> countsMap = Maps.newHashMap();

        for(BaseQualityData bqData : mBaseQualityData)
        {
            if(bqData == null)
                continue;

            if(bqData.hasIndel())
                continue;

            Map<BaseQualityKey,Integer> keyCounts = bqData.formKeyCounts(mConfig.QualityRecalibration.MaxAltCount);

            for(Map.Entry<BaseQualityKey,Integer> entry : keyCounts.entrySet())
            {
                Integer count = countsMap.get(entry.getKey());
                countsMap.put(entry.getKey(), count != null ? count + entry.getValue() : entry.getValue());
            }
        }

        for(Map.Entry<BaseQualityKey,Integer> entry : countsMap.entrySet())
        {
            QualityCounter counter = new QualityCounter(entry.getKey());
            counter.increment(entry.getValue());
            mQualityCounts.add(counter);
        }

        mPerfCounter.stop();

        mResults.addBaseQualityRegionCounter(this);
        mResults.addPerfCounter(mPerfCounter);

        mBaseQualityData = null;
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
        CigarTraversal.traverseCigar(record, this);
    }

    @Override
    public void handleInsert(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int readIndex, final int refPos)
    {
        // need to add one because indel is actually AFTER this by convention
        int indelPos = refPos + 1;
        mIndelPositions.add(indelPos);
        handleAlignment(record, SINGLE, readIndex, refPos);
        markIndelPosition(indelPos);
    }

    @Override
    public void handleDelete(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int readIndex, final int refPos)
    {
        int indelPos = refPos + 1;
        mIndelPositions.add(indelPos);
        handleAlignment(record, SINGLE, readIndex, refPos);
        markIndelPosition(indelPos);
    }

    private void markIndelPosition(int position)
    {
        if(!mRegion.containsPosition(position))
            return;

        int posIndex = position - mRegion.start();

        BaseQualityData baseQualityData = mBaseQualityData[posIndex];

        if(baseQualityData == null)
        {
            byte ref = mIndexedBases.base(position);
            byte[] trinucleotideContext = mIndexedBases.trinucleotideContext(position);

            baseQualityData = new BaseQualityData(position, ref, trinucleotideContext);
            mBaseQualityData[posIndex] = baseQualityData;
        }

        baseQualityData.setHasIndel();
    }

    @Override
    public void handleAlignment(final SAMRecord record, final CigarElement cigarElement, final int startReadIndex, final int refPos)
    {
        for(int i = 0; i < cigarElement.getLength(); i++)
        {
            int readIndex = startReadIndex + i;
            int position = refPos + i;

            if(position > mRegion.end())
                return;

            if(position < mRegion.start())
                continue;

            byte ref = mIndexedBases.base(position);
            byte alt = record.getReadBases()[readIndex];
            byte quality = record.getBaseQualities()[readIndex];
            byte[] trinucleotideContext = mIndexedBases.trinucleotideContext(position);

            if(alt == N || !isValid(trinucleotideContext))
                continue;

            // missing old key(var(A->G) cxt(CAA) qual(37)) new count(139)
            //if((char)ref != 'A' || (char)alt != 'G' || !new String(trinucleotideContext).equals("CAA") || (int)quality != 37)
            //    continue;

            BaseQualityData baseQualityData = getOrCreateBaseQualData(position, ref, trinucleotideContext);
            baseQualityData.processRead(alt, quality);

            /*
            Map<BaseQualityKey,Integer> posCounts = mQualityMap.get(position);

            if(posCounts == null)
            {
                posCounts = Maps.newHashMap();
                mQualityMap.put(position, posCounts);
            }

            boolean matched = false;
            for(Map.Entry<BaseQualityKey,Integer> entry : posCounts.entrySet())
            {
                BaseQualityKey key = entry.getKey();

                if(key.matches(ref, alt, quality, trinucleotideContext))
                {
                    entry.setValue(entry.getValue() + 1);
                    matched = true;
                    break;
                }
            }

            if(!matched)
            {
                posCounts.put(new BaseQualityKey(ref, alt, trinucleotideContext, quality), 1);
            }
            */
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
