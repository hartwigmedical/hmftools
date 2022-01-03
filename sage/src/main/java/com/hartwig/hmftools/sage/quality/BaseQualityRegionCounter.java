package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
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
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.read.IndexedBases;
import com.hartwig.hmftools.sage.ref.RefSequence;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.ReferenceSequenceFile;

class BaseQualityRegionCounter implements CigarHandler, Callable<Long>
{
    private final String mBamFile;
    private final ChrBaseRegion mRegion;
    private final ReferenceSequenceFile mRefGenome;
    private final IndexedBases mIndexedBases;
    private final SageConfig mConfig;

    private final Set<Integer> mIndelPositions = Sets.newHashSet();

    private final Set<QualityCounter> mQualityCounts; // summarised counts with position removed

    private final Map<Integer,Map<BaseQualityKey,Integer>> mQualityMap; // counts by position then location context

    private static final CigarElement SINGLE = new CigarElement(1, CigarOperator.M);
    private static final byte N = (byte) 'N';

    private int mReadCounter;
    private final PerformanceCounter mReadsPc;
    private final PerformanceCounter mSummaryPc;

    public BaseQualityRegionCounter(
            final SageConfig config, final String bamFile, final ReferenceSequenceFile refGenome, final ChrBaseRegion region)
    {
        mConfig = config;
        mBamFile = bamFile;

        mRegion = region;
        mRefGenome = refGenome;

        if(mRefGenome != null)
        {
            final RefSequence refSequence = new RefSequence(mRegion, mRefGenome);
            mIndexedBases = refSequence.alignment();
        }
        else
        {
            mIndexedBases = null;
        }

        mQualityMap = Maps.newHashMap();
        mQualityCounts = Sets.newHashSet();

        mReadCounter = 0;
        mReadsPc = new PerformanceCounter("BqrReads");
        mSummaryPc = new PerformanceCounter("BqrSummary");
    }

    @Override
    public Long call()
    {
        produceRegionCounts();
        return (long)0;
    }

    public Collection<QualityCounter> getQualityCounts() { return mQualityCounts; }

    protected Map<Integer,Map<BaseQualityKey,Integer>> getQualityMap() { return mQualityMap; }

    public void logPerfs()
    {
        SG_LOGGER.debug(String.format("region(%s) readCount(%d) time(read=%.3f summary=%.3f)",
                mRegion, mReadCounter, mReadsPc.getMaxTime(), mSummaryPc.getMaxTime()));
    }

    public void produceRegionCounts()
    {
        SG_LOGGER.trace("processing BQR region {}", mRegion);

        readBam();

        mSummaryPc.start();

        // remove locations where the alt count exceeds the configured limit
        Map<Integer,Set<BaseQualityKey>> repeatedAltLocations = findRepeatedAltLocations();

        if(!repeatedAltLocations.isEmpty())
        {
            SG_LOGGER.trace("region({}) has {} repeated alt locations", mRegion, repeatedAltLocations.size());
        }

        // form a set of counts by variant, no longer taking position into account
        Map<BaseQualityKey,Integer> countsMap = Maps.newHashMap();

        for(Map.Entry<Integer,Map<BaseQualityKey,Integer>> posEntry : mQualityMap.entrySet())
        {
            int position = posEntry.getKey();

            if(mIndelPositions.contains(position))
                continue;

            Set<BaseQualityKey> repeatedAlts = repeatedAltLocations.get(position);

            for(Map.Entry<BaseQualityKey,Integer> entry : posEntry.getValue().entrySet())
            {
                if(repeatedAlts != null && repeatedAlts.contains(altKey(entry.getKey())))
                {
                    SG_LOGGER.trace("skipped repeated alt location({})", entry.getKey());
                    continue;
                }

                BaseQualityKey key = entry.getKey();
                Integer count = countsMap.get(key);
                countsMap.put(key, count != null ? count + entry.getValue() : entry.getValue());
            }
        }

        for(Map.Entry<BaseQualityKey,Integer> entry : countsMap.entrySet())
        {
            QualityCounter counter = new QualityCounter(entry.getKey());
            counter.increment(entry.getValue());
            mQualityCounts.add(counter);
        }

        mSummaryPc.stop();
    }

    private void readBam()
    {
        if(mBamFile == null)
            return;

        mReadsPc.start();

        BamSlicer slicer = new BamSlicer(mConfig.MinMapQuality);

        try
        {
            final SamReader tumorReader = SamReaderFactory.makeDefault()
                    .validationStringency(mConfig.Stringency)
                    .referenceSource(new ReferenceSource(mRefGenome))
                    .open(new File(mBamFile));

            slicer.slice(tumorReader, Lists.newArrayList(mRegion), this::processRecord);
        }
        catch(Exception e)
        {
            throw new CompletionException(e);
        }

        mReadsPc.stop();
    }

    private Map<Integer,Set<BaseQualityKey>> findRepeatedAltLocations()
    {
        // at each position, find any alt repeated more than X times at any base quality or context
        Map<Integer,Set<BaseQualityKey>> repeatedAlts = Maps.newHashMap();

        for(Map.Entry<Integer,Map<BaseQualityKey,Integer>> posEntry : mQualityMap.entrySet())
        {
            int position = posEntry.getKey();

            Map<BaseQualityKey,Integer> altCounts = Maps.newHashMap();

            for(Map.Entry<BaseQualityKey,Integer> entry : posEntry.getValue().entrySet())
            {
                BaseQualityKey altKey = altKey(entry.getKey());

                Integer count = altCounts.get(altKey);
                int newCount = entry.getValue();
                altCounts.put(altKey, count != null ? count + newCount : newCount);
            }

            Set<BaseQualityKey> repeatedAltKeys = altCounts.entrySet().stream()
                    .filter(x -> x.getKey().Ref != x.getKey().Alt)
                    .filter(x -> x.getValue() > mConfig.QualityRecalibration.MaxAltCount).map(x -> x.getKey()).collect(Collectors.toSet());

            if(!repeatedAltKeys.isEmpty())
            {
                repeatedAlts.put(position, repeatedAltKeys);
            }
        }

        return repeatedAlts;
    }

    private static BaseQualityKey altKey(final BaseQualityKey key)
    {
        return new BaseQualityKey(key.Ref, key.Alt, null, (byte)0);
    }

    public void processRecord(@NotNull final SAMRecord record)
    {
        ++mReadCounter;
        CigarTraversal.traverseCigar(record, this);
    }

    @Override
    public void handleInsert(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int readIndex, final int refPos)
    {
        // Need to add one because indel is actually AFTER this by convention
        mIndelPositions.add(refPos + 1);
        handleAlignment(record, SINGLE, readIndex, refPos);
    }

    @Override
    public void handleDelete(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int readIndex, final int refPos)
    {
        mIndelPositions.add(refPos + 1);
        handleAlignment(record, SINGLE, readIndex, refPos);
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
        }
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
