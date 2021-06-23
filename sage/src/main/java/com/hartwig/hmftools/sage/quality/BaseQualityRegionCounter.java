package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
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
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
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

class BaseQualityRegionCounter implements CigarHandler, Callable
{
    private final String mBamFile;
    private final BaseRegion mRegion;
    private final ReferenceSequenceFile mRefGenome;
    private final IndexedBases mIndexedBases;
    private final SageConfig mConfig;

    private final Set<Integer> mIndelPositions = Sets.newHashSet();

    private final Set<QualityCounter> mQualityCounts;

    private final Map<QualityCounterKey,QualityCounter> mQualityMap = Maps.newHashMap();

    private static final CigarElement SINGLE = new CigarElement(1, CigarOperator.M);
    private static final byte N = (byte) 'N';

    public BaseQualityRegionCounter(
            final SageConfig config, final String bamFile, final ReferenceSequenceFile refGenome, final BaseRegion region)
    {
        mConfig = config;
        mBamFile = bamFile;

        mRegion = region;
        mRefGenome = refGenome;

        final RefSequence refSequence = new RefSequence(mRegion, mRefGenome);
        mIndexedBases = refSequence.alignment();

        mQualityCounts = Sets.newHashSet();
    }

    @Override
    public Long call()
    {
        produceRegionCounts();
        return (long)0;
    }

    public Collection<QualityCounter> getQualityCounts() { return mQualityCounts; }

    public void produceRegionCounts()
    {
        SG_LOGGER.trace("processing BQR region {}", mRegion);

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

        // remove locations where the alt count exceeds the configured limit
        List<QualityCounter> countsByAlt = groupByAlt(mQualityMap.values());

        final Set<QualityCounterKey> altsToRemove = countsByAlt.stream()
                .filter(x -> x.ref() != x.alt())
                .filter(x -> x.count() > mConfig.QualityRecalibration.MaxAltCount)
                .map(x -> x.Key)
                .collect(Collectors.toSet());

        for(QualityCounter count : mQualityMap.values())
        {
            final QualityCounterKey altKey = altKey(count);
            if(!mIndelPositions.contains(count.position()) && !altsToRemove.contains(altKey))
            {
                mQualityCounts.add(count);
            }
        }
    }

    public void processRecord(@NotNull final SAMRecord record)
    {
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

            if(alt != N && isValid(trinucleotideContext))
            {
                QualityCounterKey key = new QualityCounterKey(ref, alt, quality, position, trinucleotideContext);
                QualityCounter counter = mQualityMap.get(key);

                if(counter == null)
                {
                    counter = new QualityCounter(key);
                    mQualityMap.put(key, counter);
                }
                else
                {
                    if(key.compareTo(counter.Key) != 0)
                    {
                        SG_LOGGER.info("invalid: key({}) vs key2({})", key.toString(), counter.Key.toString());
                    }
                }

                counter.increment();

                // mQualityMap.computeIfAbsent(key, QualityCounter::new).increment();
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

    private static List<QualityCounter> groupByAlt(final Collection<QualityCounter> quality)
    {
        // determine counts but only comparing by ref, alt and position, not trinucleotide context
        final Map<QualityCounterKey,QualityCounter> map = Maps.newHashMap();

        for(QualityCounter count : quality)
        {
            QualityCounterKey altKey = altKey(count);
            QualityCounter counter = map.get(altKey);

            if(counter == null)
            {
                counter = new QualityCounter(altKey(count));
                map.put(altKey, counter);
            }

            counter.increment(count.count());
        }

        final List<QualityCounter> result = Lists.newArrayList(map.values());
        Collections.sort(result);

        return result;
    }

    private static QualityCounterKey altKey(final QualityCounter count)
    {
        return new QualityCounterKey(count.ref(), count.alt(), (byte)0, count.position(), null);
    }

}
