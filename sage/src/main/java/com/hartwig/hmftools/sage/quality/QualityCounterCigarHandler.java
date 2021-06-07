package com.hartwig.hmftools.sage.quality;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.samtools.CigarHandler;
import com.hartwig.hmftools.common.samtools.CigarTraversal;
import com.hartwig.hmftools.sage.read.IndexedBases;
import com.hartwig.hmftools.sage.ref.RefSequence;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

class QualityCounterCigarHandler implements CigarHandler
{

    private static final CigarElement SINGLE = new CigarElement(1, CigarOperator.M);
    private static final byte N = (byte) 'N';

    private final IndexedBases refGenome;
    private final GenomeRegion bounds;
    private final int maxAltCount;

    private final Set<Integer> indelPositions = Sets.newHashSet();
    private final Map<QualityCounterKey, QualityCounter> qualityMap = Maps.newHashMap();

    public QualityCounterCigarHandler(final RefSequence refGenome, final GenomeRegion bounds, final int maxAltCount)
    {
        this.refGenome = refGenome.alignment();
        this.bounds = bounds;
        this.maxAltCount = maxAltCount;
    }

    public void processRecord(@NotNull final SAMRecord record)
    {
        CigarTraversal.traverseCigar(record, this);
    }

    @NotNull
    public Collection<QualityCounter> counts()
    {
        final Set<QualityCounterKey> altsToRemove = groupByAlt(qualityMap.values()).stream()
                .filter(x -> x.ref() != x.alt())
                .filter(x -> x.count() > maxAltCount)
                .map(QualityCounter::key)
                .collect(Collectors.toSet());

        final Set<QualityCounter> result = Sets.newHashSet();
        for(QualityCounter count : qualityMap.values())
        {
            final QualityCounterKey altKey = altKey(count);
            if(!indelPositions.contains(count.position()) && !altsToRemove.contains(altKey))
            {
                result.add(count);
            }
        }

        return result;
    }

    @Override
    public void handleInsert(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int readIndex, final int refPos)
    {
        // Need to add one because indel is actually AFTER this by convention
        indelPositions.add(refPos + 1);
        handleAlignment(record, SINGLE, readIndex, refPos);
    }

    @Override
    public void handleDelete(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int readIndex, final int refPos)
    {
        indelPositions.add(refPos + 1);
        handleAlignment(record, SINGLE, readIndex, refPos);
    }

    @Override
    public void handleAlignment(@NotNull final SAMRecord r, @NotNull final CigarElement e, final int startReadIndex, final int refPos)
    {
        for(int i = 0; i < e.getLength(); i++)
        {
            int readIndex = startReadIndex + i;
            int position = refPos + i;

            if(position > bounds.end())
            {
                return;
            }

            if(position < bounds.start())
            {
                continue;
            }

            byte ref = refGenome.base(position);
            byte alt = r.getReadBases()[readIndex];
            byte quality = r.getBaseQualities()[readIndex];
            byte[] trinucleotideContext = refGenome.trinucleotideContext(position);

            if(alt != N && isValid(trinucleotideContext))
            {
                final QualityCounterKey key = ImmutableQualityCounterKey.builder()
                        .ref(ref)
                        .alt(alt)
                        .qual(quality)
                        .position(position)
                        .trinucleotideContext(trinucleotideContext)
                        .build();
                qualityMap.computeIfAbsent(key, QualityCounter::new).increment();
            }
        }
    }

    private static boolean isValid(byte[] trinucleotideContext)
    {
        for(byte b : trinucleotideContext)
        {
            if(b == N)
            {
                return false;
            }
        }
        return trinucleotideContext.length == 3;
    }

    @NotNull
    private static List<QualityCounter> groupByAlt(final Collection<QualityCounter> quality)
    {
        final Map<QualityCounterKey, QualityCounter> map = Maps.newHashMap();

        for(QualityCounter count : quality)
        {
            final QualityCounterKey key = altKey(count);
            map.computeIfAbsent(key, QualityCounter::new).increment(count.count());
        }

        final List<QualityCounter> result = Lists.newArrayList(map.values());
        Collections.sort(result);

        return result;
    }

    @NotNull
    public static QualityCounterKey altKey(@NotNull final QualityCounter count)
    {
        return ImmutableQualityCounterKey.builder().from(count).qual((byte) 0).trinucleotideContext().build();
    }
}
