package com.hartwig.hmftools.sage.quality;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.sage.config.SageConfig;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class QualityRecalibration
{
    private final ExecutorService mExecutorService;
    private final IndexedFastaSequenceFile mRefGenome;
    private final SageConfig mConfig;

    public QualityRecalibration(final SageConfig config, final ExecutorService executorService, final IndexedFastaSequenceFile refGenome)
    {
        mExecutorService = executorService;
        mRefGenome = refGenome;
        mConfig = config;
    }

    @NotNull
    public CompletableFuture<List<QualityRecalibrationRecord>> qualityRecalibrationRecords(@NotNull final String bamFile)
    {
        final Map<QualityCounterKey, QualityCounter> map = new ConcurrentHashMap<>(Maps.newHashMap());

        final List<CompletableFuture<Void>> doneList = Lists.newArrayList();
        final List<GenomeRegion> regions =
                new QualityRecalibrationRegions(mRefGenome).regions(mConfig.BaseQualityRecalibration.SampleSize);

        for(GenomeRegion region : regions)
        {
            for(CompletableFuture<Collection<QualityCounter>> slice : submitAllRegions(bamFile, region))
            {
                final CompletableFuture<Void> done = slice.thenAccept(counts ->
                {
                    for(QualityCounter count : counts)
                    {
                        final QualityCounterKey key = withoutPosition(count);
                        map.computeIfAbsent(key, QualityCounter::new).increment(count.count());
                    }
                });
                doneList.add(done);
            }
        }

        return CompletableFuture.allOf(doneList.toArray(new CompletableFuture[0])).thenApply(aVoid ->
        {
            final List<QualityCounter> sortedList = Lists.newArrayList(map.values());
            Collections.sort(sortedList);
            return QualityRecalibrationFactory.create(sortedList);
        });
    }

    public CompletableFuture<Collection<QualityCounter>> addRegion(String bam, String contig, int start, int end)
    {
        BaseRegion bounds = new BaseRegion(contig, start, end);
        return CompletableFuture.supplyAsync(() -> new QualityCounterFactory(mConfig, bam, mRefGenome).regionCount(bounds), mExecutorService);
    }

    public List<CompletableFuture<Collection<QualityCounter>>> submitAllRegions(@NotNull final String bam,
            @NotNull final GenomeRegion region)
    {
        return submitAllRegions(bam, region.chromosome(), (int) region.start(), (int) region.end());
    }

    public List<CompletableFuture<Collection<QualityCounter>>> submitAllRegions(@NotNull final String bam, @NotNull final String contig,
            int minPosition, int maxPosition)
    {
        final List<CompletableFuture<Collection<QualityCounter>>> result = Lists.newArrayList();

        final int regionSliceSize = 100_000;
        for(int i = 0; ; i++)
        {
            int start = minPosition + i * regionSliceSize;
            int end = start + regionSliceSize - 1;

            if(end < minPosition)
            {
                continue;
            }

            result.add(addRegion(bam, contig, Math.max(start, minPosition), Math.min(end, maxPosition)));

            if(end >= maxPosition)
            {
                break;
            }
        }

        return result;
    }

    @NotNull
    private static QualityCounterKey withoutPosition(@NotNull final QualityCounter count)
    {
        return ImmutableQualityCounterKey.builder().from(count).position(0).build();
    }
}
