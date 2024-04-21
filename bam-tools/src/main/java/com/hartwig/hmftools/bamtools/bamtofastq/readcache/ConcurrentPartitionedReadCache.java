package com.hartwig.hmftools.bamtools.bamtofastq.readcache;

import static java.lang.Math.max;

import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Supplier;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.inject.Inject;
import com.hartwig.hmftools.bamtools.bamtofastq.BamToFastqConfig;

import org.apache.commons.lang3.NotImplementedException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

// TODO NEXT: TEST
public class ConcurrentPartitionedReadCache implements ReadCacheInterface
{
    private static final String UNMAPPED_KEY = "";

    private final BamToFastqConfig mConfig;
    private final Map<String, List<LockableReadCache>> mPartitionReadCachesByChromosome;

    private final AtomicInteger mSize;

    @Inject
    public ConcurrentPartitionedReadCache(final BamToFastqConfig config, final Supplier<LockableReadCache> lockableReadCacheSupplier)
    {
        mConfig = config;
        mPartitionReadCachesByChromosome = Maps.newHashMap();
        populatePartitionCachesByChromosome(lockableReadCacheSupplier);
        mSize = new AtomicInteger();
    }

    private void populatePartitionCachesByChromosome(final Supplier<LockableReadCache> lockableReadCacheSupplier)
    {
        mPartitionReadCachesByChromosome.put(UNMAPPED_KEY, Lists.newArrayList(lockableReadCacheSupplier.get()));
        for(SAMSequenceRecord refSequence : mConfig.RefGenome.refGenomeFile().getSequenceDictionary().getSequences())
        {
            List<LockableReadCache> chrReadCaches = Lists.newArrayList();
            mPartitionReadCachesByChromosome.put(refSequence.getContig(), chrReadCaches);
            int chromosomeLength = refSequence.getSequenceLength();
            for(int start = 1; start <= chromosomeLength; start += mConfig.PartitionSize)
            {
                chrReadCaches.add(lockableReadCacheSupplier.get());
            }
        }
    }

    private void add(final LockableReadCache lockableReadCache, final SAMRecord read)
    {
        try
        {
            lockableReadCache.lock();
            int currentSize = lockableReadCache.size();
            lockableReadCache.add(read);
            int newSize = lockableReadCache.size();
            mSize.getAndAdd(newSize - currentSize);
        }
        finally
        {
            lockableReadCache.unlock();
        }
    }

    @Override
    public void add(final SAMRecord read)
    {
        if(read.getReadUnmappedFlag() && read.getMateUnmappedFlag())
        {
            add(mPartitionReadCachesByChromosome.get(UNMAPPED_KEY).get(0), read);
            return;
        }

        String chromosome;
        int alignmentStart;
        if(read.getFirstOfPairFlag())
        {
            chromosome = read.getReferenceName();
            alignmentStart = read.getAlignmentStart();
        }
        else
        {
            chromosome = read.getMateReferenceName();
            alignmentStart = read.getMateAlignmentStart();
        }

        List<LockableReadCache> partitionReadCaches = mPartitionReadCachesByChromosome.get(chromosome);
        LockableReadCache lockableReadCache = partitionReadCaches.get(max(0, alignmentStart - 1) / mConfig.PartitionSize);
        add(lockableReadCache, read);
    }

    @Override
    public void flush()
    {
        mPartitionReadCachesByChromosome.values().stream().flatMap(List::stream).forEach(LockableReadCache::flush);
    }

    @Override
    public int size()
    {
        return mSize.get();
    }

    @Override
    public boolean isEmpty()
    {
        return mSize.get() == 0;
    }

    // TODO: We do not really need this here.
    @Override
    public void mergeStats(final ReadCacheInterface readCache)
    {
        throw new NotImplementedException("mergeStats is not implemented in ConcurrentPartitionedReadCache");
    }

    @Override
    public void logStats()
    {
        Stream<LockableReadCache> allLockableReadCaches = mPartitionReadCachesByChromosome.values().stream().flatMap(List::stream);
        LockableReadCache merged = null;
        for(Iterator<LockableReadCache> it = allLockableReadCaches.iterator(); it.hasNext(); )
        {
            LockableReadCache lockableReadCache = it.next();
            if(merged == null)
            {
                merged = lockableReadCache;
                continue;
            }

            merged.mergeStats(lockableReadCache);
        }

        merged.logStats();
    }
}
