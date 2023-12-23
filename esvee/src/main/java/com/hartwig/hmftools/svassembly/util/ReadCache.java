package com.hartwig.hmftools.svassembly.util;

import java.lang.ref.SoftReference;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Supplier;
import java.util.stream.Stream;

import com.hartwig.hmftools.svassembly.models.Record;

import org.jetbrains.annotations.Nullable;

public class ReadCache
{
    private final Map<String, List<LRUNode>> mCachedByChromosome = new ConcurrentHashMap<>();

    @Nullable
    private LRUNode mHead;
    @Nullable
    private LRUNode mTail;

    private int mCachedAlignmentCount;
    private final int mMaxCachedAlignments;
    private final CacheLoader mLoader;

    public AtomicInteger CacheHits = new AtomicInteger();
    public AtomicInteger CacheMisses = new AtomicInteger();

    public ReadCache(final int maxCachedAlignments, final CacheLoader loader)
    {
        mMaxCachedAlignments = maxCachedAlignments;
        mLoader = loader;
    }

    public Stream<Record> read(final String chromosome, final int startPosition, final int endPosition)
    {
        final List<LRUNode> nodes = mCachedByChromosome.computeIfAbsent(chromosome, ignored -> new ArrayList<>());
        final CachedReadKey key;
        final CompletableFuture<CachedReadValue> future;
        synchronized (this)
        {
            for(int i = 0; i < nodes.size(); i++)
            {
                final LRUNode node = nodes.get(i);
                if(node.Key.containsCompletely(chromosome, startPosition, endPosition))
                {
                    @Nullable
                    final Stream<Record> reads = node.Value.get().extract(chromosome, startPosition, endPosition);
                    if(reads != null)
                    {
                        CacheHits.incrementAndGet();
                        touch(node);
                        return reads;
                    }

                    // This node has expired
                    removeFully(node, i--);
                }
            }

            // No node caches what we need. We'll need to load
            key = mLoader.determineCacheRange(chromosome, startPosition, endPosition);
            future = new CompletableFuture<>();
            final var newNode = new LRUNode(key, ThrowingSupplier.rethrow(future::get));
            nodes.add(newNode);
            touch(newNode);
        }

        CacheMisses.incrementAndGet();

        while(true)
        {
            final CachedReadValue value;
            @Nullable
            final List<Record> records;
            try
            {
                value = mLoader.load(key);
                records = value.Records.get(); // Hold a reference to ensure extract works.
            }
            catch(final Throwable throwable)
            {
                future.completeExceptionally(throwable);
                throw throwable;
            }
            if(records == null)
                continue;

            future.complete(value);
            mCachedAlignmentCount += value.Size;
            while(mCachedAlignmentCount > mMaxCachedAlignments)
                evictOldest();
            return value.extract(chromosome, startPosition, endPosition);
        }
    }

    private synchronized void touch(final LRUNode node)
    {
        remove(node);

        assert (mTail == null) == (mHead == null);
        if(mHead == null)
            mHead = mTail = node;
        else
        {
            mTail.Next = node;
            node.Previous = mTail;
            mTail = node;
        }
    }

    private synchronized void evictOldest()
    {
        while(mHead != null)
        {
            final boolean wasResident = removeFully(mHead, -1);
            if(wasResident)
                break;
        }
    }

    /**
     * Assumes caller holding a lock, does not update mCachedReads
     */
    private void remove(final LRUNode node)
    {
        final var prev = node.Previous;
        final var next = node.Next;

        if(prev != null)
            prev.Next = next;
        if(next != null)
            next.Previous = prev;

        node.Previous = node.Next = null;

        if(mHead == node)
            mHead = next;
        if(mTail == node)
            mTail = prev;
    }

    private boolean removeFully(final LRUNode node, final int listIndex)
    {
        final var toRemove = mHead;
        if(toRemove == null)
            return false;

        remove(toRemove);

        final List<LRUNode> list = mCachedByChromosome.get(node.Key.Chromosome);
        if(listIndex == -1)
        {
            for(int i = 0; i < list.size(); i++)
                if(list.get(i) == node)
                {
                    list.set(i, list.get(list.size() - 1));
                    list.remove(list.size() - 1);
                    break;
                }
        }
        else if(list.size() == 1)
            list.clear();
        else
        {
            list.set(listIndex, list.get(list.size() - 1));
            list.remove(list.size() - 1);
        }

        final CachedReadValue value = toRemove.Value.get();

        mCachedAlignmentCount -= value.Size;
        return value.Records.get() != null;
    }

    public interface CacheLoader
    {
        CachedReadKey determineCacheRange(final String requestedChromosome, final int requestedStartPosition,
                final int requestedEndPosition);

        CachedReadValue load(final CachedReadKey key);
    }

    private static class LRUNode
    {
        public final CachedReadKey Key;
        public final Supplier<CachedReadValue> Value;

        @Nullable
        public LRUNode Previous, Next;

        private LRUNode(final CachedReadKey key, final Supplier<CachedReadValue> value)
        {
            Key = key;
            Value = value;
        }
    }

    public static class CachedReadKey
    {
        public final String Chromosome;
        public final int StartPosition, EndPosition;

        public CachedReadKey(final String chromosome, final int startPosition, final int endPosition)
        {
            Chromosome = chromosome;
            StartPosition = startPosition;
            EndPosition = endPosition;
        }

        public boolean containsCompletely(final String chromosome, final int startPosition, final int endPosition)
        {
            return Chromosome.equals(chromosome) && StartPosition <= startPosition && EndPosition >= endPosition;
        }
    }

    public static class CachedReadValue
    {
        public final int Size;
        public final SoftReference<List<Record>> Records;

        public CachedReadValue(final List<Record> records)
        {
            Size = records.size();
            Records = new SoftReference<>(records);
        }

        @Nullable
        public Stream<Record> extract(final String chromosome, final int startPosition, final int endPosition)
        {
            @Nullable
            final List<Record> records = Records.get();
            if(records == null)
                return null;

            return records.stream()
                    .filter(r -> overlaps(r, chromosome, startPosition, endPosition));
        }

        private boolean overlaps(final Record record, final String chromosome, final int startPosition, final int endPosition)
        {
            if(!record.isUnmapped())
            {
                return record.getChromosome().equals(chromosome)
                        && record.getAlignmentStart() <= endPosition && record.getAlignmentEnd() >= startPosition;
            }
            else
            {
                if(!record.isPairedRead())
                    return false;

                final int mateStart = record.getMateAlignmentStart();
                final int mateEnd = mateStart + (record.getLength() * 2);
                return Objects.requireNonNull(record.getMateChromosome()).equals(chromosome)
                        && mateStart <= endPosition && mateEnd >= startPosition;
            }
        }
    }
}
