package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.lang.ref.SoftReference;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Supplier;
import java.util.stream.Stream;

import org.jetbrains.annotations.Nullable;

public class ReadCache
{
    private final Map<String,List<LRUNode>> mCachedByChromosome;

    public static final int MAX_CACHE_SIZE = 3_000_000;
    public static final int CACHE_QUERY_BUFFER = 2000;

    @Nullable
    private LRUNode mHead;
    @Nullable
    private LRUNode mTail;

    private int mCachedAlignmentCount;
    private final int mMaxCachedAlignments;
    private final CacheLoader mLoader;

    private final AtomicInteger mCacheHits;
    private final AtomicInteger mCacheMisses;

    public interface CacheLoader
    {
        CachedReadKey determineCacheRange(final String requestedChromosome, final int positionStart, final int positionEnd);

        CachedReadValue load(final CachedReadKey key);
    }

    public ReadCache(final int maxCachedAlignments, final CacheLoader loader)
    {
        mMaxCachedAlignments = maxCachedAlignments;
        mLoader = loader;
        mCachedByChromosome = new ConcurrentHashMap<>();

        mCacheHits = new AtomicInteger();
        mCacheMisses = new AtomicInteger();
    }

    public int totalHits() { return mCacheHits.get() + mCacheMisses.get(); }

    public void logStats()
    {
        final int hits = mCacheHits.get();
        final int misses = mCacheMisses.get();

        SV_LOGGER.info("closing cache: cache hits({}) misses({}) hitRate({}%)",
                hits, misses, String.format("%.2f", (100.0f * hits) / (hits + misses)));
    }

    public Stream<Read> read(final String chromosome, final int positionStart, final int positionEnd)
    {
        final List<LRUNode> nodes = mCachedByChromosome.computeIfAbsent(chromosome, ignored -> new ArrayList<>());
        final CachedReadKey key;
        final CompletableFuture<CachedReadValue> future;
        synchronized (this)
        {
            for(int i = 0; i < nodes.size(); i++)
            {
                final LRUNode node = nodes.get(i);
                if(node.Key.containsCompletely(chromosome, positionStart, positionEnd))
                {
                    @Nullable
                    final Stream<Read> reads = node.Value.get().extract(chromosome, positionStart, positionEnd);
                    if(reads != null)
                    {
                        mCacheHits.incrementAndGet();
                        touch(node);
                        return reads;
                    }

                    // This node has expired
                    removeFully(node, i--);
                }
            }

            // No node caches what we need. We'll need to load
            key = mLoader.determineCacheRange(chromosome, positionStart, positionEnd);
            future = new CompletableFuture<>();
            final var newNode = new LRUNode(key, ThrowingSupplier.rethrow(future::get));
            nodes.add(newNode);
            touch(newNode);
        }

        mCacheMisses.incrementAndGet();

        while(true)
        {
            final CachedReadValue value;
            @Nullable
            final List<Read> reads;
            try
            {
                value = mLoader.load(key);
                reads = value.Records.get(); // Hold a reference to ensure extract works.
            }
            catch(final Throwable throwable)
            {
                future.completeExceptionally(throwable);
                throw throwable;
            }
            if(reads == null)
                continue;

            future.complete(value);
            mCachedAlignmentCount += value.Size;
            while(mCachedAlignmentCount > mMaxCachedAlignments)
                evictOldest();
            return value.extract(chromosome, positionStart, positionEnd);
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

    // assumes caller holding a lock, does not update mCachedReads
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
            {
                if(list.get(i) == node)
                {
                    list.set(i, list.get(list.size() - 1));
                    list.remove(list.size() - 1);
                    break;
                }
            }
        }
        else if(list.size() == 1)
        {
            list.clear();
        }
        else
        {
            list.set(listIndex, list.get(list.size() - 1));
            list.remove(list.size() - 1);
        }

        final CachedReadValue value = toRemove.Value.get();

        mCachedAlignmentCount -= value.Size;
        return value.Records.get() != null;
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
        public final int PositionStart;
        public final int PositionEnd;

        public CachedReadKey(final String chromosome, final int positionStart, final int positionEnd)
        {
            Chromosome = chromosome;
            PositionStart = positionStart;
            PositionEnd = positionEnd;
        }

        public boolean containsCompletely(final String chromosome, final int startPosition, final int endPosition)
        {
            return positionsWithin(startPosition, endPosition, PositionStart, PositionEnd);
        }
    }

    public static class CachedReadValue
    {
        public final int Size;
        public final SoftReference<List<Read>> Records;

        public CachedReadValue(final List<Read> reads)
        {
            Size = reads.size();
            Records = new SoftReference<>(reads);
        }

        public Stream<Read> extract(final String chromosome, final int positionStart, final int positionEnd)
        {
            final List<Read> reads = Records.get();
            if(reads == null)
                return null;

            return reads.stream().filter(x -> overlaps(x, chromosome, positionStart, positionEnd));
        }

        private boolean overlaps(final Read read, final String chromosome, final int positionStart, final int positionEnd)
        {
            // filter should match a normal slice which the existing logic did not
            if(read.isUnmapped())
            {
                // CHECK: use the mate's position or actually are these the same anyway?
                int matePosStart = read.mateAlignmentStart();
                int matePosEnd = matePosStart + read.getLength() * 2; // CHECK: consider using actual mate end using cigar
                return positionsOverlap(matePosStart, matePosEnd, positionStart, positionEnd);
            }
            else
            {
                return positionsOverlap(read.getAlignmentStart(), read.getAlignmentEnd(), positionStart, positionEnd);
            }

            /*
            if(!read.isUnmapped())
            {
                return read.getChromosome().equals(chromosome)
                        && read.getAlignmentStart() <= positionEnd && read.getAlignmentEnd() >= positionStart;
            }
            else
            {
                if(!read.isPairedRead())
                    return false;

                final int mateStart = read.getMateAlignmentStart();
                final int mateEnd = mateStart + (read.getLength() * 2);
                return Objects.requireNonNull(read.getMateChromosome()).equals(chromosome)
                        && mateStart <= positionEnd && mateEnd >= positionStart;
            }
            */
        }
    }

    private interface ThrowingSupplier<RET> extends Supplier<RET>
    {
        static <RET> Supplier<RET> rethrow(final ThrowingSupplier<RET> supplier)
        {
            return supplier;
        }

        RET getThrowing() throws Exception;

        @Override
        default RET get()
        {
            try
            {
                return getThrowing();
            }
            catch(final Exception e)
            {
                throw new RuntimeException(e);
            }
        }
    }
}
