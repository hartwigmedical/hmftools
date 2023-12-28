package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.util.List;
import java.util.stream.Stream;

import com.hartwig.hmftools.esvee.util.ReadCache;

public class CachingSAMSource implements SAMSource
{
    private final SAMSource mInner;

    private final ReadCache mCache;

    public CachingSAMSource(final SAMSource inner)
    {
        mInner = inner;

        mCache = new ReadCache(3_000_000, new ReadCache.CacheLoader()
        {
            @Override
            public ReadCache.CachedReadKey determineCacheRange(final String chromosome, final int startPosition, final int endPosition)
            {
                final int adjustedStart = Math.max(1, startPosition - 2_000);
                final int adjustedEnd = endPosition + 2_000;

                return new ReadCache.CachedReadKey(chromosome, adjustedStart, adjustedEnd);
            }

            @Override
            public ReadCache.CachedReadValue load(final ReadCache.CachedReadKey key)
            {
                final List<Read> alignments = mInner.findReadsContaining(key.Chromosome, key.StartPosition, key.EndPosition);
                return new ReadCache.CachedReadValue(alignments);
            }
        });
    }

    @Override
    public Stream<Read> streamReadsContaining(final String chromosome, final int startPosition, final int endPosition)
    {
        try
        {
            return mCache.read(chromosome, startPosition, endPosition);
        }
        finally
        {
            final int hits = mCache.CacheHits.get();
            final int misses = mCache.CacheMisses.get();
            if((hits + misses) % 10_000 == 0)
            {
                SV_LOGGER.debug("Cache Hits: {}, Cache Misses: {}, Hit Rate: {}%",
                        hits, misses, String.format("%.2f", (100.0f * hits) / (hits + misses)));
            }
        }
    }

    @Override
    public void close()
    {
        final int hits = mCache.CacheHits.get();
        final int misses = mCache.CacheMisses.get();

        SV_LOGGER.info("closing cache: cache hits({}) misses({}) hitRate({}%)",
                hits, misses, String.format("%.2f", (100.0f * hits) / (hits + misses)));

        mInner.close();
    }
}
