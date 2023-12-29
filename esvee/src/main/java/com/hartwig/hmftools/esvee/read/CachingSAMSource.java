package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.util.List;
import java.util.stream.Stream;

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
            public ReadCache.CachedReadKey determineCacheRange(final String chromosome, final int positionStart, final int positionEnd)
            {
                final int adjustedStart = Math.max(1, positionStart - 2_000);
                final int adjustedEnd = positionEnd + 2_000;

                return new ReadCache.CachedReadKey(chromosome, adjustedStart, adjustedEnd);
            }

            @Override
            public ReadCache.CachedReadValue load(final ReadCache.CachedReadKey key)
            {
                final List<Read> alignments = mInner.findReadsContaining(key.Chromosome, key.PositionStart, key.PositionEnd);
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
            if(mCache.totalHits() % 10_000 == 0)
            {
                mCache.logStats();
            }
        }
    }

    @Override
    public void close()
    {
        mCache.logStats();
        mInner.close();
    }
}
