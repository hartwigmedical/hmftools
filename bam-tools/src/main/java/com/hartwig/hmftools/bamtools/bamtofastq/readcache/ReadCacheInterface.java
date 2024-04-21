package com.hartwig.hmftools.bamtools.bamtofastq.readcache;

import java.util.function.Consumer;

import htsjdk.samtools.SAMRecord;

public interface ReadCacheInterface extends Consumer<SAMRecord>
{
    // TODO: Only needed for EvictingCache.
    void flush();

    int size();

    boolean isEmpty();

    void add(final SAMRecord read);

    // TODO: should this be here?
    void mergeStats(final ReadCacheInterface other);

    void logStats();

    @Override
    default void accept(final SAMRecord samRecord)
    {
        add(samRecord);
    }
}
