package com.hartwig.hmftools.svassembly.sam;

import java.util.List;
import java.util.stream.Stream;

import com.hartwig.hmftools.svassembly.models.Record;

public class CompositeSAMSource implements SAMSource
{
    private final List<SAMSource> mSources;

    public CompositeSAMSource(final List<SAMSource> sources)
    {
        mSources = sources;
    }

    @Override
    public Stream<Record> unmappedReads()
    {
        return mSources.stream().flatMap(SAMSource::unmappedReads);
    }

    @Override
    public Stream<Record> streamReadsContaining(final String chromosome, final int startPosition, final int endPosition)
    {
        return mSources.stream().flatMap(source -> source.streamReadsContaining(chromosome, startPosition, endPosition));
    }

    @Override
    public void close()
    {
        mSources.forEach(SAMSource::close);
    }
}
