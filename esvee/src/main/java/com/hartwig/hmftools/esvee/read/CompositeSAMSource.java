package com.hartwig.hmftools.esvee.read;

import java.util.List;
import java.util.stream.Stream;

public class CompositeSAMSource implements SAMSource
{
    private final List<SAMSource> mSources;

    public CompositeSAMSource(final List<SAMSource> sources)
    {
        mSources = sources;
    }

    @Override
    public Stream<Read> streamReadsContaining(final String chromosome, final int startPosition, final int endPosition)
    {
        return mSources.stream().flatMap(source -> source.streamReadsContaining(chromosome, startPosition, endPosition));
    }

    @Override
    public void close()
    {
        mSources.forEach(SAMSource::close);
    }
}