package com.hartwig.hmftools.esvee.sam;

import java.util.Objects;
import java.util.stream.Stream;

import com.hartwig.hmftools.esvee.models.Record;

public class NormalisingSource implements SAMSource
{
    private final SAMSource mInner;
    private final RecordNormaliser mNormaliser;

    public NormalisingSource(final SAMSource inner, final RecordNormaliser normaliser)
    {
        mInner = inner;
        mNormaliser = normaliser;
    }

    @Override
    public Stream<Record> streamReadsContaining(final String chromosome, final int startPosition, final int endPosition)
    {
        return mInner.streamReadsContaining(chromosome, startPosition, endPosition)
                .map(mNormaliser::normalise)
                .filter(Objects::nonNull);
    }

    @Override
    public void close()
    {
        mInner.close();
    }
}
