package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.esvee.SvConstants.BAM_HEADER_SAMPLE_ID_TAG;

import java.io.File;
import java.util.Collections;
import java.util.HashSet;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.esvee.util.ThrowingConsumer;
import com.hartwig.hmftools.esvee.util.ThrowingFunction;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class DirectSAMSource implements SAMSource
{
    private final String mBamFileName;
    private final ThreadLocal<SamReader> mReader;
    private final Set<SamReader> mAllReaders = Collections.synchronizedSet(new HashSet<>());
    private final ReadRescue mReadRescue;

    @SuppressWarnings("unused")
    public DirectSAMSource(final File bamFile, final File referenceSequence)
    {
        this(bamFile, referenceSequence, null);
    }

    public DirectSAMSource(final File bamFile, final File referenceSequence, final String sampleId)
    {
        mBamFileName = bamFile.getAbsolutePath();

        mReader = ThreadLocal.withInitial(() ->
        {
            final SamReader reader = SamReaderFactory.makeDefault()
                    .referenceSequence(referenceSequence)
                    .open(bamFile);
            reader.getFileHeader().setAttribute("filename", mBamFileName);

            reader.getFileHeader().setAttribute(BAM_HEADER_SAMPLE_ID_TAG, sampleId);

            mAllReaders.add(reader);
            return reader;
        });

        mReadRescue = new ReadRescue(new RefGenomeSource(ThrowingFunction.rethrow((File file) -> new IndexedFastaSequenceFile(file)).apply(referenceSequence)));
    }

    @Override
    public Stream<Read> streamReadsContaining(final String chromosome, final int startPosition, final int endPosition)
    {
        final int sequenceIndex = mReader.get().getFileHeader().getSequenceIndex(chromosome);
        final QueryInterval interval = new QueryInterval(sequenceIndex, startPosition, endPosition);

        try (final Stream<Read> stream = mReader.get().queryOverlapping(new QueryInterval[] { interval }).stream()
                .filter(alignment -> !alignment.getDuplicateReadFlag())
                .map(Read::new)
                // .map(mReadRescue::rescueRead) // CHECK: not required
                .filter(Objects::nonNull)) {
            return stream.collect(Collectors.toList()).stream();
        }

        // private final RecordNormaliser mNormaliser;
        // then call mNormaliser::normalise() on each record
    }

    @Override
    public void close()
    {
        mAllReaders.forEach(ThrowingConsumer.rethrow(SamReader::close));
    }
}
