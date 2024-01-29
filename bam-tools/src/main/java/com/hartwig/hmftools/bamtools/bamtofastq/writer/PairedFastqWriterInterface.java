package com.hartwig.hmftools.bamtools.bamtofastq.writer;

import java.io.Closeable;
import java.io.IOException;
import java.util.function.BiConsumer;

import htsjdk.samtools.SAMRecord;

// TODO NEXT: TEST
public abstract class PairedFastqWriterInterface implements Closeable, BiConsumer<SAMRecord, SAMRecord>
{
    public abstract void writePairedFastqRecord(final SAMRecord read1, final SAMRecord read2);

    @Override
    public abstract void close() throws IOException;

    public abstract void logStats();

    @Override
    public final void accept(final SAMRecord read1, final SAMRecord read2)
    {
        writePairedFastqRecord(read1, read2);
    }
}
