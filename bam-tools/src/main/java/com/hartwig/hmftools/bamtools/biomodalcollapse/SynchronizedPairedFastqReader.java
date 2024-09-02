package com.hartwig.hmftools.bamtools.biomodalcollapse;

import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalCollapseUtil.nextFastqRecord;

import java.io.BufferedReader;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.fastq.FastqRecord;

public class SynchronizedPairedFastqReader
{
    private final BufferedReader mFastq1Reader;
    private final BufferedReader mFastq2Reader;
    private final int mMaxPairedFastqRead;

    private int mPairedFastqReadCount;
    private boolean mFinished;

    public SynchronizedPairedFastqReader(final BufferedReader fastq1Reader, final BufferedReader fastq2Reader, int maxPairedFastqRead)
    {
        mFastq1Reader = fastq1Reader;
        mFastq2Reader = fastq2Reader;
        mMaxPairedFastqRead = maxPairedFastqRead;

        mPairedFastqReadCount = 0;
        mFinished = false;
    }

    public SynchronizedPairedFastqReader(final BufferedReader fastq1Reader, final BufferedReader fastq2Reader)
    {
        this(fastq1Reader, fastq2Reader, -1);
    }

    @Nullable
    public synchronized Pair<FastqRecord, FastqRecord> getNext()
    {
        if(mFinished)
        {
            return null;
        }

        FastqRecord fastq1 = nextFastqRecord(mFastq1Reader);
        FastqRecord fastq2 = nextFastqRecord(mFastq2Reader);

        if(fastq1 == null && fastq2 == null)
        {
            mFinished = true;
            return null;
        }

        if(fastq1 == null)
        {
            throw new RuntimeException("Fastq1 reader ran out of records before Fastq2");
        }

        if(fastq2 == null)
        {
            throw new RuntimeException("Fastq2 reader ran out of records before Fastq1");
        }

        mPairedFastqReadCount++;
        if(mMaxPairedFastqRead > 0 && mPairedFastqReadCount >= mMaxPairedFastqRead)
        {
            mFinished = true;
        }

        return Pair.of(fastq1, fastq2);
    }
}
