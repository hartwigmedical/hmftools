package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;

import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.function.Consumer;

public class BufferedWriter<T> implements Consumer<T>, AutoCloseable
{
    private final BufferedWriterConsumer<T> mConsumer;
    private final Timestamp mTimestamp;
    private final List<T> mBuffer;
    private final int mBufferSize;
    private boolean mInitialised;

    public BufferedWriter(final BufferedWriterConsumer<T> consumer)
    {
        this(consumer, DB_BATCH_INSERT_SIZE);
    }

    public BufferedWriter(final BufferedWriterConsumer<T> consumer, int batchInsertSize)
    {
        mConsumer = consumer;
        mTimestamp = new Timestamp(new Date().getTime());
        mBufferSize = batchInsertSize;
        mBuffer = new ArrayList<>(batchInsertSize + 1);
    }

    public void initialise()
    {
        mInitialised = true;
        mConsumer.initialise();
    }

    @Override
    public void accept(final T entry)
    {
        if(!mInitialised)
        {
            mInitialised = true;
            mConsumer.initialise();
        }

        mBuffer.add(entry);
        if(mBuffer.size() >= mBufferSize)
        {
            writeBuffer();
        }
    }

    private void writeBuffer()
    {
        mConsumer.accept(mTimestamp, mBuffer);
        mBuffer.clear();
    }

    @Override
    public void close()
    {
        if(!mBuffer.isEmpty())
        {
            writeBuffer();
        }
    }
}
