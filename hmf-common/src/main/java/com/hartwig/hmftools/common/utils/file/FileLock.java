package com.hartwig.hmftools.common.utils.file;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.RandomAccessFile;
import java.nio.channels.FileChannel;

import org.jetbrains.annotations.Nullable;

public final class FileLock implements AutoCloseable
{
    private final RandomAccessFile mFile;
    private final java.nio.channels.FileLock mLock;

    private FileLock(final RandomAccessFile file, final java.nio.channels.FileLock lock)
    {
        mFile = file;
        mLock = lock;
    }

    @Nullable
    public static FileLock create(final File file)
    {
        try
        {
            RandomAccessFile raf = new RandomAccessFile(file, "rw");
            FileChannel channel = raf.getChannel();
            java.nio.channels.FileLock lock = channel.lock(0L, Long.MAX_VALUE, false);
            return new FileLock(raf, lock);
        }
        catch(Exception e)
        {
            return null;
        }
    }

    public BufferedReader getBufferedReader() throws IOException
    {
        return new BufferedReader(new InputStreamReader(new FileInputStream(mFile.getFD())));
    }

    public BufferedWriter getBufferedWriter() throws IOException
    {
        return new BufferedWriter(new OutputStreamWriter(new FileOutputStream(mFile.getFD())));
    }

    public void clear() throws IOException
    {
        mFile.setLength(0);
    }

    @Override
    public void close() throws Exception
    {
        mLock.close();
        mFile.close();
    }
}
