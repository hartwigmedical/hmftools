package com.hartwig.hmftools.common.bam;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.SAM_LOGGER;

import java.io.File;
import java.io.OutputStream;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.ProgressLoggerInterface;
import htsjdk.samtools.util.zip.DeflaterFactory;

public class FastBamWriter implements SAMFileWriter
{
    private final String mFilename;
    private final SAMFileHeader mHeader;

    private BinaryCodec mOutputBinaryCodec;
    private BAMRecordCodec mBamRecordCodec;
    private BlockCompressedOutputStream mBlockCompressedOutputStream;
    private boolean mClosed;

    // taken from HTS JDK
    private static final int HTSJDK_COMPRESSION_LEVEL = 5;
    private static final byte[] HTSJDK_BAM_MAGIC = "BAM\u0001".getBytes();
    private static final int BUFFER_SIZE = 131072;

    public FastBamWriter(final SAMFileHeader header, final String filename)
    {
        mHeader = header;
        mFilename = filename;

        try
        {
            DeflaterFactory deflaterFactory = BlockCompressedOutputStream.getDefaultDeflaterFactory();
            File outputFile = new File(filename);

            OutputStream os = IOUtil.maybeBufferOutputStream(Files.newOutputStream(outputFile.toPath()), BUFFER_SIZE);
            mBlockCompressedOutputStream = new BlockCompressedOutputStream(os, (Path) null, HTSJDK_COMPRESSION_LEVEL, deflaterFactory);
            mOutputBinaryCodec = new BinaryCodec(mBlockCompressedOutputStream);
            mOutputBinaryCodec.setOutputFileName(filename);

            writeHeader();

            mBamRecordCodec = new BAMRecordCodec(mHeader);
            mBamRecordCodec.setOutputStream(mOutputBinaryCodec.getOutputStream(), mFilename);
        }
        catch(Exception e)
        {
            SAM_LOGGER.error("failed to open BAM({}) for writing: {}", filename, e.toString());

            mOutputBinaryCodec = null;
            mBamRecordCodec = null;
            mBlockCompressedOutputStream = null;
        }

        mClosed = false;
    }

    public String filename() { return mFilename; }

    @Override
    public void addAlignment(SAMRecord record)
    {
        writeAlignment(record);
    }

    @Override
    public SAMFileHeader getFileHeader() { return mHeader; }

    @Override
    public void setProgressLogger(final ProgressLoggerInterface loggerInterface) {}

    @Override
    public void close()
    {
        if(mClosed)
            return;

        mOutputBinaryCodec.close();
        mClosed = true;
    }

    private void writeHeader()
    {
        StringWriter headerTextBuffer = new StringWriter();
        (new SAMTextHeaderCodec()).encode(headerTextBuffer, mHeader);

        mOutputBinaryCodec.writeBytes(HTSJDK_BAM_MAGIC);
        mOutputBinaryCodec.writeString(headerTextBuffer.toString(), true, false);
        mOutputBinaryCodec.writeInt(mHeader.getSequenceDictionary().size());

        Iterator sequenceIter = mHeader.getSequenceDictionary().getSequences().iterator();

        while(sequenceIter.hasNext())
        {
            SAMSequenceRecord sequenceRecord = (SAMSequenceRecord)sequenceIter.next();
            mOutputBinaryCodec.writeString(sequenceRecord.getSequenceName(), true, true);
            mOutputBinaryCodec.writeInt(sequenceRecord.getSequenceLength());
        }
    }

    private void writeAlignment(final SAMRecord record)
    {
        mBamRecordCodec.encode(record);
    }
}
