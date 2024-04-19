package com.hartwig.hmftools.bamtools.slice;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.bam.BamOperations;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class SliceWriter
{
    private final SliceConfig mConfig;

    private int mRecordWriteCount;
    private final SAMFileWriter mBamWriter;
    private final BufferedWriter mReadWriter;
    private String mOutputBam;

    public SliceWriter(final SliceConfig config)
    {
        mConfig = config;
        mRecordWriteCount = 0;
        mBamWriter = config.WriteBam ? initialiseBam() : null;
        mReadWriter = config.WriteReads ? initialiseReadWriter() : null;
    }

    private SAMFileWriter initialiseBam()
    {
        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));

        mOutputBam = mConfig.formFilename(WriteType.BAM);

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();

        if(mConfig.UnsortedBam)
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        else
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(mOutputBam));
    }

    public synchronized void writeRead(final SAMRecord read) { doWriteRecord(read); }
    public synchronized void writeReads(final List<SAMRecord> reads) { reads.forEach(x -> doWriteRecord(x));}

    private void doWriteRecord(final SAMRecord read)
    {
        ++mRecordWriteCount;

        if(mBamWriter != null)
            mBamWriter.addAlignment(read);

        if(mReadWriter != null)
            writeReadData(read);
    }

    private BufferedWriter initialiseReadWriter()
    {
        try
        {
            String filename = mConfig.formFilename(WriteType.READS);
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("ReadId,Chromosome,PosStart,PosEnd,Cigar");
            writer.write(",InsertSize,MateChr,MatePosStart,MapQual,SuppData,Flags");
            writer.write(",FirstInPair,ReadReversed,Proper,Unmapped,MateUnmapped,Supplementary,Duplicate");

            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            BT_LOGGER.error(" failed to create read writer: {}", e.toString());
        }

        return null;
    }

    private void writeReadData(final SAMRecord record)
    {
        if(mReadWriter == null)
            return;

        try
        {
            mReadWriter.write(format("%s,%s,%d,%d,%s",
                    record.getReadName(), record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd(), record.getCigar()));

            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));

            mReadWriter.write(format(",%d,%s,%d,%d,%s,%d",
                    abs(record.getInferredInsertSize()), record.getMateReferenceName(), record.getMateAlignmentStart(), record.getMappingQuality(),
                    suppData != null ? suppData.asCsv() : "N/A", record.getFlags()));

            mReadWriter.write(format(",%s,%s,%s,%s,%s,%s,%s",
                    record.getFirstOfPairFlag(), record.getReadNegativeStrandFlag(), record.getProperPairFlag(), record.getReadUnmappedFlag(),
                    record.getMateUnmappedFlag(), record.getSupplementaryAlignmentFlag(), record.getDuplicateReadFlag()));

            mReadWriter.newLine();
        }
        catch(IOException e)
        {
            BT_LOGGER.error(" failed to write read data: {}", e.toString());
        }
    }

    public void close()
    {
        if(mBamWriter != null)
        {
            BT_LOGGER.info("{} records written to BAM: {}", mRecordWriteCount, mOutputBam);

            if(!mConfig.UnsortedBam)
            {
                // TODO
                // BT_LOGGER.debug("indexing sliced BAM");
                // BamOperations.indexBam(bamToolName(), bamToolPath(), finalBamFilename, mConfig.Threads)
            }

            mBamWriter.close();
        }

        closeBufferedWriter(mReadWriter);
    }
}
