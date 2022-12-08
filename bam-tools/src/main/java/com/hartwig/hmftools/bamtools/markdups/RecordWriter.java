package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class RecordWriter
{
    private final MarkDupsConfig mConfig;

    private int mRecordWriteCount;
    private final SAMFileWriter mBamWriter;
    private final BufferedWriter mReadWriter;
    private String mOutputBam;

    public RecordWriter(final MarkDupsConfig config)
    {
        mConfig = config;
        mRecordWriteCount = 0;
        mBamWriter = initialiseBam();
        mReadWriter = initialiseReadWriter();
    }

    private SAMFileWriter initialiseBam()
    {
        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));

        mOutputBam = mConfig.OutputDir + mConfig.SampleId + ".slice";

        if(mConfig.OutputId != null)
            mOutputBam += "." + mConfig.OutputId;

        mOutputBam += ".bam";

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();
        fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(mOutputBam));
    }

    public synchronized void writeRecord(final SAMRecord record, boolean isDuplicate)
    {
        ++mRecordWriteCount;

        writeReadData(record, isDuplicate);

        record.setDuplicateReadFlag(isDuplicate); // overwrite any existing status

        if(mBamWriter != null)
            mBamWriter.addAlignment(record);
    }

    private BufferedWriter initialiseReadWriter()
    {
        try
        {
            String filename = mConfig.formFilename("reads");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("ReadId,Chromosome,PosStart,PosEnd,Cigar");
            writer.write(",InsertSize,MateChr,MatePosStart,MapQual,SuppData,Flags");
            writer.write(",FirstInPair,ReadReversed,Proper,Unmapped,MateUnmapped,Supplementary,Duplicate,CalcDuplicate");

            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            BM_LOGGER.error(" failed to create read writer: {}", e.toString());
        }

        return null;
    }

    private void writeReadData(final SAMRecord record, boolean isDuplicate)
    {
        if(mReadWriter == null)
            return;

        if(!record.getDuplicateReadFlag() && !isDuplicate)
            return;

        try
        {
            mReadWriter.write(format("%s,%s,%d,%d,%s",
                    record.getReadName(), record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd(), record.getCigar()));

            SupplementaryReadData suppData = SupplementaryReadData.from(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));

            mReadWriter.write(format(",%d,%s,%d,%d,%s,%d",
                    abs(record.getInferredInsertSize()), record.getMateReferenceName(), record.getMateAlignmentStart(), record.getMappingQuality(),
                    suppData != null ? suppData.asCsv() : "N/A", record.getFlags()));

            mReadWriter.write(format(",%s,%s,%s,%s,%s,%s,%s,%s",
                    record.getFirstOfPairFlag(), record.getReadNegativeStrandFlag(), record.getProperPairFlag(), record.getReadUnmappedFlag(),
                    record.getMateUnmappedFlag(), record.getSupplementaryAlignmentFlag(), record.getDuplicateReadFlag(), isDuplicate));

            mReadWriter.newLine();
        }
        catch(IOException e)
        {
            BM_LOGGER.error(" failed to write read data: {}", e.toString());
        }
    }

    public void close()
    {
        // write any non-duplicate supplementaries
        // BM_LOGGER.debug("final duplicate readIds({}) cached supplementaries({})",
        //        mDuplicateReadIds.size(), mSupplementaryReads.size());

        /*
        mSupplementaryReads.values().forEach(x -> writeRecord(x));
        mSupplementaryReads.clear();
        mDuplicateReadIds.clear();
        */

        if(mBamWriter != null)
        {
            BM_LOGGER.info("{} records written to BAM: {}", mRecordWriteCount, mOutputBam);
            mBamWriter.close();
        }

        closeBufferedWriter(mReadWriter);
    }
}
