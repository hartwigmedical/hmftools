package com.hartwig.hmftools.bamtools.slice;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bamops.BamToolName.fromPath;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.StringUtil;

public class SliceWriter
{
    private final SliceConfig mConfig;

    private int mRecordWriteCount;
    private final SAMFileWriter mBamWriter;
    private final BufferedWriter mReadWriter;
    private String mUnsortedOutputBam;
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

        boolean writeUnsortedBam = mConfig.UnsortedBam || mConfig.BamToolPath != null;
        mUnsortedOutputBam = mOutputBam.replaceAll("bam", "unsorted.bam");

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();

        String outputBamFile = mConfig.BamToolPath != null ? mUnsortedOutputBam : mOutputBam;

        if(writeUnsortedBam)
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        else
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(outputBamFile));
    }

    public synchronized void writeRead(final SAMRecord read) { doWriteRecord(read); }

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

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("ReadId").add("Chromosome").add("PosStart").add("PosEnd").add("Cigar");
            sj.add("InsertSize").add("MateChr").add("MatePosStart").add("MapQual").add("SuppData").add("Flags");
            sj.add("FirstInPair").add("ReadReversed").add("Proper").add("Unmapped").add("MateUnmapped").add("Supplementary").add("Duplicate");

            if(mConfig.WriteReadBases)
                sj.add("ReadSequence").add("BaseQual");

            writer.write(sj.toString());

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
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(record.getReadName());
            sj.add(record.getContig());
            sj.add(String.valueOf(record.getAlignmentStart()));
            sj.add(String.valueOf(record.getAlignmentEnd()));
            sj.add(record.getCigarString());

            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));

            sj.add(String.valueOf(abs(record.getInferredInsertSize())));
            sj.add(record.getMateReferenceName());
            sj.add(String.valueOf(record.getMateAlignmentStart()));
            sj.add(String.valueOf(record.getMappingQuality()));
            sj.add(suppData != null ? suppData.asDelimStr() : "N/A");
            sj.add(String.valueOf(record.getFlags()));

            sj.add(String.valueOf(record.getFirstOfPairFlag()));
            sj.add(String.valueOf(record.getReadNegativeStrandFlag()));
            sj.add(String.valueOf(record.getProperPairFlag()));
            sj.add(String.valueOf(record.getReadUnmappedFlag()));
            sj.add(String.valueOf(record.getMateUnmappedFlag()));
            sj.add(String.valueOf(record.getSupplementaryAlignmentFlag()));
            sj.add(String.valueOf(record.getDuplicateReadFlag()));

            if(mConfig.WriteReadBases)
            {
                sj.add(StringUtil.bytesToString(record.getReadBases()));
                sj.add(record.getBaseQualityString());
            }

            mReadWriter.write(sj.toString());
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
            BT_LOGGER.info("{} records written to BAM", mRecordWriteCount);

            mBamWriter.close();

            if(!mConfig.UnsortedBam && mConfig.BamToolPath != null)
            {
                BT_LOGGER.info("writing sorted BAM: {}", mOutputBam);

                BamToolName toolName = fromPath(mConfig.BamToolPath);

                boolean success = BamOperations.sortBam(toolName, mConfig.BamToolPath, mUnsortedOutputBam, mOutputBam, mConfig.Threads);

                if(success && toolName == BamToolName.SAMTOOLS)
                {
                    success = BamOperations.indexBam(toolName, mConfig.BamToolPath, mOutputBam, mConfig.Threads);
                }

                if(success)
                {
                    try
                    {
                        Files.deleteIfExists(Paths.get(mUnsortedOutputBam));
                    }
                    catch(IOException e)
                    {
                        BT_LOGGER.error("error deleting interim file: {}", e.toString());
                    }
                }
            }
        }

        closeBufferedWriter(mReadWriter);
    }
}
