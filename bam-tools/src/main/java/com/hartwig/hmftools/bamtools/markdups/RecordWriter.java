package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.bamtools.markdups.ReadOutput.DUPLICATES;
import static com.hartwig.hmftools.bamtools.markdups.ReadOutput.MISMATCHES;
import static com.hartwig.hmftools.bamtools.markdups.ReadOutput.NONE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Set;

import com.google.common.collect.Sets;
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

    private final Set<SAMRecord> mReadsWritten;

    public RecordWriter(final MarkDupsConfig config)
    {
        mConfig = config;
        mRecordWriteCount = 0;
        mBamWriter = initialiseBam();
        mReadWriter = initialiseReadWriter();
        mReadsWritten = Sets.newHashSet();
    }

    private SAMFileWriter initialiseBam()
    {
        if(!mConfig.WriteBam)
            return null;

        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));

        mOutputBam = mConfig.OutputDir + mConfig.SampleId + ".mark_dups";

        if(mConfig.OutputId != null)
            mOutputBam += "." + mConfig.OutputId;

        mOutputBam += ".bam";

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();
        fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(mOutputBam));
    }

    public int recordWriteCount() { return mRecordWriteCount; }
    public Set<SAMRecord> readsWritten() { return mReadsWritten; }

    public synchronized void writeFragment(final Fragment fragment)
    {
        if(fragment.readsWritten())
        {
            BM_LOGGER.error("fragment({}) reads already written", fragment);
            return;
        }

        fragment.setReadWritten();
        fragment.reads().forEach(x -> writeRecord(x, fragment.status()));
    }

    private void writeRecord(final SAMRecord read, FragmentStatus fragmentStatus)
    {
        ++mRecordWriteCount;

        if(mConfig.runReadChecks())
        {
            if(mReadsWritten.contains(read))
            {
                BM_LOGGER.error("read({}) coords({}:{}-{}) already written",
                        read.getReadName(), read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd());
            }
            else
            {
                mReadsWritten.add(read);
            }
        }

        writeReadData(read, fragmentStatus);

        read.setDuplicateReadFlag(fragmentStatus == DUPLICATE); // overwrite any existing status

        if(mBamWriter != null)
            mBamWriter.addAlignment(read);
    }

    private BufferedWriter initialiseReadWriter()
    {
        if(mConfig.LogReadType == NONE)
            return null;

        try
        {
            String filename = mConfig.formFilename("reads");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("ReadId,Chromosome,PosStart,PosEnd,Cigar");
            writer.write(",InsertSize,MateChr,MatePosStart,Duplicate,CalcDuplicate,MapQual,SuppData,Flags");
            writer.write(",FirstInPair,ReadReversed,Proper,Unmapped,MateUnmapped,Supplementary,Secondary");

            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            BM_LOGGER.error(" failed to create read writer: {}", e.toString());
        }

        return null;
    }

    private void writeReadData(final SAMRecord read, FragmentStatus fragmentStatus)
    {
        if(mReadWriter == null)
            return;

        if(mConfig.LogReadType == DUPLICATES)
        {
            if(!read.getDuplicateReadFlag() && fragmentStatus == FragmentStatus.NONE)
                return;
        }
        else if(mConfig.LogReadType == MISMATCHES)
        {
            if(read.getDuplicateReadFlag() == (fragmentStatus == DUPLICATE))
                return;
        }

        try
        {
            mReadWriter.write(format("%s,%s,%d,%d,%s",
                    read.getReadName(), read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(), read.getCigar()));

            SupplementaryReadData suppData = SupplementaryReadData.from(read.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));

            mReadWriter.write(format(",%d,%s,%d,%s,%s,%d,%s,%d",
                    abs(read.getInferredInsertSize()), read.getMateReferenceName(), read.getMateAlignmentStart(), read.getDuplicateReadFlag(),
                    fragmentStatus, read.getMappingQuality(), suppData != null ? suppData.asCsv() : "N/A", read.getFlags()));

            mReadWriter.write(format(",%s,%s,%s,%s,%s,%s,%s",
                    read.getFirstOfPairFlag(), read.getReadNegativeStrandFlag(), read.getProperPairFlag(), read.getReadUnmappedFlag(),
                    read.getMateUnmappedFlag(), read.getSupplementaryAlignmentFlag(), read.isSecondaryAlignment()));

            mReadWriter.newLine();
        }
        catch(IOException e)
        {
            BM_LOGGER.error(" failed to write read data: {}", e.toString());
        }
    }

    public void close()
    {
        if(mBamWriter != null)
        {
            BM_LOGGER.info("{} records written to BAM: {}", mRecordWriteCount, mOutputBam);
            mBamWriter.close();
        }

        closeBufferedWriter(mReadWriter);
    }
}
