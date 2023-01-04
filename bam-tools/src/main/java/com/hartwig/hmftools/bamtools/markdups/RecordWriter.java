package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.SUPPLEMENTARY;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.UNSET;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.readToString;
import static com.hartwig.hmftools.bamtools.markdups.ReadOutput.DUPLICATES;
import static com.hartwig.hmftools.bamtools.markdups.ReadOutput.MISMATCHES;
import static com.hartwig.hmftools.bamtools.markdups.ReadOutput.NONE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
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

    private final BamWriter mBamWriter;
    private final BamWriter mCandidateBamWriter;
    private final BamWriter mSupplementaryBamWriter;
    private final BufferedWriter mReadWriter;
    private final BufferedWriter mResolvedReadWriter;

    private final Set<SAMRecord> mReadsWritten; // debug only

    private class BamWriter
    {
        private final String mFilename;
        private final SAMFileWriter mBamWriter;
        private int mWriteCount;

        public BamWriter(final String filename)
        {
            mFilename = filename;
            mBamWriter = filename != null ? initialiseBam(filename) : null;
            mWriteCount = 0;
        }

        public void writeRecord(final SAMRecord read)
        {
            ++mWriteCount;

            if(mBamWriter != null)
                mBamWriter.addAlignment(read);
        }

        public String filename() { return mFilename; }
        public int writeCount() { return mWriteCount; }

        public void close()
        {
            if(mBamWriter != null)
                mBamWriter.close();
        }
    }

    public RecordWriter(final MarkDupsConfig config)
    {
        mConfig = config;

        if(mConfig.WriteBam)
        {
            String bamFilename = formBamFilename("mark_dups");
            BM_LOGGER.info("writing new BAM file: {}", bamFilename);
            mBamWriter = new BamWriter(bamFilename);
        }
        else
        {
            mBamWriter = new BamWriter(null);
        }

        if(mConfig.UseInterimFiles)
        {
            mCandidateBamWriter = new BamWriter(formBamFilename("candidate"));
            mSupplementaryBamWriter = new BamWriter(formBamFilename("supplementary"));
            mResolvedReadWriter = initialiseResolvedReadWriter();
        }
        else
        {
            mCandidateBamWriter = null;
            mSupplementaryBamWriter = null;
            mResolvedReadWriter = null;
        }

        mReadWriter = initialiseReadWriter();
        mReadsWritten = Sets.newHashSet();
    }

    private SAMFileWriter initialiseBam(final String filename)
    {
        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();
        fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(filename));
    }

    private String formBamFilename(final String type)
    {
        String filename = mConfig.OutputDir + mConfig.SampleId + "." + type;

        if(mConfig.OutputId != null)
            filename += "." + mConfig.OutputId;

        filename += ".bam";
        return filename;
    }

    public int recordWriteCount() { return mBamWriter.writeCount(); }

    public synchronized void writeFragments(final List<Fragment> fragments) { fragments.forEach(x -> doWriteFragment(x)); }
    public synchronized void writeFragment(final Fragment fragment) { doWriteFragment(fragment); }

    private void doWriteFragment(final Fragment fragment)
    {
        if(fragment.readsWritten())
        {
            BM_LOGGER.error("fragment({}) reads already written", fragment);
            return;
        }

        fragment.setReadWritten();
        fragment.reads().forEach(x -> writeRead(x, fragment));
    }

    public synchronized void writeCachedFragment(final Fragment fragment) { doWriteCachedFragment(fragment); }

    private void doWriteCachedFragment(final Fragment fragment)
    {
        if(!mConfig.UseInterimFiles)
            return;

        if(fragment.status() == SUPPLEMENTARY)
        {
            fragment.reads().forEach(x -> mSupplementaryBamWriter.writeRecord(x));
        }
        else
        {
            fragment.reads().forEach(x -> mCandidateBamWriter.writeRecord(x));
        }
    }

    private void writeRead(final SAMRecord read, final Fragment fragment)
    {
        if(mConfig.runReadChecks())
        {
            if(mReadsWritten.contains(read))
            {
                BM_LOGGER.error("read({}) already written", readToString(read));
            }
            else
            {
                mReadsWritten.add(read);
            }
        }

        writeReadData(read, fragment);

        read.setDuplicateReadFlag(fragment.status() == DUPLICATE); // overwrite any existing status

        mBamWriter.writeRecord(read);
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
            writer.write(",InsertSize,MateChr,MatePosStart,Duplicate,CalcDuplicate,MateCigar,Coords,AvgBaseQual,MapQual,SuppData");
            writer.write(",Flags,FirstInPair,ReadReversed,Proper,Unmapped,MateUnmapped,Supplementary,Secondary");

            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            BM_LOGGER.error(" failed to create read writer: {}", e.toString());
        }

        return null;
    }

    private void writeReadData(final SAMRecord read, final Fragment fragment)
    {
        if(mReadWriter == null)
            return;

        if(mConfig.LogReadType == DUPLICATES)
        {
            if(!read.getDuplicateReadFlag() && !fragment.status().isDuplicate())
                return;
        }
        else if(mConfig.LogReadType == MISMATCHES)
        {
            if(fragment.status() != UNSET)
            {
                if(read.getDuplicateReadFlag() == (fragment.status() == DUPLICATE))
                    return;
            }
        }

        try
        {
            mReadWriter.write(format("%s,%s,%d,%d,%s",
                    read.getReadName(), read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(), read.getCigar()));

            SupplementaryReadData suppData = SupplementaryReadData.from(read.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));

            mReadWriter.write(format(",%d,%s,%d,%s,%s,%s,%s,%.2f,%d,%s,%d",
                    abs(read.getInferredInsertSize()), read.getMateReferenceName(), read.getMateAlignmentStart(),
                    read.getDuplicateReadFlag(), fragment.status(), read.hasAttribute(MATE_CIGAR_ATTRIBUTE), fragment.coordinates().Key,
                    fragment.averageBaseQual(), read.getMappingQuality(), suppData != null ? suppData.asCsv() : "N/A", read.getFlags()));

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

    private BufferedWriter initialiseResolvedReadWriter()
    {
        try
        {
            String filename = mConfig.formFilename("resolved_fragments");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("ReadId,Status,RemotePartition");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            BM_LOGGER.error(" failed to create resolved fragments writer: {}", e.toString());
        }

        return null;
    }

    public void writeResolvedReadData(final String readId, final FragmentStatus status, final String chrPartition)
    {
        if(mResolvedReadWriter == null)
            return;

        try
        {
            mResolvedReadWriter.write(format("%s,%s,%s", readId, status, chrPartition));
            mResolvedReadWriter.newLine();
        }
        catch(IOException e)
        {
            BM_LOGGER.error(" failed to write resolved fragment: {}", e.toString());
        }
    }

    public void closeInterimFiles()
    {
        if(mCandidateBamWriter != null)
        {
            BM_LOGGER.info("{} candidate reads written", mCandidateBamWriter.writeCount());
            mCandidateBamWriter.close();
        }

        if(mSupplementaryBamWriter != null)
        {
            BM_LOGGER.info("{} supplementary reads written", mSupplementaryBamWriter.writeCount());
            mSupplementaryBamWriter.close();
        }

        closeBufferedWriter(mResolvedReadWriter);
    }

    public void close()
    {
        BM_LOGGER.info("{} records written to BAM", mBamWriter.writeCount());
        mBamWriter.close();

        closeBufferedWriter(mReadWriter);
    }

    public Set<SAMRecord> readsWritten() { return mReadsWritten; }
}
