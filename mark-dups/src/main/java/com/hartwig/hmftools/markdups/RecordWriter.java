package com.hartwig.hmftools.markdups;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UMI_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UMI_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.PRIMARY;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.UNSET;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.readToString;
import static com.hartwig.hmftools.markdups.ReadOutput.DUPLICATES;
import static com.hartwig.hmftools.markdups.ReadOutput.MISMATCHES;
import static com.hartwig.hmftools.markdups.ReadOutput.NONE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.samtools.UmiReadType;
import com.hartwig.hmftools.markdups.common.Fragment;
import com.hartwig.hmftools.markdups.common.FragmentStatus;
import com.hartwig.hmftools.markdups.consensus.DuplicateGroup;

import org.jetbrains.annotations.Nullable;

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
    private final BufferedWriter mReadWriter;

    private boolean mCacheReads;
    private final Set<String> mExpectedReads;
    private int mNonConsensusReadCount;
    private int mConsensusReadCount;

    private int mExpectedReadsSize;

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
        mCacheReads = config.runReadChecks();
        mNonConsensusReadCount = 0;
        mConsensusReadCount = 0;
        mExpectedReadsSize = 0;

        if(mConfig.WriteBam)
        {
            String bamFilename = formBamFilename("mark_dups");
            MD_LOGGER.info("writing new BAM file: {}", bamFilename);
            mBamWriter = new BamWriter(bamFilename);
        }
        else
        {
            mBamWriter = new BamWriter(null);
        }

        mReadWriter = initialiseReadWriter();
        mExpectedReads = Sets.newHashSet();
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

    public int recordWriteCount() { return mNonConsensusReadCount; }
    public int recordWriteCountConsensus() { return mConsensusReadCount; }

    public synchronized void writeFragments(final List<Fragment> fragments, boolean excludeUmis)
    {
        fragments.forEach(x -> doWriteFragment(x, excludeUmis));
    }

    public synchronized void writeFragment(final Fragment fragment) { doWriteFragment(fragment, true); }

    public synchronized void writeRead(final SAMRecord read, final FragmentStatus fragmentStatus)
    {
        writeRead(read, fragmentStatus, null);
    }

    public synchronized void writeDuplicateGroup(final DuplicateGroup group, final List<SAMRecord> completeReads)
    {
        for(SAMRecord read : completeReads)
        {
            if(read.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
            {
                mBamWriter.writeRecord(read);
                ++mConsensusReadCount;

                if(mReadWriter != null)
                {
                    writeReadData(read, PRIMARY, group.coordinatesKey(), 0, group.id());
                }

                continue;
            }

            if(mConfig.UMIs.Enabled)
                read.setAttribute(UMI_ATTRIBUTE, group.id());

            writeRead(read, DUPLICATE, group.coordinatesKey(), 0, group.id());
        }
    }

    private void doWriteFragment(final Fragment fragment, boolean excludeUmis)
    {
        if(excludeUmis && fragment.umi() != null) // reads in UMI groups are only written as a complete group
            return;

        if(fragment.readsWritten())
        {
            MD_LOGGER.error("fragment({}) reads already written", fragment);
            return;
        }

        fragment.setReadWritten();
        fragment.reads().forEach(x -> writeRead(x, fragment.status(), fragment));
    }

    private void writeRead(final SAMRecord read, final FragmentStatus fragmentStatus, @Nullable final Fragment fragment)
    {
        writeRead(
                read, fragmentStatus,
                fragment != null ? fragment.coordinates().Key : "",
                fragment != null ? fragment.averageBaseQual() : 0,
                fragment != null ? fragment.umi() : "");
    }

    private void writeRead(
            final SAMRecord read, final FragmentStatus fragmentStatus, final String fragmentCoordinates,
            final double avgBaseQual, final String umiId)
    {
        ++mNonConsensusReadCount;

        writeReadData(read, fragmentStatus, fragmentCoordinates, avgBaseQual, umiId);

        read.setDuplicateReadFlag(fragmentStatus == DUPLICATE); // overwrite any existing status

        mBamWriter.writeRecord(read);

        if(mCacheReads)
            removeWrittenRead(read);
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
            writer.write(",InsertSize,MateChr,MatePosStart,Duplicate,CalcDuplicate,MateCigar,Coords");

            if(mConfig.UMIs.Enabled)
                writer.write(",Umi,UmiType");

            writer.write(",AvgBaseQual,MapQual,SuppData,Flags,FirstInPair,ReadReversed,Unmapped,MateUnmapped,Supplementary,Secondary");

            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            MD_LOGGER.error(" failed to create read writer: {}", e.toString());
        }

        return null;
    }

    private void writeReadData(
            final SAMRecord read, final FragmentStatus fragmentStatus, final String fragmentCoordinates,
            final double avgBaseQual, final String umiId)
    {
        if(mReadWriter == null)
            return;

        if(mConfig.LogReadType == DUPLICATES)
        {
            if(!read.getDuplicateReadFlag() && !fragmentStatus.isDuplicate())
                return;
        }
        else if(mConfig.LogReadType == MISMATCHES)
        {
            if(fragmentStatus != UNSET)
            {
                if(read.getDuplicateReadFlag() == (fragmentStatus == DUPLICATE))
                    return;
            }
        }

        try
        {
            mReadWriter.write(format("%s,%s,%d,%d,%s",
                    read.getReadName(), read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(), read.getCigar()));

            SupplementaryReadData suppData = SupplementaryReadData.from(read.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));

            mReadWriter.write(format(",%d,%s,%d,%s,%s,%s,%s",
                    abs(read.getInferredInsertSize()), read.getMateReferenceName(), read.getMateAlignmentStart(),
                    read.getDuplicateReadFlag(), fragmentStatus, read.hasAttribute(MATE_CIGAR_ATTRIBUTE), fragmentCoordinates));

            if(mConfig.UMIs.Enabled)
            {
                String umiType = read.getStringAttribute(UMI_TYPE_ATTRIBUTE);
                mReadWriter.write(format(",%s,%s", umiId != null ? umiId : "", umiType != null ? umiType : UmiReadType.NONE));
            }

            mReadWriter.write(format(",%.2f,%d,%s,%d",
                    avgBaseQual, read.getMappingQuality(), suppData != null ? suppData.asCsv() : "N/A", read.getFlags()));

            mReadWriter.write(format(",%s,%s,%s,%s,%s,%s",
                    read.getFirstOfPairFlag(), read.getReadNegativeStrandFlag(), read.getReadUnmappedFlag(),
                    read.getMateUnmappedFlag(), read.getSupplementaryAlignmentFlag(), read.isSecondaryAlignment()));

            mReadWriter.newLine();
        }
        catch(IOException e)
        {
            MD_LOGGER.error(" failed to write read data: {}", e.toString());
        }
    }

    public void close()
    {
        MD_LOGGER.info("{} records written to BAM", mBamWriter.writeCount());
        mBamWriter.close();

        closeBufferedWriter(mReadWriter);
    }

    public synchronized void registerRead(final SAMRecord read)
    {
        if(!mCacheReads)
            return;

        mExpectedReads.add(formCachedReadString(read));
    }

    private static final int LOG_CACHE_COUNT = 100000;

    private void removeWrittenRead(final SAMRecord read)
    {
        if(!mCacheReads)
            return;

        String readStr = formCachedReadString(read);

        if(!mExpectedReads.contains(readStr))
        {
            MD_LOGGER.warn("writing uncached read({}): {}", readStr, readToString(read));
        }
        else
        {
            mExpectedReads.remove(readStr);
        }

        if(abs(mExpectedReads.size() - mExpectedReadsSize) >= LOG_CACHE_COUNT)
        {
            mExpectedReadsSize = mExpectedReads.size();
            MD_LOGGER.info("record-writer cached expected reads count({})", mExpectedReadsSize);
        }
    }

    public void logUnwrittenReads()
    {
        if(!mCacheReads || mExpectedReads.isEmpty())
            return;

        MD_LOGGER.warn("unwritten read count({})", mExpectedReads.size());

        for(String readStr : mExpectedReads)
        {
            MD_LOGGER.warn("unwritten read({})", readStr);
        }
    }

    private static String formCachedReadString(final SAMRecord read)
    {
        boolean unmapped = read.getReadUnmappedFlag();

        return format("%s_%s_%d_%s_%s_%s",
                read.getReadName(),
                unmapped ? "unmapped" : read.getReferenceName(),
                unmapped ? 0 : read.getAlignmentStart(),
                read.getReadNegativeStrandFlag() ? "fwd" : "rev",
                read.getFirstOfPairFlag() ? "R1" : "R2",
                read.getSupplementaryAlignmentFlag() ? "supp" : "prim");
    }

    @VisibleForTesting
    public void setCacheReads() { mCacheReads = true;}

    @VisibleForTesting
    public void resetRecordWriteCounts()
    {
        mConsensusReadCount = 0;
        mNonConsensusReadCount = 0;
    }
}
