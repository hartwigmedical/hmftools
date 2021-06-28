package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;
import static com.hartwig.hmftools.telo.TeloConstants.CANONICAL_TELOMERE_SEQ;
import static com.hartwig.hmftools.telo.TeloConstants.CANONICAL_TELOMERE_SEQ_REV;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.samtools.BamSlicer;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamReader
{
    private final TeloConfig mConfig;

    private final SamReader mSamReader;
    // private final BamSlicer mBamSlicer;

    private final BufferedWriter mReadWriter;
    private int mReadCount;

    public BamReader(final TeloConfig config)
    {
        mConfig = config;

        mSamReader = mConfig.SampleBamFile != null ?
                SamReaderFactory.makeDefault()
                        .referenceSequence(new File(mConfig.RefGenomeFile))
                        .open(new File(mConfig.SampleBamFile)) : null;

        // mBamSlicer = new BamSlicer(1, true, true, false);

        mReadWriter = createReadDataWriter(mConfig);
        mReadCount = 0;
    }

    public void findTelomereContent()
    {
        // final QueryInterval[] queryIntervals = createIntervals(regions, samReader.getFileHeader());

        try (final SAMRecordIterator iterator = mSamReader.iterator()) // or mSamReader.queryUnmapped()
        {
            while (iterator.hasNext())
            {
                final SAMRecord record = iterator.next();

                if (passesFilters(record))
                {
                    processReadRecord(record);
                }
            }
        }

        TE_LOGGER.info("processed {} reads", mReadCount);
    }

    private boolean passesFilters(final SAMRecord record)
    {
        /*
        if(record.getMappingQuality() < 10)
            return false;

        if(record.getReadUnmappedFlag())
            return false;

        if(record.isSecondaryAlignment())
            return false;

        if(record.getSupplementaryAlignmentFlag())
            return false;

        if(record.getDuplicateReadFlag())
            return false;
        */

        return true;
    }

    private void processReadRecord(final SAMRecord record)
    {
        // analyse telomeric content..
        ++mReadCount;

        writeReadData(record);
    }

    private static BufferedWriter createReadDataWriter(final TeloConfig config)
    {
        if(!config.WriteReads)
            return null;

        try
        {
            final String outputFileName = config.OutputDir + "telo_read_data.csv";

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("ReadIndex,ReadId,Chromosome,PosStart,PosEnd,MateChr,MatePosStart,HasTeloContent");
            writer.write(",Cigar,InsertSize,FirstInPair,Unmapped,ReadReversed,Flags,SuppData,ReadBases");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            TE_LOGGER.error("failed to create read data writer: {}", e.toString());
            return null;
        }
    }

    private void writeReadData(final SAMRecord record)
    {
        if(mReadWriter == null)
            return;

        try
        {
            final String readBases = record.getReadString();
            boolean hasTeloContent = readBases.contains(CANONICAL_TELOMERE_SEQ) || readBases.contains(CANONICAL_TELOMERE_SEQ_REV);

            mReadWriter.write(String.format("%d,%s,%s,%d,%d,%s,%d,%s",
                    mReadCount, record.getReadName(), record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentEnd(),
                    record.getMateReferenceName(), record.getMateAlignmentStart(), hasTeloContent));

            String suppAlignmentData = record.hasAttribute("SA") ? record.getStringAttribute("SA") : "";

            mReadWriter.write(String.format(",%s,%d,%s,%s,%s,%d",
                    record.getCigarString(), record.getInferredInsertSize(), record.getFirstOfPairFlag(),
                    record.getReadUnmappedFlag(), record.getReadNegativeStrandFlag(), record.getFlags()));

            mReadWriter.write(String.format(",%s,%s",
                    suppAlignmentData, readBases));

            mReadWriter.newLine();
        }
        catch(IOException e)
        {
            TE_LOGGER.error("failed to write read data file: {}", e.toString());
        }
    }



}
