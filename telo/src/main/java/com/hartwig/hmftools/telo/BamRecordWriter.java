package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.BlockingQueue;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

class TelBamRecord
{
    public SAMRecord samRecord;
    public boolean hasTeloContent = false;
    public boolean poison = false;
}

public class BamRecordWriter implements Runnable
{
    private final BlockingQueue<TelBamRecord> mTelBamRecordQ;

    private final Set<String> mIncompleteReadNames;
    private final Map<String, ReadGroup> mIncompleteReadGroups = new HashMap<>();
    private final Writer mCsvWriter;
    private final SAMFileWriter mBamFileWriter;
    private int mNumCompletedGroups = 0;
    private volatile boolean mProcessingMateRegions = false;

    public BamRecordWriter(final TeloConfig config, BlockingQueue<TelBamRecord> telBamRecordQ, Set<String> incompleteReadNames)
    {
        mTelBamRecordQ = telBamRecordQ;
        mIncompleteReadNames = incompleteReadNames;

        String sampleBamFileBasename = new File(config.BamFile).getName();
        final String csvOutputFile = config.OutputDir + "/" + sampleBamFileBasename + ".telo_read_data.csv.gz";
        mCsvWriter = createReadDataWriter(csvOutputFile);

        SamReader samReader = TeloUtils.openSamReader(config);

        final String bamOutputFile = config.OutputDir + "/" + sampleBamFileBasename + ".telo_bam.bam";
        mBamFileWriter = new SAMFileWriterFactory().makeBAMWriter(samReader.getFileHeader(), false, new File(bamOutputFile));
    }

    @Override
    public void run()
    {
        while(true)
        {
            TelBamRecord task;
            try
            {
                task = mTelBamRecordQ.take();
                if(task.poison)
                {
                    break;
                }
                processReadRecord(task.samRecord, task.hasTeloContent);
            }
            catch (InterruptedException e)
            {
                break;
            }
        }
    }

    public Map<String, ReadGroup> getIncompleteReadGroups() { return mIncompleteReadGroups; }
    public void setProcessingMateRegions(boolean b) { mProcessingMateRegions = b; }

    public void finish()
    {
        // we write out the final incomplete group
        mIncompleteReadGroups.values().forEach(this::writeReadGroup);

        try
        {
            mCsvWriter.close();
        }
        catch (IOException e)
        {
            throw new IllegalStateException("Could not close buffered writer: " + mCsvWriter + ": " + e.getMessage());
        }

        mBamFileWriter.close();
        TE_LOGGER.info("wrote {} read groups, complete({}), incomplete({})",
                mNumCompletedGroups + mIncompleteReadGroups.size(),
                mNumCompletedGroups, mIncompleteReadGroups.size());
    }

    private static Writer createReadDataWriter(final String outputFile)
    {
        try
        {
            Writer writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile, false)));
            writer.write("ReadId,Chromosome,PosStart,PosEnd,MateChr,MatePosStart,HasTeloContent");
            writer.write(",Cigar,InsertSize,FirstInPair,Unmapped,MateUnmapped,Flags,SuppData,CompleteFrag,ReadBases,BaseQualities");
            writer.write('\n');
            return writer;
        }
        catch (IOException e)
        {
            TE_LOGGER.error("failed to create read data writer: {}", e.toString());
            return null;
        }
    }

    private void processReadRecord(final SAMRecord record, boolean hasTelomereContent)
    {
        if(mProcessingMateRegions && !record.getReadPairedFlag())
        {
            return;
        }

        ReadGroup readGroup = mIncompleteReadGroups.get(record.getReadName());

        if(readGroup == null)
        {
            if (mProcessingMateRegions)
            {
                // do not create a new group if we are doing mate regions
                // the group would have already been created and written out
                return;
            }
            if (hasTelomereContent)
            {
                // cache if new
                readGroup = new ReadGroup(record);
                mIncompleteReadGroups.put(record.getReadName(), readGroup);
                mIncompleteReadNames.add(readGroup.id());
            }
            else
            {
                return;
            }
        }
        else
        {
            // a group already exist. This should be the mate of the record that is
            // already in the read group, check to see if that's the case
            // note: the following have to read how queryMate function is implemented in htsjdk SamReader to work out
            // what is the correct way to find the mate of a read. See htsjdk/samtools/SamReader.java
            SAMRecord mateRecord = readGroup.Reads.get(0);
            if(mateRecord.getFirstOfPairFlag() != record.getFirstOfPairFlag())
            {
                readGroup.Reads.add(record);
            }
        }

        if(readGroup.isComplete())
            processCompleteReadGroup(readGroup);
    }

    private void processCompleteReadGroup(final ReadGroup readGroup)
    {
        if (mIncompleteReadGroups.remove(readGroup.id()) != null)
        {
            mIncompleteReadNames.remove(readGroup.id());

            if (readGroup.Reads.size() > 2)
            {
                TE_LOGGER.debug("read group size: {}", readGroup.Reads.size());
            }

            writeReadGroup(readGroup);

            ++mNumCompletedGroups;
        }
    }

    public void writeReadGroup(final ReadGroup readGroup)
    {
        if(mCsvWriter != null)
        {
            try
            {
                boolean completeGroup = readGroup.isComplete();

                for(SAMRecord record : readGroup.Reads)
                {
                    ReadRecord readRecord = ReadRecord.from(record);
                    readRecord.setTeloContent(TeloUtils.hasTelomericContent(record.getReadString()));

                    mCsvWriter.write(String.format("%s,%s,%d,%d,%s,%d,%s",
                            readRecord.Id, readRecord.Chromosome, readRecord.PosStart, readRecord.PosEnd,
                            readRecord.mateChromosome(), readRecord.mateStartPosition(), readRecord.hasTeloContent()));

                    String suppAlignmentData = readRecord.hasSuppAlignment() ? readRecord.getSuppAlignment().replace(",", ";") : "";
                    suppAlignmentData = suppAlignmentData.replace(",", ";");

                    mCsvWriter.write(String.format(",%s,%d,%s,%s,%s,%d",
                            readRecord.Cigar.toString(), readRecord.fragmentInsertSize(), readRecord.isFirstOfPair(),
                            readRecord.isUnmapped(), readRecord.isMateUnmapped(), readRecord.flags()));

                    mCsvWriter.write(String.format(",%s,%s,%s,%s",
                            suppAlignmentData, completeGroup, readRecord.ReadBases, readRecord.BaseQualityString));

                    mCsvWriter.write('\n');
                }
            }
            catch(IOException e)
            {
                TE_LOGGER.error("failed to write read data file: {}", e.toString());
            }
        }

        if(mBamFileWriter != null)
        {
            for(SAMRecord samRecord : readGroup.Reads)
            {
                mBamFileWriter.addAlignment(samRecord);
            }
        }
    }
}
