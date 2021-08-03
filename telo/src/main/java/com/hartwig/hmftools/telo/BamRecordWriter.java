package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.BlockingQueue;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

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
    private final BufferedWriter mCsvWriter;
    private final SAMFileWriter mBamFileWriter;
    private int mCompletedGroups = 0;
    private volatile boolean mProcessingMateRegions = false;

    public BamRecordWriter(final TeloConfig config, BlockingQueue<TelBamRecord> telBamRecordQ, Set<String> incompleteReadNames)
    {
        mTelBamRecordQ = telBamRecordQ;
        mIncompleteReadNames = incompleteReadNames;

        String sampleBamFileBasename = new File(config.SampleBamFile).getName();
        final String csvOutputFile = config.OutputDir + "/" + sampleBamFileBasename + ".telo_read_data.csv";
        mCsvWriter = createReadDataWriter(csvOutputFile);

        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(config.RefGenomeFile)).open(new File(config.SampleBamFile));

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

    public Map<String, ReadGroup> getmIncompleteReadGroups() { return mIncompleteReadGroups; }
    public void setProcessingMateRegions(boolean b) { mProcessingMateRegions = b; }

    public void close()
    {
        closeBufferedWriter(mCsvWriter);
        mBamFileWriter.close();
    }

    private static BufferedWriter createReadDataWriter(final String outputFile)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);
            writer.write("ReadId,Chromosome,PosStart,PosEnd,MateChr,MatePosStart,HasTeloContent");
            writer.write(",Cigar,InsertSize,FirstInPair,Unmapped,MateUnmapped,Flags,SuppData,CompleteFrag,ReadBases");
            writer.newLine();
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
            readGroup.Reads.add(record);
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
                TE_LOGGER.info("read group size: {}", readGroup.Reads.size());
            }

            writeReadGroup(readGroup);

            ++mCompletedGroups;
        }
    }

    public void writeReadGroup(final ReadGroup readGroup)
    {
        if(mCsvWriter != null)
        {
            try
            {
                boolean completeGroup = readGroup.isComplete();

                for(SAMRecord samRecord : readGroup.Reads)
                {
                    ReadRecord readRecord = ReadRecord.from(samRecord);
                    readRecord.setTeloContent(TeloUtils.hasTelomericContent(samRecord.getReadString()));

                    mCsvWriter.write(String.format("%s,%s,%d,%d,%s,%d,%s",
                            readRecord.Id, readRecord.Chromosome, readRecord.PosStart, readRecord.PosEnd,
                            readRecord.mateChromosome(), readRecord.mateStartPosition(), readRecord.hasTeloContent()));

                    String suppAlignmentData = readRecord.hasSuppAlignment() ? readRecord.getSuppAlignment().replace(",", ";") : "";
                    suppAlignmentData = suppAlignmentData.replace(",", ";");

                    mCsvWriter.write(String.format(",%s,%d,%s,%s,%s,%d",
                            readRecord.Cigar.toString(), readRecord.fragmentInsertSize(), readRecord.isFirstOfPair(),
                            readRecord.isUnmapped(), readRecord.isMateUnmapped(), readRecord.flags()));

                    mCsvWriter.write(String.format(",%s,%s,%s",
                            suppAlignmentData, completeGroup, readRecord.ReadBases));

                    mCsvWriter.newLine();
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

    public void writeAllIncompleteReadGroups()
    {
        // we write out the final incomplete group
        TE_LOGGER.info("writing final {} incomplete groups", mIncompleteReadGroups.size());
        mIncompleteReadGroups.values().forEach(this::writeReadGroup);
    }
}
