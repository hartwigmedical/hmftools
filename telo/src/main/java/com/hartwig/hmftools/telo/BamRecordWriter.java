package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.BlockingQueue;
import java.util.zip.GZIPOutputStream;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.utils.FileWriterUtils;
import com.hartwig.hmftools.telo.analysers.TelomereReadsAnalyser;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
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
    private int mNumAcceptedReads = 0;
    private volatile boolean mProcessingMateRegions = false;
    private final TelomereReadsAnalyser mTelbamAnalyser = new TelomereReadsAnalyser();
    private final String mLengthCsvPath;

    public BamRecordWriter(final TeloConfig config, BlockingQueue<TelBamRecord> telBamRecordQ, Set<String> incompleteReadNames)
    {
        mTelBamRecordQ = telBamRecordQ;
        mIncompleteReadNames = incompleteReadNames;

        String sampleBamFileBasename = new File(config.BamFile).getName();
        // Remove the extension.
        int extensionIndex = sampleBamFileBasename.lastIndexOf('.');
        if (extensionIndex != -1)
            sampleBamFileBasename = sampleBamFileBasename.substring(0, extensionIndex);
        final String csvOutputFile = config.OutputDir + "/" + sampleBamFileBasename + ".telo_read_data.csv.gz";
        mCsvWriter = createReadDataWriter(csvOutputFile);

        SamReader samReader = TeloUtils.openSamReader(config);

        final String telbamPath = config.OutputDir + "/" + sampleBamFileBasename + ".telbam.bam";
        mBamFileWriter = new SAMFileWriterFactory().makeBAMWriter(samReader.getFileHeader(), false, new File(telbamPath));

        mLengthCsvPath = config.OutputDir + "/" + sampleBamFileBasename + ".tel_length.csv";
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
    public int getNumAcceptedReads() { return mNumAcceptedReads; }
    public void setProcessingMateRegions(boolean b) { mProcessingMateRegions = b; }

    public void finish()
    {
        // we write out the final incomplete group
        mIncompleteReadGroups.values().forEach(x-> { mTelbamAnalyser.onReadGroup(x); writeReadGroup(x); });

        for(ReadGroup rg : mIncompleteReadGroups.values())
        {
            TE_LOGGER.debug("incomplete read group: id={}", rg.getName());
            for(SAMRecord r : rg.Reads)
            {
                TE_LOGGER.debug("record: readid {}, suppl flag: {} suppl attri: {}, seq: {}, start: {}, end: {}, paired: {}",
                        r.getReadName(), r.getSupplementaryAlignmentFlag(), r.getStringAttribute(SAMTag.SA.name()), r.getReferenceName(),
                        r.getAlignmentStart(), r.getAlignmentEnd(), r.getReadPairedFlag());
            }
            for(SAMRecord r : rg.SupplementaryReads)
            {
                TE_LOGGER.debug("supplementary: readid {}, suppl flag: {} suppl attri: {}, seq: {}, start: {}, end: {}, paired: {}",
                        r.getReadName(), r.getSupplementaryAlignmentFlag(), r.getStringAttribute(SAMTag.SA.name()), r.getReferenceName(),
                        r.getAlignmentStart(), r.getAlignmentEnd(), r.getReadPairedFlag());
            }
        }

        try
        {
            mCsvWriter.close();
        }
        catch (IOException e)
        {
            throw new IllegalStateException("Could not close buffered writer: " + mCsvWriter + ": " + e.getMessage());
        }

        mBamFileWriter.close();
        writeTelLengthCsv();

        TE_LOGGER.info("wrote {} read groups, complete({}), incomplete({})",
                mNumCompletedGroups + mIncompleteReadGroups.size(),
                mNumCompletedGroups, mIncompleteReadGroups.size());

        TE_LOGGER.info("frag type counts: F1({}) F2({}) F4({})",
                mTelbamAnalyser.getFragmentTypeCount(ReadGroup.FragmentType.F1), mTelbamAnalyser.getFragmentTypeCount(ReadGroup.FragmentType.F2),
                mTelbamAnalyser.getFragmentTypeCount(ReadGroup.FragmentType.F4));
    }

    private static Writer createReadDataWriter(final String outputFile)
    {
        try
        {
            Writer writer = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFile, false)));
            writer.write("ReadId,Chromosome,PosStart,PosEnd,MateChr,MatePosStart,HasTeloContent");
            writer.write(",Cigar,InsertSize,FirstInPair,Unmapped,MateUnmapped,IsSupplementary,Flags");
            writer.write(",SuppData,CompleteFrag,ReadBases,BaseQualities");
            writer.write('\n');
            return writer;
        }
        catch (IOException e)
        {
            TE_LOGGER.error("failed to create read data writer: {}", e.toString());
            return null;
        }
    }

    private void processReadRecord(@NotNull final SAMRecord record, boolean hasTelomereContent)
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
                readGroup = new ReadGroup(record.getReadName());
                mIncompleteReadGroups.put(record.getReadName(), readGroup);
                mIncompleteReadNames.add(readGroup.getName());
            }
            else
            {
                return;
            }
        }
        if (!readGroup.contains(record))
        {
            acceptRead(readGroup, record);
        }
        assert(readGroup.invariant());

        if(readGroup.isComplete())
            processCompleteReadGroup(readGroup);
    }

    private void acceptRead(@NotNull final ReadGroup readGroup, @NotNull SAMRecord record)
    {
        assert(!readGroup.contains(record));

        if (record.getSupplementaryAlignmentFlag())
        {
            readGroup.SupplementaryReads.add(record);
        }
        else
        {
            readGroup.Reads.add(record);
        }

        ++mNumAcceptedReads;
    }

    private void processCompleteReadGroup(@NotNull final ReadGroup readGroup)
    {
        if (mIncompleteReadGroups.remove(readGroup.getName()) != null)
        {
            mIncompleteReadNames.remove(readGroup.getName());

            if (readGroup.Reads.size() > 2)
            {
                TE_LOGGER.debug("read group size: {}", readGroup.Reads.size());
            }

            mTelbamAnalyser.onReadGroup(readGroup);
            writeReadGroup(readGroup);

            ++mNumCompletedGroups;
        }
    }

    public void writeReadGroup(@NotNull final ReadGroup readGroup)
    {
        if(mCsvWriter != null)
        {
            try
            {
                boolean completeGroup = readGroup.isComplete();

                for(SAMRecord record : Iterables.concat(readGroup.Reads, readGroup.SupplementaryReads))
                {
                    ReadRecord readRecord = ReadRecord.from(record);
                    readRecord.setTeloContent(TeloUtils.hasTelomericContent(record.getReadString()));

                    mCsvWriter.write(String.format("\"%s\",%s,%d,%d,%s,%d,%s",
                            readRecord.Id, readRecord.Chromosome, readRecord.PosStart, readRecord.PosEnd,
                            readRecord.mateChromosome(), readRecord.mateStartPosition(), readRecord.hasTeloContent()));

                    String suppAlignmentData = readRecord.getSuppAlignment();
                    if (suppAlignmentData == null)
                    {
                        suppAlignmentData = "";
                    }

                    mCsvWriter.write(String.format(",%s,%d,%s,%s,%s,%s,%d",
                            readRecord.Cigar.toString(), readRecord.fragmentInsertSize(), readRecord.isFirstOfPair(),
                            readRecord.isUnmapped(), readRecord.isMateUnmapped(), record.getSupplementaryAlignmentFlag(),
                            readRecord.flags()));

                    mCsvWriter.write(String.format(",\"%s\",%s,\"%s\",\"%s\"",
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
            for(SAMRecord samRecord : Iterables.concat(readGroup.Reads, readGroup.SupplementaryReads))
            {
                mBamFileWriter.addAlignment(samRecord);
            }
        }
    }

    public void writeTelLengthCsv()
    {
        try
        {
            BufferedWriter lengthCsvWriter = FileWriterUtils.createBufferedWriter(mLengthCsvPath, false);
            lengthCsvWriter.write("F1,F2,F4\n");
            lengthCsvWriter.write(String.format("%d,%d,%d\n", mTelbamAnalyser.getFragmentTypeCount(ReadGroup.FragmentType.F1),
                    mTelbamAnalyser.getFragmentTypeCount(ReadGroup.FragmentType.F2),
                    mTelbamAnalyser.getFragmentTypeCount(ReadGroup.FragmentType.F4)));
            lengthCsvWriter.close();

        } catch (IOException e)
        {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
    }
}
