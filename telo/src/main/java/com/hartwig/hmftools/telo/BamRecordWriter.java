package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.BlockingQueue;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.telo.analysers.TelomereReadsAnalyser;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;

import tech.tablesaw.api.BooleanColumn;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.StringColumn;

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
    private final tech.tablesaw.api.Table mReadDataTable = tech.tablesaw.api.Table.create("ReadData");
    private final SAMFileWriter mBamFileWriter;
    private int mNumCompletedGroups = 0;
    private int mNumAcceptedReads = 0;
    private volatile boolean mProcessingMateRegions = false;
    private final TelomereReadsAnalyser mTelbamAnalyser = new TelomereReadsAnalyser();
    private final String mReadDataCsvPath;
    private final String mLengthCsvPath;

    public BamRecordWriter(final TeloConfig config, BlockingQueue<TelBamRecord> telBamRecordQ, Set<String> incompleteReadNames)
    {
        mTelBamRecordQ = telBamRecordQ;
        mIncompleteReadNames = incompleteReadNames;

        SamReader samReader = TeloUtils.openSamReader(config);

        final String telbamPath = String.format("%s/%s.telo.%s.telbam.bam", config.OutputDir, config.SampleId, config.SampleType);
        mBamFileWriter = new SAMFileWriterFactory().makeBAMWriter(samReader.getFileHeader(), false, new File(telbamPath));

        mLengthCsvPath = String.format("%s/%s.telo.%s.tel_length.tsv", config.OutputDir, config.SampleId, config.SampleType);
        mReadDataCsvPath = String.format("%s/%s.telo.%s.read_data.tsv", config.OutputDir, config.SampleId, config.SampleType);

        setReadDataTableColumns();
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
    public void setProcessingMissingReadRegions(boolean b) { mProcessingMateRegions = b; }

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

        mBamFileWriter.close();
        writeReadDataToTsv();
        writeTelLengthCsv();

        TE_LOGGER.info("wrote {} read groups, complete({}), incomplete({})",
                mNumCompletedGroups + mIncompleteReadGroups.size(),
                mNumCompletedGroups, mIncompleteReadGroups.size());

        TE_LOGGER.info("frag type counts: F1({}) F2({}) F4({})",
                mTelbamAnalyser.getFragmentTypeCount(ReadGroup.FragmentType.F1), mTelbamAnalyser.getFragmentTypeCount(ReadGroup.FragmentType.F2),
                mTelbamAnalyser.getFragmentTypeCount(ReadGroup.FragmentType.F4));
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
        if(mReadDataTable != null)
        {
            for(SAMRecord record : Iterables.concat(readGroup.Reads, readGroup.SupplementaryReads))
            {
                addReadDataTableRow(record, readGroup.isComplete());
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
        final tech.tablesaw.api.Table teloLengthTable = tech.tablesaw.api.Table.create("TeloLength");
        teloLengthTable.addColumns(
                tech.tablesaw.api.IntColumn.create("F1"),
                tech.tablesaw.api.IntColumn.create("F2"),
                tech.tablesaw.api.IntColumn.create("F4"));

        tech.tablesaw.api.Row row = teloLengthTable.appendRow();
        row.setInt("F1", mTelbamAnalyser.getFragmentTypeCount(ReadGroup.FragmentType.F1));
        row.setInt("F2", mTelbamAnalyser.getFragmentTypeCount(ReadGroup.FragmentType.F2));
        row.setInt("F4", mTelbamAnalyser.getFragmentTypeCount(ReadGroup.FragmentType.F4));

        try
        {
            tech.tablesaw.io.csv.CsvWriteOptions writeOptions = tech.tablesaw.io.csv.CsvWriteOptions.builder(mLengthCsvPath).separator('\t').build();
            teloLengthTable.write().csv(writeOptions);

        } catch (IOException e)
        {
            throw new IllegalStateException("Could not save to csv file: " + mLengthCsvPath + ": " + e.getMessage());
        }
    }

    private void setReadDataTableColumns()
    {
        // add all the columns we need for the CSV
        mReadDataTable.addColumns(StringColumn.create("ReadId"),
                                StringColumn.create("Chromosome"),
                                IntColumn.create("PosStart"),
                                IntColumn.create("PosEnd"),
                                StringColumn.create("MateChr"),
                                IntColumn.create("MatePosStart"),
                                BooleanColumn.create("HasTeloContent"),
                                StringColumn.create("Cigar"),
                                IntColumn.create("InsertSize"),
                                BooleanColumn.create("FirstInPair"),
                                BooleanColumn.create("Unmapped"),
                                BooleanColumn.create("MateUnmapped"),
                                BooleanColumn.create("IsSupplementary"),
                                IntColumn.create("Flags"),
                                StringColumn.create("SuppData"),
                                BooleanColumn.create("CompleteFrag"),
                                StringColumn.create("ReadBases"),
                                StringColumn.create("BaseQualities"));
    }

    private void addReadDataTableRow(SAMRecord record, boolean readGroupIsComplete)
    {
        ReadRecord readRecord = ReadRecord.from(record);
        readRecord.setTeloContent(TeloUtils.hasTelomericContent(record.getReadString()));
        String suppAlignmentData = readRecord.getSuppAlignment();
        if (suppAlignmentData == null)
        {
            suppAlignmentData = "";
        }

        tech.tablesaw.api.Row row = mReadDataTable.appendRow();
        row.setString("ReadId", readRecord.Id);
        row.setString("Chromosome", readRecord.Chromosome);
        row.setInt("PosStart", readRecord.PosStart);
        row.setInt("PosEnd", readRecord.PosEnd);
        row.setString("MateChr", readRecord.mateChromosome());
        row.setInt("MatePosStart", readRecord.mateStartPosition());
        row.setBoolean("HasTeloContent", readRecord.hasTeloContent());
        row.setString("Cigar", readRecord.Cigar.toString());
        row.setInt("InsertSize", readRecord.fragmentInsertSize());
        row.setBoolean("FirstInPair", readRecord.isFirstOfPair());
        row.setBoolean("Unmapped", readRecord.isUnmapped());
        row.setBoolean("MateUnmapped", readRecord.isMateUnmapped());
        row.setBoolean("IsSupplementary", record.getSupplementaryAlignmentFlag());
        row.setInt("Flags", readRecord.flags());
        row.setString("SuppData", suppAlignmentData);
        row.setBoolean("CompleteFrag", readGroupIsComplete);
        row.setString("ReadBases", readRecord.ReadBases);
        row.setString("BaseQualities", readRecord.BaseQualityString);
    }

    private void writeReadDataToTsv()
    {
        try
        {
            tech.tablesaw.io.csv.CsvWriteOptions writeOptions = tech.tablesaw.io.csv.CsvWriteOptions.builder(mReadDataCsvPath).separator('\t').build();
            mReadDataTable.write().csv(writeOptions);
        }
        catch (IOException e)
        {
            throw new IllegalStateException("Could not save to tsv file: " + mReadDataCsvPath + ": " + e.getMessage());
        }
    }
}
