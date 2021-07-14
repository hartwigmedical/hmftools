package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamReader implements Callable
{
    private final TeloConfig mConfig;

    private BaseRegion mBaseRegion;

    private final SamReader mSamReader;
    private final BamRecordWriter mWriter;

    private final Map<String,ReadGroup> mIncompleteReadGroups;
    private int mCompletedGroups;

    // metrics
    private int mReadCount;

    public BamReader(final TeloConfig config, final BamRecordWriter writer)
    {
        mConfig = config;
        mBaseRegion = null;

        mSamReader = mConfig.SampleBamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.SampleBamFile)) : null;

        mWriter = writer;

        mIncompleteReadGroups = Maps.newHashMap();
        mReadCount = 0;
        mCompletedGroups = 0;
    }

    public void setBaseRegion(final BaseRegion baseRegion)
    {
        mBaseRegion = baseRegion;
    }

    public void setIncompleteReadGroups(final Map<String,ReadGroup> incompleteReadGroups)
    {
        mIncompleteReadGroups.clear();
        incompleteReadGroups.putAll(incompleteReadGroups);
    }

    public Map<String,ReadGroup> getIncompleteReadGroups() { return mIncompleteReadGroups; }

    @Override
    public Long call()
    {
        findTelomereContent();
        return (long)0;
    }

    public void findTelomereContent()
    {
        mCompletedGroups = 0;

        if(mBaseRegion != null)
        {
            processBamByRegion();
        }
        else
        {
            processUnmappedReads();
            processIncompleteReadGroups();
        }
    }

    private void processBamByRegion()
    {
        TE_LOGGER.info("processing region({})", mBaseRegion.toString());

        BamSlicer bamSlicer = new BamSlicer(1, false, false, false);
        bamSlicer.slice(mSamReader, Lists.newArrayList(mBaseRegion), this::processReadRecord);

        TE_LOGGER.info("processed region({}) reads(complete={} incomplete={})",
                mBaseRegion.toString(), mCompletedGroups, mIncompleteReadGroups.size());
    }

    private void processUnmappedReads()
    {
        TE_LOGGER.info("processing unmapped reads");

        int readCount = 0;

        try (final SAMRecordIterator iterator = mSamReader.queryUnmapped())
        {
            while (iterator.hasNext())
            {
                final SAMRecord record = iterator.next();
                ++readCount;

                if (passesFilters(record))
                {
                    processReadRecord(record);
                }
            }
        }

        TE_LOGGER.info("processed {} unmapped reads(complete={} incomplete={})",
                readCount, mCompletedGroups, mIncompleteReadGroups.size());
    }

    private void processIncompleteReadGroups()
    {
        TE_LOGGER.info("processing {} incomplete reads", mIncompleteReadGroups.size());

        // TODO - consider organising these into chromosomes and merging any close positions into a single region
        List<BaseRegion> mateRegions = Lists.newArrayList();
        for(ReadGroup readGroup : mIncompleteReadGroups.values())
        {
            // TODO - are these supplementary? if so should merge their
            for(ReadRecord read : readGroup.Reads)
            {
                if(HumanChromosome.contains(read.mateChromosome()) && read.mateStartPosition() > 0)
                {
                    mateRegions.add(new BaseRegion(read.mateChromosome(), read.mateStartPosition(), read.mateStartPosition() + 1));
                    break;
                }
            }
        }

        if(!mateRegions.isEmpty())
        {
            TE_LOGGER.info("processing {} mate regions", mIncompleteReadGroups.size());

            BamSlicer bamSlicer = new BamSlicer(1, false, false, false);
            mateRegions.forEach(x -> bamSlicer.slice(mSamReader, Lists.newArrayList(x), this::processIncompleteReadRecord));
        }
    }

    private boolean hasTelomericContent(final SAMRecord record)
    {
        return TeloUtils.hasTelomericContent(record.getReadString());
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
        ++mReadCount;

        if(!hasTelomericContent(record))
            return;

        // look up any existing reads with the same ID
        ReadRecord readRecord = ReadRecord.from(record);
        readRecord.setTeloContent(true);

        ReadGroup readGroup = mIncompleteReadGroups.get(readRecord.Id);

        if(readGroup == null)
        {
            // cache if new
            readGroup = new ReadGroup(readRecord);
            mIncompleteReadGroups.put(readRecord.Id, readGroup);
        }
        else
        {
            // otherwise log details and discard if now complete
            readGroup.Reads.add(ReadRecord.from(record));

            if(readGroup.isComplete())
                processCompleteReadGroup(readGroup);
        }
    }

    private void processIncompleteReadRecord(final SAMRecord record)
    {
        ReadGroup readGroup = mIncompleteReadGroups.get(record.getReadName());

        if(readGroup == null)
            return;

        ReadRecord readRecord = ReadRecord.from(record);
        readRecord.setTeloContent(hasTelomericContent(record));

        readGroup.Reads.add(readRecord);

        if(readGroup.isComplete())
            processCompleteReadGroup(readGroup);
    }

    private void processCompleteReadGroup(final ReadGroup readGroup)
    {
        mIncompleteReadGroups.remove(readGroup.id());

        readGroup.Reads.forEach(x -> x.markCompleteGroup());
        mWriter.writeReadGroup(readGroup);

        ++mCompletedGroups;
    }
}
