package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

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

    // one sam reader per thread
    private ConcurrentMap<Thread, SamReader> mThreadSamReaders;
    private SamReader mSamReader;
    private final BamRecordWriter mWriter;

    private final ConcurrentMap<String, ReadGroup> mIncompleteReadGroups;
    private int mCompletedGroups;

    public BamReader(final TeloConfig config, final BamRecordWriter writer, ConcurrentMap<Thread, SamReader> threadSamReaders,
            ConcurrentMap<String, ReadGroup> incompleteReadGroups)
    {
        mConfig = config;
        mBaseRegion = null;

        mThreadSamReaders = threadSamReaders;
        mWriter = writer;

        mIncompleteReadGroups = incompleteReadGroups;
        mCompletedGroups = 0;
    }

    public void setBaseRegion(final BaseRegion baseRegion)
    {
        mBaseRegion = baseRegion;
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
        // we use a thread based sam reader
        mSamReader = mThreadSamReaders.get(Thread.currentThread());

        if (mSamReader == null)
        {
            TE_LOGGER.info("creating sam reader");
            mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.SampleBamFile));
            mThreadSamReaders.put(Thread.currentThread(), mSamReader);
        }

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
        TE_LOGGER.info("processing unmapped reads, incomplete={}", mIncompleteReadGroups.size());

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
        // TODO - consider organising these into chromosomes and merging any close positions into a single region
        // Question: What happens if the mate is not mapped?
        List<BaseRegion> mateRegions = Lists.newArrayList();
        for(ReadGroup readGroup : mIncompleteReadGroups.values())
        {
            for(SAMRecord read : readGroup.Reads)
            {
                if(!read.getReadPairedFlag())
                {
                    continue;
                }

                if(!read.getMateUnmappedFlag())
                {
                    if(HumanChromosome.contains(read.getMateReferenceName()) && read.getMateAlignmentStart() > 0)
                    {
                        mateRegions.add(new BaseRegion(read.getMateReferenceName(), read.getMateAlignmentStart(), read.getMateAlignmentStart() + 1));
                        break;
                    }
                    else
                    {
                        TE_LOGGER.info("mate region: {}:{} not found in chromosome", read.getMateReferenceName(), read.getMateAlignmentStart());
                    }
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
        // we discard any supplementary / secondary alignments
        if(record.isSecondaryOrSupplementary())
            return;

        // duplicate reads are fine, cause they telomere all look similar we should just keep them.
        // most tools probably will not be able to correctly align them.

        ReadGroup readGroup = mIncompleteReadGroups.get(record.getReadName());
        boolean hasTeloContent = hasTelomericContent(record);

        if(!hasTeloContent && readGroup == null)
            return;

        if(readGroup == null)
        {
            // cache if new
            readGroup = new ReadGroup(record);
            mIncompleteReadGroups.put(record.getReadName(), readGroup);
        }
        else
        {
            // otherwise log details and discard if now complete
            //TE_LOGGER.info("add to existing read group {}", record.getReadName());
            readGroup.Reads.add(record);

        }
        if(readGroup.isComplete())
            processCompleteReadGroup(readGroup);
    }

    private void processIncompleteReadRecord(final SAMRecord record)
    {
        ReadGroup readGroup = mIncompleteReadGroups.get(record.getReadName());

        if(readGroup == null)
            return;

        readGroup.Reads.add(record);

        if(readGroup.isComplete())
            processCompleteReadGroup(readGroup);
    }

    private void processCompleteReadGroup(final ReadGroup readGroup)
    {
        if (mIncompleteReadGroups.remove(readGroup.id()) != null)
        {
            if(readGroup.id().equals("HSQ1008:208:C0VH6ACXX:7:1212:6586:41485"))
            {
                TE_LOGGER.info("here");
            }

            if (readGroup.Reads.size() > 2)
            {
                TE_LOGGER.info("read group size: {}", readGroup.Reads.size());
            }

            // readGroup.Reads.forEach(x -> x.markCompleteGroup());
            mWriter.writeReadGroup(readGroup);

            ++mCompletedGroups;
        }
    }
}
