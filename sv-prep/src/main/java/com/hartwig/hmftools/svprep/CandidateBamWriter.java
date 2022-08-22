package com.hartwig.hmftools.svprep;

import static java.lang.String.format;

import static com.hartwig.hmftools.svprep.SpanningReadCache.chrFromChrPartition;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.WriteType.CACHE_BAM;
import static com.hartwig.hmftools.svprep.reads.ReadRecord.findReadIdTrimIndex;
import static com.hartwig.hmftools.svprep.reads.ReadRecord.trimReadId;
import static com.hartwig.hmftools.svprep.reads.ReadType.CANDIDATE_SUPPORT;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.svprep.reads.ReadGroup;
import com.hartwig.hmftools.svprep.reads.ReadGroupStatus;
import com.hartwig.hmftools.svprep.reads.ReadRecord;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class CandidateBamWriter
{
    private final SvConfig mConfig;
    private final Map<String,SAMFileWriter> mCandidatesWriters;
    private final Map<String,String> mCandidatesWriterBamFiles;

    private final Map<String,Set<String>> mChrJunctionReadIds;
    private int mReadIdTrimIndex;

    public CandidateBamWriter(final SvConfig config)
    {
        mConfig = config;
        mCandidatesWriters = Maps.newHashMap();
        mCandidatesWriterBamFiles = Maps.newHashMap();
        mChrJunctionReadIds = Maps.newHashMap();
        mReadIdTrimIndex = -1;
    }

    public void addJunctionReadId(final Set<String> remotePartitions, final String readId)
    {
        for(String remotePartition : remotePartitions)
        {
            String chromosome = chrFromChrPartition(remotePartition);

            Set<String> readIds = mChrJunctionReadIds.get(chromosome);

            if(readIds == null)
            {
                readIds = Sets.newHashSet();
                mChrJunctionReadIds.put(chromosome, readIds);
            }

            readIds.add(readId);
        }
    }

    public void writeCandidateRead(final ReadRecord read)
    {
        if(!mConfig.UseCacheBam)
            return;

        SAMFileWriter writer = mCandidatesWriters.get(read.Chromosome);

        if(writer == null)
        {
            SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));
            String bamFile = format("%s_%s.bam", mConfig.formFilename(CACHE_BAM), read.Chromosome);
            mCandidatesWriterBamFiles.put(read.Chromosome, bamFile);

            SAMFileHeader fileHeader = samReader.getFileHeader().clone();
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            writer = new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(bamFile));
            mCandidatesWriters.put(read.Chromosome, writer);
        }

        writer.addAlignment(read.record());
    }

    public void assignCandidateReads(final ResultsWriter resultsWriter)
    {
        if(mCandidatesWriters.isEmpty())
            return;

        SV_LOGGER.info("assigning candidate reads for {} chromosomes", mChrJunctionReadIds.size());

        // close before re-accessing
        mCandidatesWriters.values().forEach(x -> x.close());

        List<CandidateReadMatchTask> chromosomeTasks = Lists.newArrayList();

        for(Map.Entry<String,Set<String>> entry : mChrJunctionReadIds.entrySet())
        {
            String chromosome = entry.getKey();
            Set<String> junctionReadIds = entry.getValue();

            if(!mCandidatesWriterBamFiles.containsKey(chromosome))
            {
                junctionReadIds.clear();
                continue;
            }

            String bamFile = mCandidatesWriterBamFiles.get(chromosome);

            final SamReader samReader = SamReaderFactory.makeDefault()
                    .referenceSequence(new File(mConfig.RefGenomeFile))
                    .open(new File(bamFile));

            CandidateReadMatchTask chrTask = new CandidateReadMatchTask(chromosome, samReader, resultsWriter, junctionReadIds);
            chromosomeTasks.add(chrTask);
        }

        final List<Callable> callableList = chromosomeTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mConfig.Threads);

        SV_LOGGER.info("candidate reads assignment complete");
    }

    private class CandidateReadMatchTask implements Callable
    {
        private final String mChromosome;
        private final SamReader mSamReader;
        private final Set<String> mJunctionReadIds;
        private final ResultsWriter mResultsWriter;

        public CandidateReadMatchTask(
                final String chromosome, final SamReader samReader, final ResultsWriter resultsWriter, final Set<String> readIds)
        {
            mChromosome = chromosome;
            mSamReader = samReader;
            mResultsWriter = resultsWriter;
            mJunctionReadIds = readIds;
        }

        private static final int READ_GROUP_FLUSH = 1000;

        @Override
        public Long call()
        {
            SV_LOGGER.info("chr({}) assigning candidates from {} junction fragments", mChromosome, mJunctionReadIds.size());

            int matchedCandidates = 0;
            int recordCount = 0;

            SAMRecordIterator iter = mSamReader.iterator();

            List<ReadGroup> readGroups = Lists.newArrayList();

            while(iter.hasNext())
            {
                SAMRecord record = iter.next();
                ++recordCount;

                if(record.getDuplicateReadFlag() || record.isSecondaryAlignment())
                    continue;

                if(checkJunctionRead(record.getContig(), record.getReadName()))
                {
                    ++matchedCandidates;
                    ReadRecord read = ReadRecord.from(record);
                    read.setReadType(CANDIDATE_SUPPORT);
                    ReadGroup readGroup = new ReadGroup(read);
                    readGroup.markHasRemoteJunctionReads();
                    readGroup.setGroupState(ReadGroupStatus.EXPECTED);
                    readGroups.add(readGroup);

                    if(readGroups.size() >= READ_GROUP_FLUSH)
                    {
                        mResultsWriter.writeReadGroup(readGroups);
                        readGroups.clear();
                    }
                }

                if((recordCount % 100000) == 0)
                {
                    SV_LOGGER.debug("chr({}) processed {} candidate reads, matched({})",
                            mChromosome, recordCount, matchedCandidates);
                }
            }

            mResultsWriter.writeReadGroup(readGroups);

            SV_LOGGER.info("chr({}) matched and wrote {} candidate reads", mChromosome, matchedCandidates);

            return (long)0;
        }
    }

    private boolean checkJunctionRead(final String chromosome, final String readName)
    {
        String readId;
        if(mConfig.TrimReadId)
        {
            if(mReadIdTrimIndex < 0)
                mReadIdTrimIndex = findReadIdTrimIndex(readName);

            readId = trimReadId(readName, mReadIdTrimIndex);
        }
        else
        {
            readId = readName;
        }

        Set<String> readIds = mChrJunctionReadIds.get(chromosome);
        if(readIds == null)
            return false;

        return readIds != null && readIds.contains(readId);
    }
}
