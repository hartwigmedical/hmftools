package com.hartwig.hmftools.svprep;

import static java.lang.String.format;

import static com.hartwig.hmftools.svprep.SpanningReadCache.chrFromChrPartition;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.WriteType.CACHE_BAM;
import static com.hartwig.hmftools.svprep.reads.ReadType.CANDIDATE_SUPPORT;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.svprep.reads.ReadGroup;
import com.hartwig.hmftools.svprep.reads.ReadGroupStatus;
import com.hartwig.hmftools.svprep.reads.ReadIdTrimmer;
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
    private ReadIdTrimmer mReadIdTrimmer;

    public CandidateBamWriter(final SvConfig config)
    {
        mConfig = config;
        mCandidatesWriters = Maps.newHashMap();
        mCandidatesWriterBamFiles = Maps.newHashMap();
        mChrJunctionReadIds = Maps.newHashMap();
        mReadIdTrimmer = new ReadIdTrimmer(mConfig.TrimReadId);
    }

    public void addJunctionReadId(final Set<String> remotePartitions, final String readId)
    {
        for(String remotePartition : remotePartitions)
        {
            String chromosome = chrFromChrPartition(remotePartition);

            Set<String> readIds = mChrJunctionReadIds.computeIfAbsent(chromosome, k -> Sets.newHashSet());

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

        SV_LOGGER.info("assigning candidate mate read to junctions");

        // close before re-accessing
        mCandidatesWriters.values().forEach(SAMFileWriter::close);

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

        TaskExecutor.executeTasks(chromosomeTasks, mConfig.Threads);

        SV_LOGGER.info("candidate reads assignment complete");

        // clean-up cache files
        if(!mConfig.NoCleanUp)
        {
            try
            {
                SV_LOGGER.debug("deleting {} candidate BAMs", mCandidatesWriterBamFiles.size());

                for(String cacheFile : mCandidatesWriterBamFiles.values())
                {
                    Files.delete(Paths.get(cacheFile));
                }
            }
            catch(IOException e)
            {
                SV_LOGGER.error("failed to delete candidate BAM: {}", e.toString());
            }
        }
    }

    private class CandidateReadMatchTask implements Callable<Long>
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
            SV_LOGGER.debug("chr({}) assigning candidates from {} junction fragments", mChromosome, mJunctionReadIds.size());

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

            SV_LOGGER.debug("chr({}) matched and wrote {} candidate reads", mChromosome, matchedCandidates);

            return (long)0;
        }
    }

    private boolean checkJunctionRead(final String chromosome, final String readName)
    {
        String readId = mReadIdTrimmer.trim(readName);

        Set<String> readIds = mChrJunctionReadIds.get(chromosome);
        if(readIds == null)
            return false;

        return readIds != null && readIds.contains(readId);
    }
}
