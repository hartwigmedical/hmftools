package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.esvee.prep.SpanningReadCache.chrFromChrPartition;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.CACHE_BAM;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.esvee.prep.types.ReadGroup;
import com.hartwig.hmftools.esvee.prep.types.ReadGroupStatus;
import com.hartwig.hmftools.esvee.common.ReadIdTrimmer;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;
import com.hartwig.hmftools.esvee.prep.types.ReadType;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class CandidateBamWriter
{
    private final PrepConfig mConfig;
    private final Map<String,SAMFileWriter> mCandidatesWriters;
    private final Map<String,String> mCandidatesWriterBamFiles;

    private final Map<String,Set<String>> mChrJunctionReadIds;
    private ReadIdTrimmer mReadIdTrimmer;

    public CandidateBamWriter(final PrepConfig config)
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

            Set<String> readIds = mChrJunctionReadIds.get(chromosome);

            if(readIds == null)
            {
                readIds = Sets.newHashSet();
                mChrJunctionReadIds.put(chromosome, readIds);
            }

            readIds.add(readId);
        }
    }

    public void writeCandidateRead(final PrepRead read)
    {
        if(!mConfig.UseCacheBam)
            return;

        SAMFileWriter writer = mCandidatesWriters.get(read.Chromosome);

        if(writer == null)
        {
            SamReader samReader = SamReaderFactory.makeDefault()
                    .referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.bamFile()));

            String bamFile = String.format("%s_%s.bam", mConfig.formFilename(CACHE_BAM), read.Chromosome);
            mCandidatesWriterBamFiles.put(read.Chromosome, bamFile);

            SAMFileHeader fileHeader = samReader.getFileHeader().clone();
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            // add read group info if reads are from multiple input BAMs
            if(mConfig.BamFiles.size() > 1)
            {
                for(int i = 1; i < mConfig.BamFiles.size(); ++i)
                {
                    SamReader nextReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile))
                            .open(new File(mConfig.BamFiles.get(i)));

                    for(SAMReadGroupRecord readGroupRecord : nextReader.getFileHeader().getReadGroups())
                    {
                        if(!fileHeader.getReadGroups().contains(readGroupRecord))
                            fileHeader.addReadGroup(readGroupRecord);
                    }

                    final SAMProgramRecord nextProgramRecord = nextReader.getFileHeader().getProgramRecords().get(0);
                    String newProgramId = String.format("%s.%d", nextProgramRecord.getId(), i);

                    fileHeader.addProgramRecord(new SAMProgramRecord(newProgramId, nextProgramRecord));
                }
            }

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

        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
            System.exit(1);

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
                    PrepRead read = PrepRead.from(record);
                    read.setReadType(ReadType.CANDIDATE_SUPPORT);
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
