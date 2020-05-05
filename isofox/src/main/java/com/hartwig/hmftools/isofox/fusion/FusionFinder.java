package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.abs;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.fusion.FusionFragment.isRealignedFragmentCandidate;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formLocation;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BamSlicer;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.TransExonRef;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class FusionFinder
{
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final Map<String,List<ReadRecord>> mReadsMap;
    private final Map<String,Map<Integer,List<EnsemblGeneData>>> mChrGeneCollectionMap;
    private final Map<String,Map<Integer,BaseDepth>> mChrGeneDepthMap;

    private List<FusionTask> mFusionTasks;
    private final FusionWriter mFusionWriter;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final PerformanceCounter mPerfCounter;

    public FusionFinder(final IsofoxConfig config, final EnsemblDataCache geneTransCache)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mReadsMap = Maps.newHashMap();
        mChrGeneCollectionMap = Maps.newHashMap();
        mChrGeneDepthMap = Maps.newHashMap();
        mFusionTasks = Lists.newArrayList();

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, true, true);

        mPerfCounter = new PerformanceCounter("FusionTotal");
        mFusionWriter = new FusionWriter(mConfig);
    }

    public void addChimericReads(final Map<String,List<ReadRecord>> chimericReadMap)
    {
        mergeChimericReadMaps(mReadsMap, chimericReadMap);
    }

    public static void addChimericReads(final Map<String,List<ReadRecord>> chimericReadMap, final ReadRecord read)
    {
        List<ReadRecord> chimericReads = chimericReadMap.get(read.Id);
        if (chimericReads == null)
            chimericReadMap.put(read.Id, Lists.newArrayList(read));
        else
            chimericReads.add(read);
    }

    public void loadChimericReads()
    {
        mReadsMap.putAll(FusionWriter.loadChimericReads(mConfig.Fusions.ReadsFile));
    }

    public void addChromosomeGeneCollections(final String chromosome, final Map<Integer,List<EnsemblGeneData>> geneCollectionMap)
    {
        mChrGeneCollectionMap.put(chromosome, geneCollectionMap);
    }

    public void addChromosomeGeneDepth(final String chromosome, final Map<Integer,BaseDepth> geneDepthMap)
    {
        mChrGeneDepthMap.put(chromosome, geneDepthMap);
    }

    public void findFusions()
    {
        ISF_LOGGER.info("processing {} chimeric read groups", mReadsMap.size());

        mPerfCounter.start();

        int unpairedReads = 0;
        int filteredFragments = 0;
        int duplicates = 0;
        int skipped = 0;
        int fragments = 0;

        Map<String,List<FusionFragment>> chrPairFragments = Maps.newHashMap();

        for(Map.Entry<String,List<ReadRecord>> entry : mReadsMap.entrySet())
        {
            final List<ReadRecord> reads = entry.getValue();

            if(reads.stream().anyMatch(x -> x.isDuplicate()))
            {
                ++duplicates;
                continue;
            }

            recoverSupplementaryReads(reads);

            String readGroupStatus = "";

            if(reads.size() == 1)
            {
                if(skipUnpairedRead(reads.get(0).mateChromosome(), reads.get(0).mateStartPosition()))
                {
                    ++skipped;
                    continue;
                }
            }

            if(reads.size() == 1)
            {
                ++unpairedReads;
                readGroupStatus = "UNPAIRED";
            }
            else if(!isValidFragment(reads))
            {
                ++filteredFragments;
                readGroupStatus = "INVALID";
            }
            else
            {
                FusionFragment fragment = new FusionFragment(reads);
                ++fragments;

                final String chrPair = formChromosomePair(fragment.chromosomes()[SE_START], fragment.chromosomes()[SE_END]);
                List<FusionFragment> fragmentList = chrPairFragments.get(chrPair);

                if(fragmentList == null)
                {
                    chrPairFragments.put(chrPair, Lists.newArrayList(fragment));
                }
                else
                {
                    fragmentList.add(fragment);
                }

                // will write read data after further evaluation of fusions
                continue;
            }

            mFusionWriter.writeReadData(reads, readGroupStatus);
        }

        int chrPairCount = 0;
        int taskId = 0;
        int pairsPerThread = mConfig.Threads > 1 ? round(chrPairFragments.size() / mConfig.Threads) : chrPairFragments.size();

        List<FusionFragment> allFragments = Lists.newArrayList();
        mFusionTasks.clear();

        for(Map.Entry<String,List<FusionFragment>> entry : chrPairFragments.entrySet())
        {
            allFragments.addAll(entry.getValue());

            ++chrPairCount;
            if(chrPairCount >= pairsPerThread || chrPairCount == chrPairFragments.size())
            {
                mFusionTasks.add(new FusionTask(taskId++, mConfig, mGeneTransCache, mChrGeneDepthMap, allFragments, mFusionWriter));
                allFragments = Lists.newArrayList();
                chrPairCount = 0;
            }
        }

        ISF_LOGGER.info("chimeric fragments({} unpaired={} dups={} skip={} filtered={} candidates={}) chrPairs({}) tasks({})",
                mReadsMap.size(), unpairedReads, skipped, duplicates, filteredFragments, fragments,
                chrPairFragments.size(), mFusionTasks.size());

        if(mFusionTasks.isEmpty())
        {
            ISF_LOGGER.warn("no fusion tasks created");
            return;
        }
        else
        {
            executeFusionTasks();
            logPerformanceStats();
        }

        mFusionWriter.close();

        mPerfCounter.stop();

        if(mConfig.Fusions.PerformanceStats)
            mPerfCounter.logStats();

        ISF_LOGGER.info("fusion calling complete");
    }

    private boolean executeFusionTasks()
    {
        if(mConfig.Threads <= 1)
        {
            mFusionTasks.forEach(x -> x.call());
            return true;
        }

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("IsofoxFusions-%d").build();

        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);
        List<FutureTask> threadTaskList = new ArrayList<FutureTask>();

        for(FusionTask fusionTask : mFusionTasks)
        {
            FutureTask futureTask = new FutureTask(fusionTask);

            threadTaskList.add(futureTask);
            executorService.execute(futureTask);
        }

        if(!checkThreadCompletion(threadTaskList))
        {
            return false;
        }

        executorService.shutdown();
        return true;
    }

    private boolean checkThreadCompletion(final List<FutureTask> taskList)
    {
        try
        {
            for (FutureTask futureTask : taskList)
            {
                futureTask.get();
            }
        }
        catch (Exception e)
        {
            ISF_LOGGER.error("task execution error: {}", e.toString());
            e.printStackTrace();
            return false;
        }

        return true;
    }

    private void recoverSupplementaryReads(final List<ReadRecord> reads)
    {
        if(reads.size() >= 3)
            return;

        for(final ReadRecord read : reads)
        {
            if(read.getSuppAlignment() == null)
                continue;

            SupplementaryReadData suppData = SupplementaryReadData.from(read.getSuppAlignment());

            if(suppData == null || skipUnpairedRead(suppData.Chromosome, suppData.Position))
                return;

            final List<ReadRecord> missingReads = findMissingReads(read.Id, suppData.Chromosome, suppData.Position, reads.size() == 1);

            if(!missingReads.isEmpty())
            {
                ISF_LOGGER.debug("id({}) recovered {} missing read: {}", read.Id, missingReads.size());
                reads.addAll(missingReads);
            }

            return;
        }
    }

    private List<ReadRecord> findMissingReads(final String readId, final String chromosome, int position, boolean expectPair)
    {
        final QueryInterval[] queryIntervals = new QueryInterval[1];
        int chrIndex = mSamReader.getFileHeader().getSequenceIndex(chromosome);

        // build a buffer around this position to pick up a second read
        if(expectPair)
        {
            int buffer = mConfig.MaxFragmentLength;
            queryIntervals[0] = new QueryInterval(chrIndex, position - buffer, position + buffer);
        }
        else
        {
            queryIntervals[0] = new QueryInterval(chrIndex, position, position);
        }

        final List<SAMRecord> records = mBamSlicer.slice(mSamReader, queryIntervals);

        return records.stream().filter(x -> x.getReadName().equals(readId)).map(x -> ReadRecord.from(x)).collect(Collectors.toList());
    }

    private boolean isValidFragment(final List<ReadRecord> reads)
    {
        if(!mConfig.ExcludedGeneIds.isEmpty())
        {
            for(ReadRecord read : reads)
            {
                for(List<TransExonRef> transExonList : read.getTransExonRefs().values())
                {
                    if(transExonList.stream().anyMatch(x -> mConfig.ExcludedGeneIds.contains(x.GeneId)))
                        return false;
                }
            }
        }

        Set<String> chrGeneSet = Sets.newHashSetWithExpectedSize(3);

        for(ReadRecord read : reads)
        {
            chrGeneSet.add(formLocation(read.Chromosome, read.getGeneCollecton()));
            if(chrGeneSet.size() == 3)
                return false;
        }

        if(chrGeneSet.size() == 2)
            return true;

        if(reads.stream().anyMatch(x -> isRealignedFragmentCandidate(x)))
            return true;

        return false;
    }

    private boolean skipUnpairedRead(final String otherChromosome, int otherPosition)
    {
        if(!HumanChromosome.contains(otherChromosome))
            return true;

        if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(otherChromosome))
            return true;

        if(!mConfig.RestrictedGeneIds.isEmpty())
        {
            final List<EnsemblGeneData> chrGenes = mGeneTransCache.getChrGeneDataMap().get(otherChromosome);
            if(chrGenes == null || !chrGenes.stream().anyMatch(x -> positionWithin(otherPosition, x.GeneStart, x.GeneEnd)))
                return true;
        }

        return false;
    }

    private List<EnsemblGeneData> findGeneCollection(final String chromosome, int geneCollectionId)
    {
        final Map<Integer,List<EnsemblGeneData>> geneCollectionMap = mChrGeneCollectionMap.get(chromosome);
        return geneCollectionMap != null && geneCollectionId >= 0 ? geneCollectionMap.get(geneCollectionId) : Lists.newArrayList();
    }

    public static void mergeChimericReadMaps(final Map<String,List<ReadRecord>> destMap, final Map<String,List<ReadRecord>> sourceMap)
    {
        for(Map.Entry<String,List<ReadRecord>> entry :  sourceMap.entrySet())
        {
            List<ReadRecord> readsById = destMap.get(entry.getKey());

            if(readsById == null)
            {
                destMap.put(entry.getKey(), entry.getValue());
            }
            else
            {
                readsById.addAll(entry.getValue());
            }
        }
    }

    private void logPerformanceStats()
    {
        if(!mConfig.Fusions.PerformanceStats)
            return;

        if(!ISF_LOGGER.isDebugEnabled() && mConfig.Functions.size() > 1)
            return;

        final PerformanceCounter[] perfCounters = mFusionTasks.get(0).getPerfCounters();

        for (int i = 1; i < mFusionTasks.size(); ++i)
        {
            final PerformanceCounter[] taskPCs = mFusionTasks.get(i).getPerfCounters();

            for (int j = 0; j < perfCounters.length; ++j)
            {
                perfCounters[j].merge(taskPCs[j]);
            }
        }

        Arrays.stream(perfCounters).forEach(x -> x.logStats());
    }

    @VisibleForTesting
    public final Map<String,List<FusionReadData>> getFusionCandidates()
    {
        return mFusionTasks.isEmpty() ? Maps.newHashMap() : mFusionTasks.get(0).getFusionCandidates();
    }

    public final Map<String,List<FusionFragment>> getUnfusedFragments()
    {
        return mFusionTasks.isEmpty() ? Maps.newHashMap() : mFusionTasks.get(0).getUnfusedFragments();
    }

    public void clearState()
    {
        mReadsMap.clear();
        mFusionTasks.get(0).clearState();
    }

}
