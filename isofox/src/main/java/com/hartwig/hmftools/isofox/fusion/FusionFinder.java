package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.ENRICHED_GENE_BUFFER;
import static com.hartwig.hmftools.isofox.common.ReadRecord.NO_GENE_ID;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsOverlap;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsWithin;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentBuilder.isValidFragment;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.fusion.ReadGroup.mergeChimericReadMaps;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.ThreadFactory;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionFinder
{
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final List<ReadGroup> mChimericReadGroups;
    private final Map<String,ReadGroup> mChimericPartialReadGroups;
    private final Set<String> mDuplicateReadIds;
    private final Map<String,Map<Integer,List<EnsemblGeneData>>> mChrGeneCollectionMap;
    private final Map<String,Map<Integer,BaseDepth>> mChrGeneDepthMap;

    private List<FusionTask> mFusionTasks;
    private final FusionWriter mFusionWriter;

    private final List<GenomeRegion> mRestrictedGeneRegions;
    private final List<GenomeRegion> mExcludedGeneRegions;

    private final PerformanceCounter mPerfCounter;

    public FusionFinder(final IsofoxConfig config, final EnsemblDataCache geneTransCache)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mChimericPartialReadGroups = Maps.newHashMap();
        mChimericReadGroups = Lists.newArrayList();
        mDuplicateReadIds = Sets.newHashSet();
        mChrGeneCollectionMap = Maps.newHashMap();
        mChrGeneDepthMap = Maps.newHashMap();
        mFusionTasks = Lists.newArrayList();

        mExcludedGeneRegions = Lists.newArrayList();
        mRestrictedGeneRegions = Lists.newArrayList();
        buildGeneRegions();

        mPerfCounter = new PerformanceCounter("Fusions");
        mFusionWriter = new FusionWriter(mConfig);
    }

    public void addChimericReads(final Map<String,ReadGroup> partialReadGroups, final List<ReadGroup> readGroups)
    {
        mergeChimericReadMaps(mChimericPartialReadGroups, mChimericReadGroups, partialReadGroups);
        mChimericReadGroups.addAll(readGroups);
    }

    public void addChimericReads(final List<ReadGroup> readGroups)
    {
        mChimericReadGroups.addAll(readGroups);
    }

    public void addDuplicateReadIds(final Set<String> readIds)
    {
        mergeDuplicateReadIds(mDuplicateReadIds, readIds);
    }

    public static void addChimericReads(final Map<String,ReadGroup> chimericReadMap, final ReadRecord read)
    {
        ReadGroup chimericReads = chimericReadMap.get(read.Id);
        if (chimericReads == null)
            chimericReadMap.put(read.Id, new ReadGroup(read));
        else
            chimericReads.Reads.add(read);
    }

    public void addChromosomeGeneCollections(final String chromosome, final Map<Integer,List<EnsemblGeneData>> geneCollectionMap)
    {
        mChrGeneCollectionMap.put(chromosome, geneCollectionMap);
    }

    public void addChromosomeGeneDepth(final String chromosome, final Map<Integer,BaseDepth> geneDepthMap)
    {
        mChrGeneDepthMap.put(chromosome, geneDepthMap);
    }

    private static final int LOG_COUNT = 100000;

    private static final String LOG_READ_ID = "";
    // private static final String LOG_READ_ID = "NB500901:18:HTYNHBGX2:2:23308:18394:18413";

    public void findFusions()
    {
        // convert any set of valid reads into a fragment, and then process these in groups by chromosomal pair
        ISF_LOGGER.info("processing {} chimeric read groups", mChimericReadGroups.size());

        mPerfCounter.start();

        int invalidFragments = 0;
        int hasMissingReads = 0;
        int hasExcesssReads = 0;
        int duplicates = 0;
        int skipped = 0;
        int fragments = 0;
        int missingSuppReads = 0;

        int partialDups = 0;
        int partialSkipped = 0;

        int readGroupCount = 0;
        int nextLog = LOG_COUNT;

        final Map<String,List<FusionFragment>> chrPairFragments = Maps.newHashMap();

        final List<ReadGroup> partialGroups = Lists.newArrayList();

        mChimericReadGroups.addAll(mChimericPartialReadGroups.values());

        for(ReadGroup readGroup : mChimericReadGroups)
        {
            ++readGroupCount;

            if(readGroupCount >= nextLog)
            {
                nextLog += LOG_COUNT;
                ISF_LOGGER.info("processed {} chimeric read groups", readGroupCount);
            }

            boolean isComplete = readGroup.isComplete();
            final List<ReadRecord> reads = readGroup.Reads;

            if(reads.get(0).Id.equals(LOG_READ_ID))
            {
                ISF_LOGGER.debug("specific read: {}", reads.get(0));
            }

            if(mDuplicateReadIds.contains(reads.get(0).Id) || reads.stream().anyMatch(x -> x.isDuplicate()))
            {
                ++duplicates;

                if(!isComplete)
                    ++partialDups;

                continue;
            }

            if(reads.stream().anyMatch(x -> skipRead(x.mateChromosome(), x.mateStartPosition())))
            {
                ++skipped;

                if(!isComplete)
                    ++partialSkipped;

                continue;
            }

            if(!isComplete)
            {
                if(readGroup.hasSuppAlignment())
                {
                    if(skipMissingReads(reads))
                    {
                        ++skipped;
                        continue;
                    }

                    ++missingSuppReads;
                }

                ++invalidFragments;

                if(reads.size() > 3 || (!readGroup.hasSuppAlignment() && reads.size() > 2))
                    ++hasExcesssReads;
                else
                    ++hasMissingReads;

                mFusionWriter.writeReadData(reads, "INVALID_READ_COUNT");

                partialGroups.add(readGroup);
            }
            else
            {
                if(mChimericPartialReadGroups.containsKey(readGroup.id()))
                {
                    ISF_LOGGER.error("partial read({}) group marked as complete", readGroup.id());
                    partialGroups.add(readGroup);

                    if(mChimericReadGroups.stream().anyMatch(x -> x.id().equals(readGroup.id())))
                    {
                        ISF_LOGGER.error("partial read({}) group also in complete list", readGroup.id());
                    }
                }

                reads.forEach(x -> checkMissingGeneData(x));

                FusionFragment fragment = new FusionFragment(readGroup);

                if(fragment.type() == FusionFragmentType.UNKNOWN)
                {
                    ++invalidFragments;
                    mFusionWriter.writeReadData(reads, "INVALID_FRAG");
                    continue;
                }

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
            }
        }

        // allocate fusion pairs evenly amongst threads (if multi-thread)
        mFusionTasks.clear();

        for(int taskId = 0; taskId < max(mConfig.Threads, 1); ++taskId)
        {
            mFusionTasks.add(new FusionTask(taskId, mConfig, mGeneTransCache, mChrGeneDepthMap, mFusionWriter));
        }

        for(List<FusionFragment> chrPairFrags : chrPairFragments.values())
        {
            // allocate the next chr-pair fragment batch to the task with the least
            FusionTask leastAllocated = null;
            int minAllocated = 0;

            for(FusionTask fusionTask : mFusionTasks)
            {
                if(minAllocated == 0 || fusionTask.getFragments().size() < minAllocated)
                {
                    leastAllocated = fusionTask;
                    minAllocated = fusionTask.getFragments().size();
                    if(minAllocated == 0)
                        break;
                }
            }

            leastAllocated.getFragments().addAll(chrPairFrags);
        }

        ISF_LOGGER.info("chimeric groups({} skipped={} dups=({} existing={}) invalid={} miss={} candidates={}) chrPairs({}) tasks({})",
                mChimericReadGroups.size(), skipped, duplicates, mDuplicateReadIds.size(), invalidFragments, missingSuppReads, fragments,
                chrPairFragments.size(), mFusionTasks.size());

        if(!mChimericPartialReadGroups.isEmpty())
        {
            int complete = (int)mChimericPartialReadGroups.values().stream().filter(x -> x.isComplete()).count();

            ISF_LOGGER.info("partial groups({} complete={} dups={} skip={} excess={} miss={})",
                    mChimericPartialReadGroups.size(), complete, partialDups, partialSkipped, hasExcesssReads, hasMissingReads);
        }

        mChimericPartialReadGroups.clear();
        mChimericReadGroups.clear();
        mDuplicateReadIds.clear();

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

    private boolean skipMissingReads(final List<ReadRecord> reads)
    {
        for(final ReadRecord read : reads)
        {
            if(read.hasSuppAlignment())
            {
                SupplementaryReadData suppData = SupplementaryReadData.from(read.getSuppAlignment());

                if(skipRead(suppData.Chromosome, suppData.Position))
                    return true;

                ISF_LOGGER.info("read({}) missing supp({})", read, suppData);
            }
        }

        return false;
    }

    private boolean skipRead(final String otherChromosome, int otherPosition)
    {
        if(!HumanChromosome.contains(otherChromosome))
            return true;

        if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(otherChromosome))
            return true;

        if(!mRestrictedGeneRegions.isEmpty() && !mRestrictedGeneRegions.stream().filter(x -> x.chromosome().equals(otherChromosome))
                .anyMatch(x -> positionWithin(otherPosition, (int)x.start(), (int)x.end())))
        {
            return true;
        }

        if(mExcludedGeneRegions.stream().filter(x -> x.chromosome().equals(otherChromosome))
                .anyMatch(x -> positionWithin(otherPosition, (int)x.start(), (int)x.end())))
        {
            return true;
        }

        return false;
    }

    private void setMissingGeneCollection(final ReadRecord read)
    {
        if(read.getGeneCollectons()[SE_START] != NO_GENE_ID && read.getGeneCollectons()[SE_END] == NO_GENE_ID)
        {
            final Map<Integer,List<EnsemblGeneData>> geneCollectionMap = mChrGeneCollectionMap.get(read.Chromosome);

            int readEnd = read.getCoordsBoundary(SE_END);

            for(Map.Entry<Integer,List<EnsemblGeneData>> entry : geneCollectionMap.entrySet())
            {
                final int[] genesRange = new int[] {
                        entry.getValue().stream().mapToInt(x -> x.GeneStart).min().orElse(0),
                        entry.getValue().stream().mapToInt(x -> x.GeneEnd).max().orElse(0) };

                if(positionWithin(readEnd, genesRange[SE_START], genesRange[SE_END]))
                {
                    read.setGeneCollection(SE_END, entry.getKey(), true);
                    return;
                }
                else if(readEnd < genesRange[SE_START])
                {
                    read.setGeneCollection(SE_END, entry.getKey(), false);
                    return;
                }
            }
        }
    }

    private void checkMissingGeneData(final ReadRecord read)
    {
        if(!read.spansGeneCollections())
            return;

        if(read.getGeneCollectons()[SE_END] == NO_GENE_ID)
        {
            setMissingGeneCollection(read);
        }

        if(!read.getIsGenicRegion()[SE_END])
            return;

        // due to the way the BAM fragment allocator processes reads per gene collection, the upper gene collection will have missed its
        // transcript exon data, so populate this now

        final List<EnsemblGeneData> geneDataList = findGeneCollection(read.Chromosome, read.getGeneCollectons()[SE_END]);

        if(geneDataList == null)
            return;

        int upperCoordIndex = read.getMappedRegionCoords().size() - 1;
        final int[] upperCoords = read.getMappedRegionCoords().get(upperCoordIndex);
        final Map<RegionMatchType,List<TransExonRef>> transExonRefMap = read.getTransExonRefs(SE_END);

        for(EnsemblGeneData geneData : geneDataList)
        {
            final List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);

            for(TranscriptData transData : transDataList)
            {
                if(!positionsWithin(upperCoords[SE_START], upperCoords[SE_END], transData.TransStart, transData.TransEnd))
                    continue;

                for(ExonData exonData : transData.exons())
                {
                    if(!positionsOverlap(upperCoords[SE_START], upperCoords[SE_END], exonData.ExonStart, exonData.ExonEnd))
                        continue;

                    RegionMatchType matchType;
                    if(upperCoords[SE_START] == exonData.ExonStart || upperCoords[SE_END] == exonData.ExonEnd)
                    {
                        matchType = RegionMatchType.EXON_BOUNDARY;
                    }
                    else if(positionsWithin(upperCoords[SE_START], upperCoords[SE_END], exonData.ExonStart, exonData.ExonEnd))
                    {
                        matchType = RegionMatchType.WITHIN_EXON;
                    }
                    else
                    {
                        matchType = RegionMatchType.EXON_INTRON;
                    }

                    TransExonRef teRef = new TransExonRef(transData.GeneId, transData.TransId, transData.TransName, exonData.ExonRank);

                    final List<TransExonRef> transExonRefs = transExonRefMap.get(matchType);

                    if(transExonRefs == null)
                        transExonRefMap.put(matchType, Lists.newArrayList(teRef));
                    else
                        transExonRefs.add(teRef);

                    break;
                }
            }
        }
    }

    private void buildGeneRegions()
    {
        mConfig.EnrichedGeneIds.stream()
                .map(x -> mGeneTransCache.getGeneDataById(x))
                .filter(x -> x != null)
                .forEach(x -> mExcludedGeneRegions.add(GenomeRegions.create(
                        x.Chromosome, x.GeneStart - ENRICHED_GENE_BUFFER, x.GeneEnd + ENRICHED_GENE_BUFFER)));

        mConfig.ExcludedGeneIds.stream()
                .map(x -> mGeneTransCache.getGeneDataById(x))
                .filter(x -> x != null)
                .forEach(x -> mExcludedGeneRegions.add(GenomeRegions.create(
                        x.Chromosome, x.GeneStart, x.GeneEnd)));

        mConfig.RestrictedGeneIds.stream()
                .map(x -> mGeneTransCache.getGeneDataById(x))
                .filter(x -> x != null)
                .forEach(x -> mRestrictedGeneRegions.add(GenomeRegions.create(
                        x.Chromosome, x.GeneStart - 1000, x.GeneEnd + 1000)));
    }

    private List<EnsemblGeneData> findGeneCollection(final String chromosome, int geneCollectionId)
    {
        final Map<Integer,List<EnsemblGeneData>> geneCollectionMap = mChrGeneCollectionMap.get(chromosome);
        return geneCollectionMap != null && geneCollectionId >= 0 ? geneCollectionMap.get(geneCollectionId) : Lists.newArrayList();
    }

    public static void mergeDuplicateReadIds(final Set<String> destSet, final Set<String> sourceSet)
    {
        sourceSet.forEach(x -> destSet.add(x));
    }

    private void logPerformanceStats()
    {
        if(!mConfig.Fusions.PerformanceStats)
            return;

        if(!ISF_LOGGER.isDebugEnabled() && mConfig.Functions.size() > 1)
            return;

        final PerformanceCounter perfCounter = mFusionTasks.get(0).getPerfCounter();

        for (int i = 1; i < mFusionTasks.size(); ++i)
        {
            final PerformanceCounter taskPC = mFusionTasks.get(i).getPerfCounter();
            perfCounter.merge(mFusionTasks.get(i).getPerfCounter());
        }

        perfCounter.logStats();
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
        mChimericPartialReadGroups.clear();
        mChimericReadGroups.clear();
        mFusionTasks.get(0).clearState();
    }

}
