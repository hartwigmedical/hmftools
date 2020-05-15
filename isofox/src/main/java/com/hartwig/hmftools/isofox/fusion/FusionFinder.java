package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.abs;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.isofox.IsofoxConstants.ENRICHED_GENE_BUFFER;
import static com.hartwig.hmftools.isofox.common.ReadRecord.NO_GENE_ID;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsOverlap;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsWithin;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentBuilder.isValidFragment;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.FutureTask;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
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
import com.hartwig.hmftools.isofox.common.BamSlicer;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
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

    private final List<GenomeRegion> mRestrictedGeneRegions;
    private final List<GenomeRegion> mExcludedGeneRegions;

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
        mExcludedGeneRegions = Lists.newArrayList();
        mRestrictedGeneRegions = Lists.newArrayList();
        buildGeneRegions();

        mPerfCounter = new PerformanceCounter("Fusions");
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

    private static final int LOG_COUNT = 100000;

    private static final String LOG_READ_ID = "";
    // private static final String LOG_READ_ID = "NB500901:18:HTYNHBGX2:2:11303:18581:15228";

    public void findFusions()
    {
        ISF_LOGGER.info("processing {} chimeric read groups", mReadsMap.size());

        mPerfCounter.start();

        int invalidFragments = 0;
        int duplicates = 0;
        int skipped = 0;
        int fragments = 0;
        int missingSuppReads = 0;
        int recoveredSuppReads = 0;
        int readGroupCount = 0;
        int nextLog = LOG_COUNT;
        final boolean[] recoveryStatus = new boolean[SKIP+1];

        final Map<String,List<FusionFragment>> chrPairFragments = Maps.newHashMap();

        for(Map.Entry<String,List<ReadRecord>> entry : mReadsMap.entrySet())
        {
            ++readGroupCount;

            if(readGroupCount >= nextLog)
            {
                nextLog += LOG_COUNT;
                ISF_LOGGER.info("processed {} chimeric read groups", readGroupCount);
            }

            final List<ReadRecord> reads = entry.getValue();

            if(reads.stream().anyMatch(x -> x.isDuplicate()))
            {
                ++duplicates;
                continue;
            }

            if(reads.stream().anyMatch(x -> skipRead(x.mateChromosome(), x.mateStartPosition())))
            {
                ++skipped;
                continue;
            }

            // chimeric supplementary reads may be missed - attempt to find them
            recoverSupplementaryReads(reads, recoveryStatus);

            if(recoveryStatus[SKIP])
            {
                ++skipped;
                continue;
            }

            if(recoveryStatus[RECOVERED])
                ++recoveredSuppReads;
            else if(recoveryStatus[MISSING])
                ++missingSuppReads;

            if(reads.get(0).Id.equals(LOG_READ_ID))
            {
                ISF_LOGGER.debug("specific read: {}", reads.get(0));
            }

            if(!isValidFragment(reads))
            {
                ++invalidFragments;
                mFusionWriter.writeReadData(reads, "INVALID_READ_COUNT");
            }
            else
            {
                reads.forEach(x -> checkMissingGeneData(x));

                FusionFragment fragment = new FusionFragment(reads);

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

        ISF_LOGGER.info("chimeric groups({} skipped={} dups={} invalid={} recov={} miss={} candidates={}) chrPairs({}) tasks({})",
                mReadsMap.size(), skipped, duplicates, invalidFragments, recoveredSuppReads, missingSuppReads, fragments,
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

    private static final int RECOVERED = 0;
    private static final int MISSING = 1;
    private static final int SKIP = 2;

    private void recoverSupplementaryReads(final List<ReadRecord> reads, final boolean[] status)
    {
        status[RECOVERED] = false;
        status[MISSING] = false;
        status[SKIP] = false;

        if(reads.size() >= 3 || !reads.stream().anyMatch(x -> x.hasSuppAlignment()))
            return;

        for(final ReadRecord read : reads)
        {
            if(!read.hasSuppAlignment())
                continue;

            SupplementaryReadData suppData = SupplementaryReadData.from(read.getSuppAlignment());

            if(skipRead(suppData.Chromosome, suppData.Position))
            {
                status[SKIP] = true;
                return;
            }

            if(!mConfig.Fusions.RecoverMissingReads)
            {
                ISF_LOGGER.info("read({}) missing supp({}), recovery disabled", read, suppData);
                status[MISSING] = true;
                return;
            }

            final List<ReadRecord> missingReads = findMissingReads(read.Id, suppData.Chromosome, suppData.Position, reads.size() == 1);

            if(!missingReads.isEmpty())
            {
                ISF_LOGGER.info("id({}) recovered {} missing reads", read.Id, missingReads.size());

                if(ISF_LOGGER.isDebugEnabled())
                {
                    missingReads.stream().forEach(x -> ISF_LOGGER.info("recovered read: {}", x));
                }

                reads.addAll(missingReads);
                status[RECOVERED] = true;
            }
            else
            {
                ISF_LOGGER.info("read({}) missing supp({}) not found", read, suppData);
                status[MISSING] = true;
            }

            break;
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
