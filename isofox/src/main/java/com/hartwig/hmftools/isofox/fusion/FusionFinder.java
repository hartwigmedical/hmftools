package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.isofox.common.ReadRecord.generateMappedCoords;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.NONE;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.exonBoundary;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.getHighestMatchType;
import static com.hartwig.hmftools.isofox.common.RnaUtils.deriveCommonRegions;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.BOTH_JUNCTIONS;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.ONE_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGNED;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.FS_DOWNSTREAM;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.FS_UPSTREAM;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.fusionId;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.hasTranscriptExonMatch;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.hasTranscriptNextExonMatch;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BamSlicer;
import com.hartwig.hmftools.isofox.common.FragmentTracker;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class FusionFinder
{
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final Map<String,List<ReadRecord>> mReadsMap;
    private final Map<String,Map<Integer,List<EnsemblGeneData>>> mChrGeneCollectionMap;

    private int mNextFusionId;
    private final Map<String,List<FusionReadData>> mFusionCandidates; // keyed by the chromosome pair
    private final Map<String,List<FusionFragment>> mUnfusedFragments;

    private final FragmentTracker mReadDepthTracker;
    private final List<FusionReadData> mReadDepthFusions;
    private final FusionWriter mFusionWriter;

    private final PerformanceCounter mPerfCounter;

    public FusionFinder(final IsofoxConfig config, final EnsemblDataCache geneTransCache)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, true, true);

        mNextFusionId = 0;
        mReadsMap = Maps.newHashMap();
        mFusionCandidates = Maps.newHashMap();
        mChrGeneCollectionMap = Maps.newHashMap();

        mUnfusedFragments = Maps.newHashMap();

        mReadDepthTracker = new FragmentTracker();
        mReadDepthFusions = Lists.newArrayList();

        mPerfCounter = new PerformanceCounter("Fusions");
        mFusionWriter = new FusionWriter(mConfig);
    }

    public final Map<String,List<FusionReadData>> getFusionCandidates() { return mFusionCandidates; }

    public void clearState()
    {
        mNextFusionId = 0;
        mFusionCandidates.clear();
        mUnfusedFragments.clear();
        mReadsMap.clear();
    }

    public void addChimericReads(final Map<String,List<ReadRecord>> chimericReadMap)
    {
        mergeChimericReadMaps(mReadsMap, chimericReadMap);
    }

    public static void addChimericReads(final Map<String,List<ReadRecord>> chimericReadMap, final ReadRecord read)
    {
        List<ReadRecord> chimericReads = chimericReadMap.get(read.Id);
        if (chimericReads == null)
        {
            chimericReads = Lists.newArrayList();
            chimericReadMap.put(read.Id, chimericReads);
        }

        chimericReads.add(read);
    }

    public void addChromosomeGeneCollections(final String chromosome, final Map<Integer,List<EnsemblGeneData>> geneCollectionMap)
    {
        mChrGeneCollectionMap.put(chromosome, geneCollectionMap);
    }

    public void findFusions()
    {
        ISF_LOGGER.info("processing {} chimeric read groups", mReadsMap.size());

        mPerfCounter.start();

        int unpairedReads = 0;
        int filteredFragments = 0;
        int duplicates = 0;
        int skipped = 0;
        int junctioned = 0;

        for(Map.Entry<String,List<ReadRecord>> entry : mReadsMap.entrySet())
        {
            final List<ReadRecord> reads = entry.getValue();

            if(reads.stream().anyMatch(x -> x.isDuplicate()))
            {
                ++duplicates;
                continue;
            }

            String readGroupStatus = "";

            if(reads.size() == 1)
            {
                if(skipUnpairedRead(reads.get(0)))
                {
                    ++skipped;
                    continue;
                }

                reads.addAll(findMissingReads(reads.get(0)));

                if(reads.size() > 1)
                {
                    ISF_LOGGER.warn("read({}) found missing reads", reads.get(0).Id);
                }
            }

            if(reads.size() == 1)
            {
                ++unpairedReads;
                readGroupStatus = "UNPAIRED";
            }
            else if(isInvalidFragment(reads))
            {
                ++filteredFragments;
                readGroupStatus = "SINGLE_GENE";
            }
            else
            {
                FusionFragment fragment = new FusionFragment(reads);

                if(fragment.type() == BOTH_JUNCTIONS)
                {
                    ++junctioned;
                    createOrUpdateFusion(fragment);
                }
                else
                {
                    cacheUnfusedFragment(fragment);
                }

                // will write after further evaluation of fusions
                continue;
            }

            mFusionWriter.writeReadData(reads, readGroupStatus);
        }

        ISF_LOGGER.info("chimeric fragments({} unpaired={} dups={} skip={} filtered={} unfused={} junc={}) fusions(loc={} total={})",
                mReadsMap.size(), unpairedReads, skipped, duplicates, filteredFragments,
                mUnfusedFragments.values().stream().mapToInt(x -> x.size()).sum(), junctioned,
                mFusionCandidates.size(), mFusionCandidates.values().stream().mapToInt(x -> x.size()).sum());

        // mReadsMap.clear();

        // classify / analyse fusions
        for(List<FusionReadData> fusions : mFusionCandidates.values())
        {
            fusions.forEach(x -> setGeneData(x));
            calcFusionReadDepth(fusions);
        }

        markRelatedFusions();

        // assign any discordant reads
        assignUnfusedFragments();

        // write results
        mFusionWriter.writeFusionData(mFusionCandidates);
        mFusionWriter.writeUnfusedFragments(mUnfusedFragments);

        mPerfCounter.stop();
        mPerfCounter.logStats();
    }

    private boolean isInvalidFragment(final List<ReadRecord> reads)
    {
        final Set<String> chrGeneCollections = Sets.newHashSet();
        reads.forEach(x -> chrGeneCollections.add(String.format("%s_%d", x.Chromosome, x.getGeneCollecton())));

        if(chrGeneCollections.size() != 2)
            return true;

        return false;
    }

    private void cacheUnfusedFragment(final FusionFragment fragment)
    {
        List<FusionFragment> fragments = mUnfusedFragments.get(fragment.locationPair());

        if(fragments == null)
        {
            fragments = Lists.newArrayList();
            mUnfusedFragments.put(fragment.locationPair(), fragments);
        }

        fragments.add(fragment);
    }

    private boolean skipUnpairedRead(final ReadRecord read)
    {
        if(!HumanChromosome.contains(read.mateChromosome()))
            return true;

        if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(read.mateChromosome()))
            return true;

        if(!mConfig.RestrictedGeneIds.isEmpty())
        {
            final List<EnsemblGeneData> chrGenes = mGeneTransCache.getChrGeneDataMap().get(read.mateChromosome());
            if(chrGenes == null || !chrGenes.stream().anyMatch(x -> positionWithin(read.mateStartPosition(), x.GeneStart, x.GeneEnd)))
                return true;
        }

        return false;
    }

    private List<ReadRecord> findMissingReads(final ReadRecord read)
    {
        int chrSeqIndex = mSamReader.getFileHeader().getSequenceIndex(read.mateChromosome());

        QueryInterval[] queryInterval = new QueryInterval[1];
        queryInterval[0] = new QueryInterval(chrSeqIndex, (int)read.mateStartPosition(), (int)read.mateStartPosition());

        List<SAMRecord> records = mBamSlicer.slice(mSamReader, queryInterval);

        return records.stream().filter(x -> x.getReadName().equals(read.Id)).map(x -> ReadRecord.from(x)).collect(Collectors.toList());
    }

    private void createOrUpdateFusion(final FusionFragment fragment)
    {
        // scenarios:
        // 1. New fusion with correct splice-junction support - may or may not match a known transcript and exon
        //  -
        // 2. Additional fragment supporting the same junction
        // 3. Potential discordant read
        // 4. Invalid fragments for various reasons

        // fusions will be stored in a map keyed by their location pair (chromosome + geneCollectionId)
        List<FusionReadData> fusions = mFusionCandidates.get(fragment.locationPair());

        if(fusions == null)
        {
            fusions = Lists.newArrayList();
            mFusionCandidates.put(fragment.locationPair(), fusions);
        }

        for(final FusionReadData fusionData : fusions)
        {
            if(fusionData.junctionMatch(fragment))
            {
                fusionData.addFusionFragment(fragment);
                return;
            }
        }

        final FusionReadData fusionData = new FusionReadData(mNextFusionId++, fragment);
        fusions.add(fusionData);
    }

    private List<EnsemblGeneData> findGeneCollection(final String chromosome, int geneCollectionId)
    {
        final Map<Integer,List<EnsemblGeneData>> geneCollectionMap = mChrGeneCollectionMap.get(chromosome);
        return geneCollectionMap != null && geneCollectionId >= 0 ? geneCollectionMap.get(geneCollectionId) : Lists.newArrayList();
    }

    private void setGeneData(final FusionReadData fusionData)
    {
        // get the genes supporting the splice junction in the terms of an SV (ie lower chromosome and lower position first)
        final List<List<EnsemblGeneData>> genesByPosition = Lists.newArrayList(Lists.newArrayList(), Lists.newArrayList());

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final List<String> spliceGeneIds = fusionData.getSampleFragment().getGeneIds(se);

            if(!spliceGeneIds.isEmpty())
            {
                genesByPosition.set(se, spliceGeneIds.stream()
                        .map(x -> mGeneTransCache.getGeneDataById(x)).collect(Collectors.toList()));
            }

            // purge any invalid transcript exons when the splice junction is known
            final List<TranscriptData> transDataList = Lists.newArrayList();
            spliceGeneIds.forEach(x -> transDataList.addAll(mGeneTransCache.getTranscripts(x)));

            for(FusionFragment fragment : fusionData.getAllFragments())
            {
                fragment.validateTranscriptExons(transDataList, se);
            }
        }

        // organise genes by strand based on the orientations around the splice junction
        // a positive orientation implies either an upstream +ve strand gene or a downstream -ve strand gene
        final byte[] sjOrientations = fusionData.junctionOrientations();

        boolean foundCandidates = false;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            int upstreamIndex = se;
            int downstreamIndex = switchIndex(se);

            // for the start being the upstream gene, if the orientation is +1 then require a +ve strand gene,
            // and so the end is the downstream gene and will set the orientation as required

            final List<EnsemblGeneData> upstreamGenes = genesByPosition.get(upstreamIndex).stream()
                    .filter(x -> x.Strand == sjOrientations[upstreamIndex]).collect(Collectors.toList());

            final List<EnsemblGeneData> downstreamGenes = genesByPosition.get(downstreamIndex).stream()
                    .filter(x -> x.Strand == -sjOrientations[downstreamIndex]).collect(Collectors.toList());

            if(!upstreamGenes.isEmpty() && !downstreamGenes.isEmpty())
            {
                if(foundCandidates)
                {
                    // both combinations have possible gene-pairings
                    ISF_LOGGER.warn("fusion({}) has multiple gene pairings by strand and orientation", fusionData.toString());
                    fusionData.setIncompleteData();
                    break;
                }

                foundCandidates = true;
                fusionData.setStreamData(upstreamGenes, downstreamGenes, se == SE_START);
            }
        }

        // mark donor-acceptor types whether strands are known or not
        final byte[] geneStrands = fusionData.getGeneStrands();

        for(FusionFragment fragment : fusionData.getAllFragments())
        {
            fragment.setJunctionTypes(mConfig.RefFastaSeqFile, geneStrands) ;
        }

        fusionData.cacheTranscriptData();
    }

    private void markRelatedFusions()
    {
        for( Map.Entry<String,List<FusionReadData>> entry : mFusionCandidates.entrySet())
        {
            // final String locationId = entry.getKey();
            final List<FusionReadData> fusions = entry.getValue();

            if(fusions.size() == 1)
                continue;

            // annotate similar fusions for post-run analysis
            for(int i = 0; i < fusions.size() - 1; ++i)
            {
                FusionReadData fusion1 = fusions.get(i);

                // firstly the special case of a spliced fusion matching its unspliced fusion
                boolean isSpliced = fusion1.isKnownSpliced();
                boolean isUnspliced = fusion1.isUnspliced();

                for(int j = i + 1; j < fusions.size() - 1; ++j)
                {
                    FusionReadData fusion2 = fusions.get(j);

                    if(isSpliced && fusion2.isUnspliced())
                    {
                        if(hasTranscriptExonMatch(fusion1.getTransExonRefsByStream(FS_UPSTREAM), fusion2.getTransExonRefsByStream(FS_UPSTREAM))
                        && hasTranscriptNextExonMatch(fusion2.getTransExonRefsByStream(FS_DOWNSTREAM), fusion1.getTransExonRefsByStream(FS_DOWNSTREAM)))
                        {
                            fusion1.addRelatedFusion(fusion2.id());
                            fusion2.addRelatedFusion(fusion1.id());
                            continue;
                        }
                    }
                    else if(isUnspliced && fusion2.isKnownSpliced())
                    {
                        if(hasTranscriptExonMatch(fusion1.getTransExonRefsByStream(FS_UPSTREAM), fusion2.getTransExonRefsByStream(FS_UPSTREAM))
                        && hasTranscriptNextExonMatch(fusion1.getTransExonRefsByStream(FS_DOWNSTREAM), fusion2.getTransExonRefsByStream(FS_DOWNSTREAM)))
                        {
                            fusion1.addRelatedFusion(fusion2.id());
                            fusion2.addRelatedFusion(fusion1.id());
                            continue;
                        }
                    }

                    boolean isSimilar = false;

                    for(int se = SE_START; se <= SE_END; ++se)
                    {
                        if (abs(fusion1.junctionPositions()[se] - fusion2.junctionPositions()[se]) <= POSITION_REALIGN_DISTANCE)
                        {
                            isSimilar = true;
                            break;
                        }
                    }

                    if(isSimilar)
                    {
                        fusion1.addRelatedFusion(fusion2.id());
                        fusion2.addRelatedFusion(fusion1.id());
                    }
                }
            }
        }
    }

    private void calcFusionReadDepth(final List<FusionReadData> fusions)
    {
        mReadDepthTracker.clear();
        mReadDepthFusions.clear();
        mReadDepthFusions.addAll(fusions);

        QueryInterval[] queryInterval = new QueryInterval[SE_PAIR];

        for(int se = SE_START; se <= SE_END; ++se)
        {
            int[] fusionJuncRange = {-1, -1};

            for(final FusionReadData fusionData : fusions)
            {
                fusionJuncRange[SE_START] = fusionJuncRange[se] == -1 ?
                        fusionData.junctionPositions()[se] : min(fusionJuncRange[SE_START], fusionData.junctionPositions()[se]);

                fusionJuncRange[SE_END] = max(fusionJuncRange[SE_END], fusionData.junctionPositions()[se]);

                // set depth from existing fragments
            }

            int chrSeqIndex = mSamReader.getFileHeader().getSequenceIndex(fusions.get(0).chromosomes()[se]);
            queryInterval[se] = new QueryInterval(chrSeqIndex, fusionJuncRange[SE_START], fusionJuncRange[SE_END]);
        }

        mBamSlicer.slice(mSamReader, queryInterval, this::processReadDepthRecord);
    }

    private final static int SOFT_CLIP_JUNCTION_DISTANCE = 5;

    private void processReadDepthRecord(@NotNull final SAMRecord record)
    {
        if(mReadsMap.containsKey(record.getReadName()))
            return;

        final ReadRecord read1 = ReadRecord.from(record);
        final ReadRecord read2 = mReadDepthTracker.checkRead(read1);

        if(read2 == null)
            return;

        if(!read1.containsSoftClipping() && !read2.containsSoftClipping())
            return;

        if(!read1.Chromosome.equals(read2.Chromosome))
        {
            ISF_LOGGER.warn("read-depth read({}) spans chromosomes({} & {})", read1.Id, read1.Chromosome, read2.Chromosome);
            return;
        }

        // test these pairs against each fusion junction
        FusionFragment fragment = new FusionFragment(Lists.newArrayList(read1, read2));

        for(final FusionReadData fusionData : mReadDepthFusions)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                boolean hasReadSupport = false;
                boolean addedFragment = false;

                for(ReadRecord read : fragment.getReads())
                {
                    final List<int[]> readCoords = read.getMappedRegionCoords();
                    final int[] readBoundaries = { readCoords.get(0)[SE_START], readCoords.get(readCoords.size() - 1)[SE_END] };

                    if(!read.Chromosome.equals(fusionData.chromosomes()[se]))
                        continue;

                    final int seIndex = se;

                    if(!hasReadSupport && readCoords.stream()
                            .anyMatch(x -> positionWithin(fusionData.junctionPositions()[seIndex], x[SE_START], x[SE_END])))
                    {
                        hasReadSupport = true;
                    }

                    if(fusionData.junctionOrientations()[se] == 1 && read.isSoftClipped(SE_END))
                    {
                        if(positionWithin(
                                fusionData.junctionPositions()[se],
                                readBoundaries[SE_END], readBoundaries[SE_END] + SOFT_CLIP_JUNCTION_DISTANCE))
                        {
                            fragment.setType(ONE_JUNCTION);
                            fragment.setGeneData(fusionData.geneCollections()[se], fusionData.getTransExonRefsByPos(se));
                            fusionData.addFusionFragment(fragment);
                            addedFragment = true;
                            break;
                        }
                    }
                    else if(fusionData.junctionOrientations()[se] == -1 && read.isSoftClipped(SE_START))
                    {
                        if(positionWithin(
                                fusionData.junctionPositions()[se],
                                readBoundaries[SE_START] - SOFT_CLIP_JUNCTION_DISTANCE, readBoundaries[SE_START]))
                        {
                            fragment.setType(ONE_JUNCTION);
                            fragment.setGeneData(fusionData.geneCollections()[se], fusionData.getTransExonRefsByPos(se));
                            fusionData.addFusionFragment(fragment);
                            addedFragment = true;
                            break;
                        }
                    }
                }

                if(hasReadSupport && !addedFragment)
                    ++fusionData.getReadDepth()[se];
            }
        }
    }

    private static final int POSITION_REALIGN_DISTANCE = 20;

    private void assignUnfusedFragments()
    {
        for(Map.Entry<String,List<FusionFragment>> entry : mUnfusedFragments.entrySet())
        {
            final List<FusionReadData> fusions = mFusionCandidates.get(entry.getKey());

            if(fusions == null)
                continue;

            final List<FusionFragment> fragments = entry.getValue();

            final List<FusionFragment> allocatedFragments = Lists.newArrayList();

            for (FusionFragment fragment : fragments)
            {
                if(fragment.getTransExonRefs()[SE_START].isEmpty() || fragment.getTransExonRefs()[SE_END].isEmpty())
                    continue;

                for(FusionReadData fusion : fusions)
                {
                    if(!fusion.isValid())
                        continue;

                    if(fusion.canAddDiscordantFragment(fragment))
                    {
                        for(int se = SE_START; se <= SE_END; ++se)
                        {
                            if (fragment.junctionValid()[se] && abs(fusion.junctionPositions()[se] - fragment.junctionPositions()[se]) <= POSITION_REALIGN_DISTANCE)
                            {
                                fragment.setType(REALIGNED);
                                break;
                            }
                        }

                        fusion.addFusionFragment(fragment);
                        allocatedFragments.add(fragment);
                    }
                }
            }

            allocatedFragments.forEach(x -> fragments.remove(x));
        }
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

}
