package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.NONE;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.getHighestMatchType;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.BOTH_JUNCTIONS;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.FS_DOWNSTREAM;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.FS_UPSTREAM;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.fusionId;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.hasTranscriptExonMatch;

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

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final Map<String,List<ReadRecord>> mReadsMap;
    private final Map<String,Map<Integer,List<EnsemblGeneData>>> mChrGeneCollectionMap;

    private int mNextFusionId;
    private final Map<String,List<FusionReadData>> mFusionCandidates; // keyed by the chromosome pair
    private final Map<String,List<FusionFragment>> mUnfusedFragments;

    private BufferedWriter mReadWriter;
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

        mPerfCounter = new PerformanceCounter("Fusions");
        mReadWriter = null;
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

            if(mConfig.WriteChimericReads)
                writeReadData(reads, readGroupStatus);
        }

        ISF_LOGGER.info("chimeric fragments({} unpaired={} dups={} skip={} filtered={} unspliced={} junc={}) fusions(loc={} total={})",
                mReadsMap.size(), unpairedReads, skipped, duplicates, filteredFragments,
                mUnfusedFragments.values().stream().count(), junctioned,
                mFusionCandidates.size(), mFusionCandidates.values().stream().mapToInt(x -> x.size()).sum());

        mReadsMap.clear();

        // classify / analyse fusions
        for(List<FusionReadData> fusions : mFusionCandidates.values())
        {
            fusions.forEach(x -> setGeneData(x));
        }

        markRelatedFusions();

        // assign any discordant reads
        assignUnfusedFragments();

        // write results
        writeFusionData();

        if(mConfig.WriteChimericReads)
            writeUnfilteredFragments();

        mPerfCounter.stop();
        mPerfCounter.logStats();

        closeBufferedWriter(mReadWriter);
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
            if(fusionData.spliceJunctionMatch(fragment))
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
        final FusionFragment fragment = fusionData.getFragments().get(0);

        final List<List<EnsemblGeneData>> genesByPosition = Lists.newArrayList(Lists.newArrayList(), Lists.newArrayList());

        for(int se = SE_START; se <= SE_END; ++se)
        {
            /*
            // first extract gene info from any matched exons
            final int seIndex = se;
            fusionData.getFragments().forEach(x -> x.setSplicedTransExonRefs(seIndex));

            if(fragment.getTransExonRefs().get(se).isEmpty())
            {
                // for intronic / unspliced junctions, use orientation to find the next splice acceptor or donor
                final List<EnsemblGeneData> geneDataList = findGeneCollection(fragment.chromosomes()[se], fragment.geneCollections()[se]);

                if(geneDataList != null && !geneDataList.isEmpty())
                {
                    final List<TranscriptData> transDataList = Lists.newArrayList();
                    geneDataList.forEach(x -> transDataList.addAll(mGeneTransCache.getTranscripts(x.GeneId)));

                    fusionData.getFragments().forEach(x -> x.populateUnsplicedTransExonRefs(transDataList, seIndex));
                }
            }
            */

            final List<String> spliceGeneIds = fragment.getGeneIds(se);

            if(!spliceGeneIds.isEmpty())
            {
                genesByPosition.set(se, spliceGeneIds.stream()
                        .map(x -> mGeneTransCache.getGeneDataById(x)).collect(Collectors.toList()));
            }
        }

        // organise genes by strand based on the orientations around the splice junction
        // a positive orientation implies either an upstream +ve strand gene or a downstream -ve strand gene
        final byte[] sjOrientations = fusionData.spliceOrientations();

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

        fusionData.cacheTranscriptData();
    }

    private void markRelatedFusions()
    {
        for( Map.Entry<String,List<FusionReadData>> entry : mFusionCandidates.entrySet())
        {
            final String locationId = entry.getKey();
            final List<FusionReadData> fusions = entry.getValue();

            if(fusions.size() == 1)
                continue;

            // annotate similar fusions for post-run analysis
            for(int i = 0; i < fusions.size() - 1; ++i)
            {
                FusionReadData fusion1 = fusions.get(i);
                boolean isSpliced = fusion1.hasSplicedFragments();
                boolean isUnspliced = fusion1.hasUnsplicedFragments();

                for(int j = i + 1; j < fusions.size() - 1; ++j)
                {
                    FusionReadData fusion2 = fusions.get(j);

                    if(isSpliced == fusion2.hasUnsplicedFragments() || isUnspliced == fusion2.hasSplicedFragments())
                    {
                        if(hasTranscriptExonMatch(fusion1.getTransExonRefsByStream(FS_UPSTREAM), fusion2.getTransExonRefsByStream(FS_UPSTREAM))
                        && hasTranscriptExonMatch(fusion1.getTransExonRefsByStream(FS_DOWNSTREAM), fusion2.getTransExonRefsByStream(FS_DOWNSTREAM)))
                        {
                            fusion1.addRelatedFusion(fusion2.id());
                            fusion2.addRelatedFusion(fusion1.id());
                            continue;
                        }
                    }

                    boolean isSimilar = false;

                    for(int se = SE_START; se <= SE_END; ++se)
                    {
                        if (abs(fusion1.splicePositions()[se] - fusion2.splicePositions()[se]) <= POSITION_REALIGN_DISTANCE)
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

    private static final int POSITION_REALIGN_DISTANCE = 20;

    private void assignUnfusedFragments()
    {
        for(Map.Entry<String,List<FusionFragment>> entry : mUnfusedFragments.entrySet())
        {
            final List<FusionReadData> fusions = mFusionCandidates.get(entry.getKey());

            if(fusions == null)
                continue;

            final List<FusionFragment> fragments = entry.getValue();

            /*
            for(int se = SE_START; se <= SE_END; ++se)
            {
                int geneCollectionId = fragments.get(0).geneCollections()[se];

                for (FusionFragment fragment : entry.getValue())
                {
                    final int seIndex = se;
                    fragment.populateDiscordantTransExonRefs(geneCollectionId, seIndex);
                }
            }
            */

            final List<FusionFragment> allocatedFragments = Lists.newArrayList();

            for (FusionFragment fragment : fragments)
            {
                // discordant reads don't have assigned
                if(fragment.getTransExonRefs()[SE_START].isEmpty() || fragment.getTransExonRefs()[SE_END].isEmpty())
                    continue;

                // if the fusion has unspliced reads, then these ought to set bounds for an discordant reads

                for(FusionReadData fusion : fusions)
                {
                    if(!fusion.isValid())
                        continue;

                    if(fusion.canAddDiscordantFragment(fragment))
                    {
                        fusion.addFusionFragment(fragment);
                        allocatedFragments.add(fragment);
                        break;
                    }
                }
            }

            allocatedFragments.forEach(x -> fragments.remove(x));
        }
    }

    public static final String FUSION_FILE_ID = "fusions.csv";

    private void writeFusionData()
    {
        if(mConfig.OutputDir == null)
            return;

        try
        {
            final String outputFileName = mConfig.formOutputFile(FUSION_FILE_ID);

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write(FusionReadData.csvHeader());
            writer.newLine();

            for(Map.Entry<String,List<FusionReadData>> entry : mFusionCandidates.entrySet())
            {
                for (final FusionReadData fusion : entry.getValue())
                {
                    writer.write(fusion.toCsv());
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fusions file: {}", e.toString());
        }
    }

    private void writeUnfilteredFragments()
    {
        for(List<FusionFragment> fragments : mUnfusedFragments.values())
        {
            for(FusionFragment fragment : fragments)
            {
                writeReadData(fragment.getReads(), "UNFUSED");
            }
        }

        for(List<FusionReadData> fusions : mFusionCandidates.values())
        {
            for(FusionReadData fusion : fusions)
            {
                for(FusionFragment fragment : fusion.getFragments())
                {
                    writeReadData(fragment.getReads(), fusionId(fusion.id()));
                }
            }
        }
    }

    private void writeReadData(final List<ReadRecord> reads, final String groupStatus)
    {
        try
        {
            if(mReadWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile("chimeric_reads.csv");

                mReadWriter = createBufferedWriter(outputFileName, false);
                mReadWriter.write("ReadSetCount,ReadId,FusionGroup,Chromosome,PosStart,PosEnd,Cigar,InsertSize");
                mReadWriter.write(",Secondary,Supplementary,NegStrand,ProperPair,SuppAlign,TransExons,BestMatch,TransExonData");
                mReadWriter.newLine();
            }

            for(final ReadRecord read : reads)
            {
                mReadWriter.write(String.format("%s,%s,%s,%s,%d,%d,%s,%d",
                        reads.size(), read.Id, groupStatus, read.Chromosome,
                        read.PosStart, read.PosEnd, read.Cigar.toString(), read.fragmentInsertSize()));

                mReadWriter.write(String.format(",%s,%s,%s,%s,%s",
                        read.isSecondaryAlignment(), read.isSupplementaryAlignment(), read.isNegStrand(), read.isProperPair(),
                        read.getSuppAlignment() != null ? read.getSuppAlignment().replaceAll(",", ";") : "NONE"));

                // log the transcript exons affected, and the highest matching transcript
                String transExonData = "";

                int transExonRefCount = 0;
                RegionMatchType highestTransMatchType = getHighestMatchType(read.getTransExonRefs().keySet());

                if(highestTransMatchType != NONE)
                {
                    for (final TransExonRef transExonRef : read.getTransExonRefs().get(highestTransMatchType))
                    {
                        ++transExonRefCount;
                        transExonData = appendStr(transExonData, String.format("%s:%d",
                                transExonRef.TransName, transExonRef.ExonRank), ';');
                    }
                }

                mReadWriter.write(String.format(",%d,%s,%s",
                        transExonRefCount, highestTransMatchType, transExonRefCount == 0 ? "NONE" : transExonData));
                mReadWriter.newLine();
            }

        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric read data: {}", e.toString());
            return;
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
