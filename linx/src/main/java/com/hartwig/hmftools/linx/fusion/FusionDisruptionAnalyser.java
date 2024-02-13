package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.IG_PROMISCUOUS;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.NONE;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.fusion.DisruptionFinder.getUndisruptedCopyNumber;
import static com.hartwig.hmftools.linx.fusion.FusionConfig.WRITE_NEO_EPITOPES;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.FUSION_MAX_CHAIN_LENGTH;
import static com.hartwig.hmftools.linx.fusion.FusionFinder.validFusionTranscript;
import static com.hartwig.hmftools.linx.fusion.FusionReportability.allowSuspectChains;
import static com.hartwig.hmftools.linx.fusion.FusionReportability.findTopPriorityFusion;
import static com.hartwig.hmftools.linx.fusion.FusionReportability.isReportable;
import static com.hartwig.hmftools.linx.fusion.FusionWriter.convertBreakendsAndFusions;
import static com.hartwig.hmftools.linx.fusion.rna.RnaFusionMapper.RNA_FUSIONS_FILE;
import static com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType.FUSION;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.linx.gene.BreakendGeneData;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.linx.gene.BreakendTransData;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.CohortDataWriter;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.fusion.rna.RnaFusionMapper;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.visualiser.file.VisFusion;
import com.hartwig.hmftools.linx.visualiser.file.VisSampleData;

public class FusionDisruptionAnalyser
{
    private FusionFinder mFusionFinder;
    private DisruptionFinder mDisruptionFinder;
    private FusionWriter mFusionWriter;
    private NeoEpitopeWriter mNeoEpitopeWriter;

    private String mSampleId;
    private final String mOutputDir;
    private final EnsemblDataCache mGeneDataCache;
    private final LinxConfig mConfig;

    private final boolean mRunFusions;
    private final FusionConfig mFusionConfig;

    private final List<GeneFusion> mFusions; // all possible valid transcript-pair fusions
    private final List<GeneFusion> mUniqueFusions; // top-priority fusions from within each unique gene and SV pair
    private final Map<GeneFusion,String> mInvalidFusions;
    private final SpecialFusions mSpecialFusions;

    private RnaFusionMapper mRnaFusionMapper;
    private final VisSampleData mVisSampleData;

    private PerformanceCounter mPerfCounter;

    public FusionDisruptionAnalyser(
            final LinxConfig config, final EnsemblDataCache ensemblDataCache, final FusionResources fusionResources,
            final CohortDataWriter cohortDataWriter, final VisSampleData visSampleData)
    {
        mOutputDir = config.OutputDataPath;

        mConfig = config;
        mGeneDataCache = ensemblDataCache;

        mFusionConfig = new FusionConfig(config.CmdLineConfig);
        mRunFusions = mConfig.RunFusions;

        mFusionFinder = new FusionFinder(mFusionConfig, ensemblDataCache, fusionResources.knownFusionCache());
        mFusionWriter = new FusionWriter(mOutputDir, cohortDataWriter);
        mDisruptionFinder = new DisruptionFinder(config, ensemblDataCache, cohortDataWriter, visSampleData);
        mVisSampleData = visSampleData;

        mNeoEpitopeWriter = null;

        mFusions = Lists.newArrayList();
        mUniqueFusions = Lists.newArrayList();
        mInvalidFusions = Maps.newHashMap();

        mSpecialFusions = new SpecialFusions(mGeneDataCache, mFusionFinder, mFusionConfig);

        mRnaFusionMapper = null;

        mPerfCounter = new PerformanceCounter("Fusions");

        if(config.CmdLineConfig.hasValue(RNA_FUSIONS_FILE))
        {
            if(mConfig.Threads > 1 && mConfig.hasMultipleSamples())
            {
                LNX_LOGGER.error("RNA fusions mapping not yet supported in multi-threading mode");
            }
            else
            {
                mRnaFusionMapper = new RnaFusionMapper(
                        mOutputDir, config.CmdLineConfig, mGeneDataCache, mFusionFinder, mUniqueFusions, mInvalidFusions);
            }
        }

        if(config.CmdLineConfig.hasFlag(WRITE_NEO_EPITOPES))
        {
            mNeoEpitopeWriter = new NeoEpitopeWriter(mOutputDir, mGeneDataCache, mFusionFinder.getKnownFusionCache());
        }

        if(mRunFusions)
        {
            LNX_LOGGER.debug("fusion config: requirePhaseMatch({}) allowExonSkipping({}) requireUpstreamBiotypes({})",
                    mFusionConfig.RequirePhaseMatch, mFusionConfig.AllowExonSkipping, mFusionConfig.RequireUpstreamBiotypes);

            mSpecialFusions.cacheSpecialFusionGenes();
        }
    }

    public PerformanceCounter getPerfCounter() { return mPerfCounter; }

    public boolean hasRnaSampleData() { return mRnaFusionMapper != null; }
    public final List<GeneFusion> getFusions() { return mFusions; }
    public final List<GeneFusion> getUniqueFusions() { return mUniqueFusions; }
    public final Map<GeneFusion,String> getInvalidFusions() { return mInvalidFusions; }
    public final FusionFinder getFusionFinder() { return mFusionFinder; }
    public final DisruptionFinder getDisruptionFinder() { return mDisruptionFinder; }
    public final SpecialFusions getSpecialFusions() { return mSpecialFusions; }

    public void annotateTranscripts(final List<SvVarData> svList, boolean purgeInvalidTranscripts)
    {
        // mark any transcripts as not disruptive prior to running any fusion logic
        mDisruptionFinder.markTranscriptsDisruptive(svList);

        for(final SvVarData var : svList)
        {
            // now that transcripts have been marked as disruptive it is safe to purge any which cannot make viable fusions
            if(purgeInvalidTranscripts)
            {
                for(int be = SE_START; be <= SE_END; ++be)
                {
                    if(be == SE_END && var.isSglBreakend())
                        continue;

                    boolean isStart = isStart(be);

                    List<BreakendGeneData> genesList = var.getGenesList(isStart);

                    for(BreakendGeneData gene : genesList)
                    {
                        int transIndex = 0;
                        while (transIndex < gene.transcripts().size())
                        {
                            BreakendTransData transcript = gene.transcripts().get(transIndex);

                            // only retain transcript which are potential fusion candidates (with exception for canonical)
                            if(!transcript.isDisruptive() && !validFusionTranscript(transcript) && !transcript.isCanonical())
                            {
                                gene.transcripts().remove(transIndex);
                            }
                            else
                            {
                                ++transIndex;
                            }
                        }
                    }
                }
            }
        }
    }

    public void run(
            final String sampleId, final List<SvVarData> fullSvList, final List<SvCluster> clusters,
            final Map<String,List<SvBreakend>> chrBreakendMap)
    {
        mPerfCounter.start();

        mSampleId = sampleId;

        mUniqueFusions.clear();
        mFusionFinder.reset();

        final List<SvVarData> svList = fullSvList.stream()
                .filter(x -> x.getLinkedSVs() == null || !x.isSglBreakend()).collect(Collectors.toList());

        if(mConfig.IsGermline && mConfig.isSingleSample())
        {
            mDisruptionFinder.findReportableDisruptions(svList, clusters);
            mDisruptionFinder.writeGermlineDisruptions(sampleId, mConfig.OutputDataPath);
            return;
        }

        if(mRunFusions && mFusionFinder.hasValidConfigData())
            findFusions(svList, clusters);

        mDisruptionFinder.findReportableDisruptions(svList, clusters);

        if(mRunFusions)
        {
            mUniqueFusions.addAll(extractUniqueFusions());

            // add protein information which won't have been set for unreported fusions
            mUniqueFusions.stream().filter(x -> x.knownType() != NONE).forEach(x -> mFusionFinder.setFusionProteinFeatures(x));
        }

        final List<BreakendTransData> transcripts = buildTranscriptList(svList, mUniqueFusions);

        final List<LinxBreakend> breakends = Lists.newArrayList();
        final List<LinxFusion> fusions = Lists.newArrayList();
        convertBreakendsAndFusions(mUniqueFusions, transcripts, fusions, breakends);

        if(mConfig.isSingleSample())
        {
            mFusionWriter.writeSampleData(mSampleId, fusions, breakends);

            if(mFusionConfig.LogAllPotentials)
            {
                mFusions.forEach(x -> mFusionWriter.writeVerboseFusionData(x, mSampleId));
            }
        }
        else
        {
            // write fusions in detail when in batch mode
            final List<GeneFusion> fusionList = mFusionConfig.LogAllPotentials ? mFusions : mUniqueFusions;

            fusionList.stream()
                    .filter(x -> x.reportable() || !mFusionConfig.LogReportableOnly)
                    .forEach(x -> mFusionWriter.writeVerboseFusionData(x, mSampleId));
        }

        if(mConfig.hasMultipleSamples() || mConfig.Output.WriteCohortFiles)
        {
            mDisruptionFinder.writeCohortData(mSampleId, svList);
        }

        if(mRunFusions)
        {
            addVisualisationData(mUniqueFusions);

            if(LNX_LOGGER.isDebugEnabled())
            {
                for(final GeneFusion fusion : mUniqueFusions)
                {
                    if(fusion.knownType() != KnownFusionType.NONE)
                    {
                        LNX_LOGGER.debug("fusion({}:{}) reportable({}) knownType({}) cluster({} sv={} chain={}) SVs({} & {})",
                                fusion.id(), fusion.name(), fusion.reportable(),
                                fusion.knownTypeStr(), fusion.getAnnotations().clusterId(), fusion.getAnnotations().clusterCount(),
                                fusion.getAnnotations().chainInfo() != null ? fusion.getAnnotations().chainInfo().chainId() : -1,
                                fusion.upstreamTrans().gene().id(), fusion.downstreamTrans().gene().id());
                    }
                }
            }
        }

        if(mRnaFusionMapper != null)
            mRnaFusionMapper.assessRnaFusions(sampleId, chrBreakendMap);

        mPerfCounter.stop();
    }

    private void findFusions(final List<SvVarData> svList, final List<SvCluster> clusters)
    {
        if(mSampleId.isEmpty() || mFusionFinder == null)
            return;

        mFusions.clear();
        mInvalidFusions.clear();
        mSpecialFusions.clear();

        if(mNeoEpitopeWriter != null)
            mNeoEpitopeWriter.initialiseSample(mSampleId);

        findSingleSVFusions(svList);
        findChainedFusions(clusters);

        mFusions.addAll(mSpecialFusions.findFusions(svList));
    }

    private void findSingleSVFusions(final List<SvVarData> svList)
    {
        // always report SVs by themselves
        for(final SvVarData var : svList)
        {
            if(var.isSglBreakend() && var.getSglMappings().isEmpty())
                continue;

            // skip SVs which have been chained unless they are in an IG region
            if(var.getCluster().getSvCount() > 1 && var.getCluster().findChain(var) != null)
            {
                boolean igCandidate = !var.isSglBreakend()
                        && (mFusionFinder.getKnownFusionCache().withinIgRegion(var.chromosome(true), var.position(true))
                        != mFusionFinder.getKnownFusionCache().withinIgRegion(var.chromosome(false), var.position(false)));

                if(!igCandidate)
                    continue;
            }

            final List<BreakendGeneData> genesListStart = getBreakendGeneList(var, true);
            final List<BreakendGeneData> genesListEnd =  getBreakendGeneList(var, false);

            if(genesListStart.isEmpty() || genesListEnd.isEmpty())
                continue;

            List<GeneFusion> fusions = mFusionFinder.findFusions(genesListStart, genesListEnd);

            if(mNeoEpitopeWriter != null)
            {
                double svCopyNumber = var.isSglBreakend() ?
                    var.copyNumber(true) : (var.copyNumber(true) + var.copyNumber(false)) * 0.5;

                mNeoEpitopeWriter.processFusionCandidate(
                        genesListStart, genesListEnd, null, null, null, null,
                        svCopyNumber);
            }

            if(fusions.isEmpty())
                continue;

            if(mFusionConfig.LogReportableOnly)
            {
                fusions = fusions.stream().filter(x -> isReportable(x)).collect(Collectors.toList());
            }

            final SvCluster cluster = var.getCluster();

            // check transcript disruptions
            for(final GeneFusion fusion : fusions)
            {
                FusionAnnotations annotations = ImmutableFusionAnnotations.builder()
                        .clusterId(cluster.id())
                        .clusterCount(cluster.getSvCount())
                        .resolvedType(cluster.getResolvedType().toString())
                        .chainInfo(null)
                        .terminatedUp(false)
                        .terminatedDown(false)
                        .build();

                fusion.setAnnotations(annotations);
                mFusions.add(fusion);
            }
        }
    }

    private void findChainedFusions(final List<SvCluster> clusters)
    {
        // for now only consider simple SVs and resolved small clusters
        for(final SvCluster cluster : clusters)
        {
            if(cluster.getSvCount() == 1) // simple clusters already checked
                continue;

            if(cluster.getChains().isEmpty())
                continue;

            final List<GeneFusion> chainFusions = Lists.newArrayList();
            final List<ValidTraversalData> validPairs = Lists.newArrayList();

            for(final SvChain chain : cluster.getChains())
            {
                findChainedFusions(cluster, chain, chainFusions, validPairs);
            }

            if(chainFusions.isEmpty())
                continue;

            // now all fusions have been gathered from this chain, set the reportable one (if any)
            LNX_LOGGER.trace("cluster({}) found {} chained fusions", cluster.id(), chainFusions.size());

            // consider fusions from amongst unique gene-pairings
            final Map<String,List<GeneFusion>> genePairFusions = Maps.newHashMap();

            // allocate fusions to unique SV pairings
            for(GeneFusion fusion : chainFusions)
            {
                final String name = fusion.name();

                List<GeneFusion> fusions = genePairFusions.get(name);

                if(fusions == null)
                    genePairFusions.put(name, Lists.newArrayList(fusion));
                else
                    fusions.add(fusion);
            }

            mFusions.addAll(chainFusions.stream()
                    .filter(x -> !mFusionConfig.LogReportableOnly || x.reportable())
                    .filter(x -> x.getAnnotations() != null)
                    .collect(Collectors.toList()));
        }
    }

    private void findChainedFusions(
            final SvCluster cluster, final SvChain chain, final List<GeneFusion> chainFusions, final List<ValidTraversalData> validPairs)
    {
        // look for fusions formed by breakends connected in a chain

        // given a chain sAe - sBe - sCe - sDe, where As and De are the open ends of the chain, the following fusions need to be tested:
        // a) each SV in isolation, ie as a single SV
        // b) the left-most / lower breakend of each SV with the right-most / upper breakend of SVs higher up in the chain

        // whenever a linked pair is traversed by a fusion, it cannot touch or traverse genic regions without disrupting the fusion
        final List<LinkedPair> linkedPairs = chain.getLinkedPairs();

        for(int lpIndex1 = 0; lpIndex1 <= linkedPairs.size(); ++lpIndex1)
        {
            SvVarData lowerSV = null;
            SvBreakend lowerBreakend = null;

            // the lower link takes the other breakend of the current linked pair's 'first' SV
            // and in order for it to also test the last SV in isolation, also takes the last pair's second (upper) breakend
            if(lpIndex1 < linkedPairs.size())
            {
                LinkedPair pair = linkedPairs.get(lpIndex1);
                lowerSV = pair.first();
                lowerBreakend = pair.firstBreakend().getOtherBreakend();
            }
            else
            {
                LinkedPair prevPair = linkedPairs.get(lpIndex1 - 1);
                lowerSV = prevPair.second();
                lowerBreakend = prevPair.secondBreakend();
            }

            if(lowerSV.isSglBreakend() && lowerSV.getSglMappings().isEmpty())
                continue;

            // handle breakends from a SGL's mapping to known pair genes
            int lowerBreakendPos;
            List<BreakendGeneData> genesListLower;

            if(lowerBreakend != null)
            {
                lowerBreakendPos = lowerBreakend.position();
                genesListLower = getBreakendGeneList(lowerSV, lowerBreakend.usesStart());
            }
            else
            {
                genesListLower = getBreakendGeneList(lowerSV, false);

                if(genesListLower.isEmpty())
                    continue;

                lowerBreakendPos = genesListLower.get(0).position();
            }

            if(genesListLower.isEmpty())
                continue;

            final List<LinkedPair> traversedPairs = Lists.newArrayList();

            for(int lpIndex2 = lpIndex1; lpIndex2 <= linkedPairs.size(); ++lpIndex2)
            {
                SvVarData upperSV = null;
                SvBreakend upperBreakend = null;

                // the upper link takes the breakend of the current linked pair's 'first' SV
                // which for the first index will just be the curent SV's 2 breakends - ie testing a single SV in a chain
                // and beyond all the links, it must also test fusions with the chain's upper open breakend
                if(lpIndex2 < linkedPairs.size())
                {
                    LinkedPair pair = linkedPairs.get(lpIndex2);
                    upperSV = pair.first();
                    upperBreakend = pair.firstBreakend();
                }
                else if(chain.isClosedLoop())
                {
                    // take the breakend that loops around to the start
                    LinkedPair pair = linkedPairs.get(lpIndex2 - 1);
                    upperSV = pair.second();
                    upperBreakend = pair.secondBreakend().getOtherBreakend();
                }
                else
                {
                    // at the last link, take the open breakend of the chain
                    upperBreakend = chain.getOpenBreakend(false);
                    upperSV = chain.getChainEndSV(false);
                }

                List<BreakendGeneData> genesListUpper;

                if(upperBreakend != null)
                {
                    genesListUpper = getBreakendGeneList(upperSV, upperBreakend.usesStart());
                }
                else
                {
                    genesListUpper = getBreakendGeneList(upperSV, false);

                    if(genesListUpper.isEmpty())
                        continue;
                }

                // if a new linked pair has been traversed to reach this upper breakend, record its length
                // and test whether it traverses any genic region, and if so invalidate the fusion
                if(lpIndex2 > lpIndex1)
                {
                    LinkedPair lastPair = linkedPairs.get(lpIndex2 - 1);
                    traversedPairs.add(lastPair);
                }

                if(genesListUpper.isEmpty())
                {
                    // skip past this link and breakend to the next one, keeping the possibility of a fusion with the lower breakend open
                    continue;
                }

                // test the fusion between these 2 breakends
                List<GeneFusion> fusions = mFusionFinder.findFusions(genesListLower, genesListUpper);

                if(mNeoEpitopeWriter != null)
                {
                    final LinkedPair lowerLink = lpIndex1 == 0 ? null : linkedPairs.get(lpIndex1 - 1);
                    final LinkedPair upperLink = lpIndex2 < linkedPairs.size() ? linkedPairs.get(lpIndex2) : null;

                    double svCopyNumber = upperBreakend != null && lowerBreakend != null ?
                            (lowerBreakend.copyNumber() + upperBreakend.copyNumber()) * 0.5 : lowerSV.copyNumber(true);

                    mNeoEpitopeWriter.processFusionCandidate(
                            genesListLower, genesListUpper, traversedPairs, mDisruptionFinder, lowerLink, upperLink, svCopyNumber);
                }

                if(fusions.isEmpty())
                    continue;

                if(lpIndex2 > lpIndex1)
                {
                    // a chain cannot be an exon-exon fusion, so cull any of these
                    fusions = fusions.stream().filter(x -> !x.isExonic()).collect(Collectors.toList());
                }

                if(mFusionConfig.LogReportableOnly)
                {
                    fusions = fusions.stream().filter(x -> isReportable(x)).collect(Collectors.toList());

                    if(fusions.isEmpty())
                        continue;
                }

                int validTraversalFusionCount = 0; // between these 2 SVs

                for(GeneFusion fusion : fusions)
                {
                    // if the fusion from the upstream gene is on the positive strand, then it will have a fusion direction of +1
                    // whenever it goes through a subsequent linked pair by joining to the first (lower) breakend in the pair
                    int upGeneStrand = fusion.upstreamTrans().gene().strand();
                    boolean isPrecodingUpstream = fusion.upstreamTrans().preCoding();
                    boolean fusionLowerToUpper = fusion.upstreamTrans().gene().position() == lowerBreakendPos;

                    final List<BreakendGeneData> downstreamGenes = lowerBreakendPos == fusion.downstreamTrans().gene().position()
                            ? genesListLower : genesListUpper;

                    boolean nonDisruptiveChain = !traversedPairs.isEmpty() && !downstreamGenes.isEmpty() ?
                            mDisruptionFinder.isNonDisruptiveChainedPair(fusion.upstreamTrans(), downstreamGenes) : false;

                    // check any traversed genes
                    long totalLinkLength = 0;
                    boolean validTraversal = true;
                    boolean allowInvalidTraversal = allowSuspectChains(fusion.knownType());
                    boolean allTraversalAssembled = true;

                    for(LinkedPair pair : traversedPairs)
                    {
                        totalLinkLength += pair.baseLength();

                        if(pair.isInferred())
                            allTraversalAssembled = false;

                        // if going lower to upper, if the orientation of the first breakend in the pair is opposite to the strand of
                        // the upstream gene, then the fusion direction for that pair is the same as a the upstream gene
                        // otherwise it needs to be switched
                        int fusionDirection = 0;

                        if(fusionLowerToUpper)
                        {
                            fusionDirection = pair.firstBreakend().orientation() != upGeneStrand ? upGeneStrand : -upGeneStrand;
                        }
                        else
                        {
                            fusionDirection = pair.secondBreakend().orientation() != upGeneStrand ? upGeneStrand : -upGeneStrand;
                        }

                        // any invalid traversal causes this fusion to be entirely skipped from further analysis
                        if(validTraversal && !hasValidTraversal(validPairs, pair, fusionDirection, isPrecodingUpstream))
                        {
                            validTraversal = false;

                            if(!allowInvalidTraversal)
                                break;
                        }
                    }

                    if(!validTraversal && !allowInvalidTraversal)
                    {
                        recordInvalidFusion(fusion, "InvalidTraversal");
                        continue;
                    }

                    if(validTraversal)
                        ++validTraversalFusionCount;

                    boolean[] transTerminated = { false, false};

                    // check that the chain does not disrupt (ie terminate) the rest of the gene, checking up then downstream in turn
                    for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
                    {
                        BreakendTransData transcript = fs == FS_UP ? fusion.upstreamTrans() : fusion.downstreamTrans();
                        BreakendGeneData gene = fs == FS_UP ? fusion.upstreamTrans().gene() : fusion.downstreamTrans().gene();

                        boolean isLowerBreakend = lowerBreakendPos == gene.position();

                        boolean isChainEnd = (isLowerBreakend && lpIndex1 == 0) || (!isLowerBreakend && lpIndex2 == linkedPairs.size());

                        if(isChainEnd)
                        {
                            transTerminated[fs] = false;
                        }
                        else
                        {
                            SvBreakend breakend = isLowerBreakend ? lowerBreakend : upperBreakend;

                            if(breakend == null)
                            {
                                transTerminated[fs] = false;
                            }
                            else
                            {
                                // looks from this link outwards past the end of the transcript for any invalidation of the transcript
                                int linkIndex = isLowerBreakend ? lpIndex1 - 1 : lpIndex2;

                                if(linkIndex >= chain.getLinkedPairs().size())
                                {
                                    LNX_LOGGER.error("cluster({}) link index({}) exceeds chain link size({})",
                                            cluster.id(), linkIndex, chain.getLinkedPairs().size());
                                    break;
                                }

                                boolean logSkip = false; // fusion.knownType() != NONE;
                                transTerminated[fs] = checkTranscriptDisruptionInfo(breakend, transcript, chain, linkIndex, logSkip);
                            }
                        }
                    }

                    int linksCount = lpIndex2 - lpIndex1;

                    FusionChainInfo chainInfo = ImmutableFusionChainInfo.builder()
                            .chainId(chain.id())
                            .chainLinks(linksCount)
                            .chainLength((int)totalLinkLength)
                            .traversalAssembled(allTraversalAssembled)
                            .validTraversal(validTraversal)
                            .nonDisruptive(nonDisruptiveChain)
                            .build();

                    FusionAnnotations annotations = ImmutableFusionAnnotations.builder()
                            .clusterId(cluster.id())
                            .clusterCount(cluster.getSvCount())
                            .resolvedType(cluster.getResolvedType().toString())
                            .chainInfo(chainInfo)
                            .terminatedUp(transTerminated[SE_START])
                            .terminatedDown(transTerminated[SE_END])
                            .build();

                    fusion.setAnnotations(annotations);

                    // accept invalidated chains and transcripts for known fusions
                    boolean chainLengthOk = totalLinkLength <= FUSION_MAX_CHAIN_LENGTH;
                    boolean traversalOk = validTraversal || allowInvalidTraversal;

                    if(chainLengthOk && traversalOk && !nonDisruptiveChain)
                    {
                        if(!hasIdenticalFusion(fusion, chainFusions))
                        {
                            chainFusions.add(fusion);
                        }
                    }
                    else
                    {
                        String invalidReason = "Unknown";

                        if(!validTraversal)
                            invalidReason = "TraversesSPA";
                        else if(!chainLengthOk)
                            invalidReason = "LongChain";
                        else if(nonDisruptiveChain)
                            invalidReason = "NonDisruptiveChain";

                        recordInvalidFusion(fusion, invalidReason);
                    }
                }

                if(validTraversalFusionCount == 0 && lpIndex2 > lpIndex1 && mRnaFusionMapper == null)
                {
                    // there were fusions between these 2 SV breakends but none had valid traversal
                    // this means than any any link further along the chain will also not be valid
                    LNX_LOGGER.trace("cluster({}) chain({}) no valid traversals between({} -> {})",
                            cluster.id(), chain.id(), lpIndex1, lpIndex2);

                    break;
                }
            }
        }
    }

    private boolean hasValidTraversal(
            final List<ValidTraversalData> validPairs, final LinkedPair pair, int fusionDirection, boolean isPrecodingUpstream)
    {
        final ValidTraversalData existingData = validPairs.stream()
                .filter(x -> x.matches(pair, fusionDirection, isPrecodingUpstream)).findFirst().orElse(null);

        if(existingData != null)
            return existingData.IsValid;

        boolean validTraversal = !mDisruptionFinder.pairTraversesGene(pair, fusionDirection, isPrecodingUpstream);
        validPairs.add(new ValidTraversalData(pair, validTraversal, fusionDirection, isPrecodingUpstream));
        return validTraversal;
    }

    private void recordInvalidFusion(final GeneFusion fusion, final String reason)
    {
        if(mInvalidFusions.keySet().stream().anyMatch(x -> fusion.name().equals(x.name())))
            return;

        mInvalidFusions.put(fusion, reason);
    }

    private boolean hasIdenticalFusion(final GeneFusion newFusion, final List<GeneFusion> fusions)
    {
        for(GeneFusion fusion : fusions)
        {
            if(newFusion.upstreamTrans().gene().id() != fusion.upstreamTrans().gene().id())
                continue;

            if(newFusion.downstreamTrans().gene().id() != fusion.downstreamTrans().gene().id())
                continue;

            if(!newFusion.upstreamTrans().transName().equals(fusion.upstreamTrans().transName()))
                continue;

            if(!newFusion.downstreamTrans().transName().equals(fusion.downstreamTrans().transName()))
                continue;

            if(newFusion.downstreamTrans().hasNegativePrevSpliceAcceptorDistance() != fusion.downstreamTrans().hasNegativePrevSpliceAcceptorDistance())
                continue;

            return true;
        }

        return false;
    }

    private List<BreakendGeneData> getBreakendGeneList(final SvVarData var, boolean isStart)
    {
        if(var.isSglBreakend() && !isStart)
        {
            // limit to known fusion genes
            final List<BreakendGeneData> genesList = var.getGenesList(false);

            // take the genes if any are in an IG region
            if(var.getSglMappings().stream().anyMatch(x -> mFusionFinder.getKnownFusionCache().withinIgRegion(x.Chromosome, x.Position)))
                return genesList;

            // otherwise check known pairs (including IG pair 3' genes)
            return genesList.stream()
                    .filter(x -> mFusionFinder.getKnownFusionCache().isSingleBreakendCandidate(x.geneName(), x.isUpstream()))
                    .collect(Collectors.toList());

            // mFusionFinder.getKnownFusionCache().withinIgRegion(var.chromosome(false), var.position(false))
        }

        return var.getGenesList(isStart);
    }

    private boolean checkTranscriptDisruptionInfo(
            final SvBreakend breakend, final BreakendTransData transcript, final SvChain chain, int linkIndex, boolean logSkip)
    {
        // return true if the transcript is disrupted before the chain leaves it

        // starting with this breakend and working onwards from it in the chain, check for any disruptions to the transcript
        LinkedPair startPair = chain.getLinkedPairs().get(linkIndex);
        boolean traverseUp = startPair.firstBreakend() == breakend; // whether to search up or down the chain

        boolean transcriptTerminated = false;

        boolean isUpstream = transcript.isUpstream();

        while(linkIndex >= 0 && linkIndex <= chain.getLinkedPairs().size() - 1)
        {
            LinkedPair pair = chain.getLinkedPairs().get(linkIndex);

            SvBreakend nextBreakend = traverseUp ? pair.secondBreakend() : pair.firstBreakend();

            // the chain does not terminate the transcript if the next breakend is now past the end of the transcript or
            // if the breakend is down-stream of coding
            if(nextBreakend.orientation() == 1)
            {
                if(nextBreakend.position() > transcript.transEnd())
                    break;
                else if(!isUpstream && transcript.codingEnd() > 0 && nextBreakend.position() > transcript.codingEnd())
                    break;
            }
            else
            {
                if(nextBreakend.position() < transcript.transStart())
                    break;
                else if(!isUpstream && transcript.codingStart() > 0 && nextBreakend.position() < transcript.codingStart())
                    break;
            }

            // continue on if the next SV (ie both its breakends) is non-disruptive in this transcript
            if(DisruptionFinder.isDisruptiveInTranscript(nextBreakend.getSV(), transcript.TransData))
            {
                transcriptTerminated = true;
                break;
            }

            if(logSkip)
            {
                LNX_LOGGER.info("sample({}) cluster({}) chain({}) skipping non-disruptive SV({})",
                        mSampleId, nextBreakend.getSV().getCluster().id(), chain.id(), nextBreakend.getSV().toString());
            }

            linkIndex += traverseUp ? 1 : -1;
        }

        return transcriptTerminated;
    }

    private boolean persistFusion(final GeneFusion fusion)
    {
        if(!fusion.validChainTraversal() && !allowSuspectChains(fusion.knownType()))
            return false;

        if(!fusion.isIG() && fusion.downstreamTrans().hasNegativePrevSpliceAcceptorDistance())
            return false;

        return true;
    }

    private List<GeneFusion> extractUniqueFusions()
    {
        final Map<String,List<GeneFusion>> svIdPairFusions = Maps.newHashMap();
        final Map<String,List<GeneFusion>> genePairFusions = Maps.newHashMap();

        // allocate fusions to unique SV pairings
        for(final GeneFusion fusion : mFusions)
        {
            if(!persistFusion(fusion))
                continue;

            final String svIdPair = fusion.svIdPair();

            List<GeneFusion> fusions = svIdPairFusions.get(svIdPair);

            if(fusions == null)
                svIdPairFusions.put(svIdPair, Lists.newArrayList(fusion));
            else
                fusions.add(fusion);
        }

        // find top priority fusion amongst these fusions by SV-pair, setting reportability at the same time
        for(Map.Entry<String,List<GeneFusion>> entry : svIdPairFusions.entrySet())
        {
            List<GeneFusion> similarFusions = entry.getValue();

            // IG promiscuous fusions are not reportable but need to all be recorded
            List<GeneFusion> candidates = similarFusions.stream()
                    .filter(x -> x.knownType() == IG_PROMISCUOUS)
                    .filter(x -> x.downstreamTrans().TransData.IsCanonical)
                    .collect(Collectors.toList());

            final GeneFusion topFusion = mFusionFinder.findTopReportableFusion(similarFusions);

            if(topFusion != null && !candidates.contains(topFusion))
            {
                candidates.add(topFusion);
            }

            for(GeneFusion fusion : candidates)
            {
                List<GeneFusion> fusionsByName = genePairFusions.get(fusion.name());

                if(fusionsByName == null)
                    genePairFusions.put(fusion.name(), Lists.newArrayList(fusion));
                else
                    fusionsByName.add(fusion);
            }
        }

        List<GeneFusion> uniqueFusions = Lists.newArrayList();

        // if any gene-pair has multiple candidates, once again take the top one and
        for(Map.Entry<String,List<GeneFusion>> entry : genePairFusions.entrySet())
        {
            final List<GeneFusion> fusions = entry.getValue();
            if(fusions.size() == 1)
            {
                uniqueFusions.add(fusions.get(0));
            }
            else
            {
                GeneFusion topFusion = findTopPriorityFusion(fusions);
                if(topFusion != null)
                    uniqueFusions.add(topFusion);
            }
        }

        return uniqueFusions;
    }

    private List<BreakendTransData> buildTranscriptList(final List<SvVarData> svList, final List<GeneFusion> fusions)
    {
        // add all canonical or otherwise reportable transcript and then add any additional transcripts from the fusions
        List<BreakendTransData> transcripts = Lists.newArrayList();

        for(SvVarData var : svList)
        {
            for(int be = SE_START; be <= SE_END; ++be)
            {
                // by default don't include SGL alt-mappings - they'll be added if in a fusion
                if(var.isSglBreakend() && be == SE_END)
                    continue;

                for(BreakendGeneData geneAnnotation : var.getGenesList(isStart(be)))
                {
                    for(BreakendTransData transcript : geneAnnotation.transcripts())
                    {
                        if(transcript.isCanonical())
                            transcripts.add(transcript);
                        else if(mDisruptionFinder.matchesDisruptionTranscript(geneAnnotation.geneId(), transcript.TransData))
                            transcripts.add(transcript);
                    }
                }
            }
        }

        for(GeneFusion fusion : fusions)
        {
            if(!transcripts.contains(fusion.upstreamTrans()))
                transcripts.add(fusion.upstreamTrans());

            if(!transcripts.contains(fusion.downstreamTrans()))
                transcripts.add(fusion.downstreamTrans());
        }

        // set undisrupted copy number for all persisted breakends
        for(BreakendTransData transcript : transcripts)
        {
            if(transcript.undisruptedCopyNumberSet())
                continue;

            SvVarData var = svList.stream().filter(x -> x.id() == transcript.gene().id()).findFirst().orElse(null);

            if(var != null)
            {
                final SvBreakend breakend = var.getBreakend(transcript.gene().isStart());

                if(breakend != null)
                    transcript.setUndisruptedCopyNumber(getUndisruptedCopyNumber(breakend));
            }
        }

        return transcripts;
    }

    private void addVisualisationData(final List<GeneFusion> fusionList)
    {
        if(mVisSampleData == null || !mConfig.Output.WriteVisualisationData)
            return;

        final List<VisFusion> visFusions = Lists.newArrayList();

        for(final GeneFusion fusion : fusionList)
        {
            if(fusion.reportable() || mFusionConfig.WriteAllVisFusions)
            {
                int clusterId = fusion.getAnnotations() != null ? fusion.getAnnotations().clusterId() : -1;

                final BreakendTransData transUp = fusion.upstreamTrans();
                final BreakendTransData transDown = fusion.downstreamTrans();

                mVisSampleData.addGeneExonData(clusterId, transUp.gene().geneId(), transUp.geneName(),
                        transUp.transName(), transUp.transId(), transUp.gene().chromosome(), FUSION);

                mVisSampleData.addGeneExonData(clusterId, transDown.gene().geneId(), transDown.geneName(),
                        transDown.transName(), transDown.transId(), transDown.gene().chromosome(), FUSION);

                visFusions.add(new VisFusion(
                        mSampleId, clusterId, fusion.reportable(),
                        transUp.geneName(), transUp.transName(), transUp.gene().chromosome(), transUp.gene().position(),
                        transUp.gene().strand(), transUp.regionType().toString(), fusion.getFusedExon(true),
                        transDown.geneName(), transDown.transName(), transDown.gene().chromosome(), transDown.gene().position(),
                        transDown.gene().strand(), transDown.regionType().toString(), fusion.getFusedExon(false)));
            }
        }

        mVisSampleData.addFusions(visFusions);
    }

    public void close()
    {
        if(mRnaFusionMapper != null)
            mRnaFusionMapper.close();

        if(mNeoEpitopeWriter != null)
            mNeoEpitopeWriter.close();
    }
}
