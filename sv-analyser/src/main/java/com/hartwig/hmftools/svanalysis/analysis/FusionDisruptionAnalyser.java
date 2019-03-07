package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_LENGTH;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_LINK_COUNT;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.EXON_RANK_MIN;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.svannotation.analysis.RnaFusionData.RNA_SPLICE_TYPE_ONLY_REF;
import static com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser.checkFusionLogic;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.svanalysis.annotators.VisualiserWriter;
import com.hartwig.hmftools.svanalysis.types.SvFusion;
import com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svannotation.analysis.RnaFusionData;
import com.hartwig.hmftools.svannotation.analysis.SvDisruptionAnalyser;
import com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class FusionDisruptionAnalyser
{
    private SvFusionAnalyser mFusionFinder;
    private SvDisruptionAnalyser mDisruptionFinder;

    private String mSampleId;
    private String mOutputDir;
    private SvGeneTranscriptCollection mEnsemblDataCache;
    private Map<String, List<SvBreakend>> mChrBreakendMap;

    private boolean mSkipFusionOutput;
    private List<GeneFusion> mFusions;
    private List<SvFusion> mSvFusionList;

    ListMultimap<Chromosome, HmfTranscriptRegion> mChromosomeTranscriptMap;

    private VisualiserWriter mVisWriter;

    private static final Logger LOGGER = LogManager.getLogger(FusionDisruptionAnalyser.class);

    public FusionDisruptionAnalyser()
    {
        mFusionFinder = null;
        mDisruptionFinder = null;
        mEnsemblDataCache = new SvGeneTranscriptCollection();
        mChromosomeTranscriptMap = null;
        mOutputDir = "";
        mFusions = Lists.newArrayList();
        mSkipFusionOutput = false;
        mSvFusionList = Lists.newArrayList();
        mVisWriter = null;
        mChrBreakendMap = null;
    }

    public void skipFusionOutput(boolean toggle) { mSkipFusionOutput = toggle; }

    public void loadFusionReferenceData(final CommandLine cmdLineArgs, final String outputDir, SvGeneTranscriptCollection ensemblDataCache)
    {
        mOutputDir = outputDir;

        mEnsemblDataCache = ensemblDataCache;
        mFusionFinder = new SvFusionAnalyser(cmdLineArgs, ensemblDataCache, mOutputDir);

        List<HmfTranscriptRegion> transcriptRegions = HmfGenePanelSupplier.allGeneList37();
        mChromosomeTranscriptMap = Multimaps.fromRegions(transcriptRegions);
    }

    public final Set<String> getRnaSampleIds() { return mFusionFinder.getSampleRnaData().keySet(); }
    public final List<GeneFusion> getFusions() { return mFusions; }
    public void setVisWriter(VisualiserWriter writer) { mVisWriter = writer; }

    private void setSvGenesList(final SvVarData var, boolean applyPromotorDistance)
    {
        List<GeneAnnotation> genesList = Lists.newArrayList();

        int upstreamDistance = applyPromotorDistance ? PRE_GENE_PROMOTOR_DISTANCE : 0;

        for(int be = SVI_START; be <= SVI_END; ++be)
        {
            if(be == SVI_END && var.isNullBreakend())
                continue;

            boolean isStart = isStart(be);

            genesList.addAll(mEnsemblDataCache.findGeneAnnotationsBySv(
                    var.dbId(), isStart, var.chromosome(isStart), var.position(isStart), upstreamDistance));
        }

        if(genesList.isEmpty())
            return;

        List<GeneAnnotation> startGenes = Lists.newArrayList();
        List<GeneAnnotation> endGenes = Lists.newArrayList();

        for(GeneAnnotation gene : genesList)
        {
            gene.setSvData(var.getSvData());

            if(gene.isStart())
                startGenes.add(gene);
            else
                endGenes.add(gene);
        }

        var.setGenesList(startGenes, true);
        var.setGenesList(endGenes, false);
    }

    public void setSvGeneData(final String sampleId, final List<SvVarData> svList, boolean applyPromotorDistance,
            Map<String, List<SvBreakend>> chrBreakendMap)
    {
        mSampleId = sampleId;
        mChrBreakendMap = chrBreakendMap;

        for(final SvVarData var : svList)
        {
            if (var.isReplicatedSv())
                continue;

            // cache transcript info against each SV
            setSvGenesList(var, applyPromotorDistance);
        }
    }

    public void run(final List<SvVarData> svList, final List<SvCluster> clusters)
    {
        findFusions(svList, clusters);

        assessRnaFusions();
    }

    private void findFusions(final List<SvVarData> svList, final List<SvCluster> clusters)
    {
        if(mSampleId.isEmpty() || mFusionFinder == null)
            return;

        mFusions.clear();

        // always report SVs by themselves
        for(final SvVarData var : svList)
        {
            if(var.isReplicatedSv())
                continue;

            if(var.isNullBreakend())
                continue;

            checkFusions(var.getGenesList(true), var.getGenesList(false), var.getCluster());
        }

        boolean checkClusters = false;

        if(checkClusters)
        {
            int maxClusterSize = 50;

            // for now only consider simple SVs and resolved small clusters
            for (final SvCluster cluster : clusters)
            {
                if (cluster.getSvCount() == 1) // simple clusters already checked
                    continue;

                // if(cluster.hasReplicatedSVs() || !cluster.isFullyChained() || cluster.getTypeCount(SGL) > 0)
                //    continue;

                if (cluster.getSvCount() > maxClusterSize)
                    continue;

                for (final SvChain chain : cluster.getChains())
                {
                    findChainedFusion(cluster, chain);
                }
            }
        }
    }

    private void findChainedFusion(final SvCluster cluster, final SvChain chain)
    {
        final List<SvLinkedPair> linkedPairs = chain.getLinkedPairs();

        // the aim here is to find any 2 fused genes from any 2 links in the chain
        // if say SV index 1 is fused with SV index 4, going through linked pairs 1-2, 2-3 and 3-4,
        // the fusion will come from the other breakend of the 1-2 link (being the section fused in)
        // with the other breakend of the 3-4 link

        // so given that each linked pair's 'first' SV points up and 'second' points down the chain,
        // a fusion will come from the lower link's 'second' breakend region and the upper link's first' breakend section
        boolean inPossibleFusion = false;
        List<SvLinkedPair> traversedLinks = Lists.newArrayList();
        int startIndex = -1;
        int endIndex = -1;
        SvBreakend beStart = null;
        SvBreakend beEnd = null;
        List<GeneAnnotation> genesListStart = null;

        for(int lpIndex = 0; lpIndex < linkedPairs.size(); ++lpIndex)
        {
            final SvLinkedPair linkedPair = linkedPairs.get(lpIndex);

            final SvVarData varStart = linkedPair.first();
            final SvVarData varEnd = linkedPair.second();

            if (!inPossibleFusion)
            {
                // check whether the first breakend's other end falls within a gene
                genesListStart = varStart.getGenesList(linkedPair.firstUnlinkedOnStart());

                if (genesListStart.isEmpty())
                    continue;

                startIndex = lpIndex;
                beStart = varStart.getBreakend(linkedPair.firstUnlinkedOnStart());
                inPossibleFusion = true;
            }

            // test the end link
            List<GeneAnnotation> genesListEnd = varEnd.getGenesList(linkedPair.secondUnlinkedOnStart());

            if (genesListEnd.isEmpty())
            {
                // add this link to be checked for genes if a fusions is found
                traversedLinks.add(linkedPair);
                continue;
            }

            endIndex = lpIndex;
            beEnd = varEnd.getBreakend(linkedPair.secondUnlinkedOnStart());

            // check any traversed section before continuing on
            boolean traversesGene = false;
            long traversalLength = 0;

            for(final SvLinkedPair pair : traversedLinks)
            {
                final String chr = pair.first().chromosome(pair.firstLinkOnStart());
                long firstPos = pair.first().position(pair.firstLinkOnStart());
                long secondPos = pair.second().position(pair.secondLinkOnStart());
                long startPos = firstPos < secondPos ? firstPos : secondPos;
                long endPos = firstPos > secondPos ? firstPos : secondPos;

                traversalLength += endPos - startPos;

                List<HmfTranscriptRegion> transcripts = findGenesForRegion(chr, startPos, endPos);

                if(transcripts.isEmpty())
                    continue;

                traversesGene = true;

                LOGGER.info("cluster({}) chain({}) potential fusion: be1({} {}) & be2({} {}) indices({} -> {}) section({}: {} - {}) traverses {} genes (eg {})",
                        cluster.id(), chain.id(), beStart.toString(), genesListStart.get(0).GeneName, beEnd.toString(), genesListEnd.get(0).GeneName,
                        startIndex, endIndex, chr, startPos, endPos, transcripts.size(), transcripts.get(0).geneID());

                break;
            }

            if(!traversesGene)
            {
                if (startIndex < endIndex)
                {
                    LOGGER.info("cluster({}) chain({}) potential fusion: be1({} {}) & be2({} {}) link indices({} -> {}) traversalLen({})",
                            cluster.id(), chain.id(), beStart.toString(), genesListStart.get(0).GeneName,
                            beEnd.toString(), genesListEnd.get(0).GeneName, startIndex, endIndex, traversalLength);
                }

                checkFusions(genesListStart, genesListEnd, cluster);
            }

            // reset state
            inPossibleFusion = false;
            traversedLinks.clear();
            startIndex = -1;
            beStart = null;
            genesListStart = null;
        }
    }

    private List<HmfTranscriptRegion> findGenesForRegion(final String chromosome, long startPos, long endPos)
    {
        List<HmfTranscriptRegion> coveredTranscripts = Lists.newArrayList();
        final List<HmfTranscriptRegion> allTranscripts = mChromosomeTranscriptMap.get(HumanChromosome.fromString(chromosome));

        if(allTranscripts == null)
            return coveredTranscripts;

        for(final HmfTranscriptRegion transRegion : allTranscripts)
        {
            if(startPos <= transRegion.geneStart() && endPos >= transRegion.geneEnd())
            {
                coveredTranscripts.add(transRegion);
            }
        }

        return coveredTranscripts;
    }

    private void checkFusions(List<GeneAnnotation> breakendGenes1, List<GeneAnnotation> breakendGenes2, final SvCluster cluster)
    {
        if (breakendGenes1.isEmpty() || breakendGenes2.isEmpty())
            return;

        List<GeneFusion> fusions = mFusionFinder.findFusions(breakendGenes1, breakendGenes2);

        if (fusions.isEmpty())
            return;

        // mFusions.addAll(fusions);

        if(LOGGER.isDebugEnabled())
        {
            for (final GeneFusion fusion : fusions)
            {
                if (fusion.reportable())
                {
                    final Transcript upstream = fusion.upstreamTrans();
                    final Transcript downstream = fusion.downstreamTrans();
                    final GeneAnnotation upGene = upstream.parent();
                    final GeneAnnotation downGene = downstream.parent();

                    LOGGER.debug("sample({}) fusion: up({} {} {} {} ph={}) upSV({}: {}:{}:{} start={} strand={}) down({} {} {} {} ph={}) downSV({}: {}:{}:{} start={} strand={})",
                            mSampleId, upstream.geneName(), upstream.StableId, upstream.regionType(), upstream.codingType(), upstream.exonUpstreamPhase(),
                            upGene.id(), upGene.chromosome(), upGene.position(), upGene.orientation(), upGene.isStart(), upGene.Strand,
                            downstream.geneName(), downstream.StableId, downstream.regionType(), downstream.codingType(), downstream.exonDownstreamPhase(),
                            downGene.id(), downGene.chromosome(), downGene.position(), downGene.orientation(), downGene.isStart(), downGene.Strand);
                }
            }
        }

        String clusterInfo = String.format("%d,%d,%s", cluster.id(), cluster.getSvCount(), cluster.getResolvedType());

        if(!mSkipFusionOutput)
        {
            mFusionFinder.writeFusions(fusions, mSampleId, clusterInfo, true);
        }

        if(mVisWriter != null)
        {
            for (final GeneFusion fusion : fusions)
            {
                if(fusion.reportable())
                {
                    mVisWriter.addGeneExonData(cluster.id(),
                            fusion.upstreamTrans().parent().StableId, fusion.upstreamTrans()
                                    .parent().GeneName, fusion.upstreamTrans().StableId,
                            fusion.upstreamTrans().parent().chromosome(), "FUSION");

                    mVisWriter.addGeneExonData(cluster.id(),
                            fusion.downstreamTrans().parent().StableId, fusion.downstreamTrans()
                                    .parent().GeneName, fusion.downstreamTrans().StableId,
                            fusion.downstreamTrans().parent().chromosome(), "FUSION");
                }
            }
        }
    }

    private void assessRnaFusions()
    {
        final List<RnaFusionData> rnaFusionList = mFusionFinder.getSampleRnaData().get(mSampleId);

        if (rnaFusionList == null || rnaFusionList.isEmpty())
            return;

        for (final RnaFusionData rnaFusion : rnaFusionList)
        {
            mFusionFinder.setRnaFusionData(rnaFusion);

            annotateRnaFusions(rnaFusion);

            mFusionFinder.writeRnaMatchData(mSampleId, rnaFusion);
        }
    }

    public void annotateRnaFusions(final RnaFusionData rnaFusion)
    {

        /* Matching and annotation logic:
            - find all breakends in the RNA up and down gene
            - for them, find the any transcripts which a) have the exon boundary in the RNA position AND
            - b) are in the correct relative position:
                - upstream: at or after the RNA boundary down to the start of the next exon
                - downstream: at or before the RNA bounday up to the start of the preceding exon
                - If a transcript on the downstream gene starts on the 2nd exon, the fusion is allowed to match up to the nearer
                of a splice acceptor site with same orientation on a previous gene OR 100k bases upstream of the transcript.
                (is the distance up available for these??)
                - if multiple transcripts exist for the same breakend, take the canonical or longest (but it should make no difference)
            - if multiple breakends meet these criteria at either end, prioritise in the following order
                - both breakends are either end of the same structural variant
                - both breakends are in the same chain
                - both breakends are in the same cluster
                - otherwise take the nearest breakend to the RNA position
            - if no breakend is found on either upstream or downstream gene meeting the above criteria then record the nearest ID,
            distance and min number of skipped splice sites.
        */


        // viable breakends and their matching transcript
        List<SvBreakend> upstreamBreakends = Lists.newArrayList();
        List<SvBreakend> downstreamBreakends = Lists.newArrayList();
        List<Transcript> upstreamTranscripts = Lists.newArrayList();
        List<Transcript> downstreamTranscripts = Lists.newArrayList();

        // non-viable transcripts to be used if non viable ones are found
        List<Transcript> nearUpstreamTranscripts = Lists.newArrayList();
        List<Transcript> nearDownstreamTranscripts = Lists.newArrayList();
        List<SvBreakend> nearUpstreamBreakends = Lists.newArrayList();
        List<SvBreakend> nearDownstreamBreakends = Lists.newArrayList();

        boolean isExactRnaExon = rnaFusion.SpliceType.equals(RNA_SPLICE_TYPE_ONLY_REF);

        for(int i = 0; i <= 1 ; ++i)
        {
            boolean isUpstream = (i == 0);
            String chromosome = isUpstream ? rnaFusion.ChrUp : rnaFusion.ChrDown;
            long rnaPosition = isUpstream ? rnaFusion.PositionUp : rnaFusion.PositionDown;
            byte geneStrand = isUpstream ? rnaFusion.StrandUp : rnaFusion.StrandDown;
            List<SvBreakend> viableBreakends = isUpstream ? upstreamBreakends : downstreamBreakends;
            List<SvBreakend> nearBreakends = isUpstream ? nearUpstreamBreakends : nearDownstreamBreakends;
            List<Transcript> viableTranscripts = isUpstream ? upstreamTranscripts : downstreamTranscripts;
            List<Transcript> nearTranscripts = isUpstream ? nearUpstreamTranscripts : nearDownstreamTranscripts;
            String geneName = isUpstream ? rnaFusion.GeneUp : rnaFusion.GeneDown;

            final List<SvBreakend> breakendList = mChrBreakendMap.get(chromosome);

            if(breakendList == null)
                continue;

            for(final SvBreakend breakend : breakendList)
            {
                final SvVarData var = breakend.getSV();

                if(var.isNoneSegment())
                    continue;

                // check breakend falls in genic region
                List<GeneAnnotation> genesList = var.getGenesList(breakend.usesStart())
                        .stream()
                        .filter(x -> x.GeneName.equals(geneName))
                        .collect(Collectors.toList());

                if(genesList.isEmpty())
                    continue;

                // check that breakend has correct orientation and position relative to RNA breakend
                if(!isViableBreakend(breakend, rnaPosition, geneStrand, isUpstream))
                    continue;

                // check whether any of the breakend's transcripts falls within the nearest exon of the RNA fusion breakpoint
                for(final Transcript trans : genesList.get(0).transcripts())
                {
                    if(trans.isCanonical())
                    {
                        nearBreakends.add(breakend);
                        nearTranscripts.add(trans);
                    }

                    if(mFusionFinder.isTranscriptBreakendViableForRnaBoundary(
                            trans, isUpstream,  breakend.position(), rnaPosition, isExactRnaExon))
                    {
                        viableBreakends.add(breakend);
                        viableTranscripts.add(trans);
                        break;
                    }
                }
            }
        }

        // run them through fusion logic (ie a pair of breakend lists), but don't require phase matching
        if(!upstreamBreakends.isEmpty() && !downstreamBreakends.isEmpty())
        {
            GeneFusion topCandidateFusion = null;
            SvBreakend topUpBreakend = null;
            SvBreakend topDownBreakend = null;

            for (int i = 0; i < upstreamBreakends.size(); ++i)
            {
                final SvBreakend upBreakend = upstreamBreakends.get(i);
                final Transcript upTrans = upstreamTranscripts.get(i);

                if(upBreakend.getSV().isNullBreakend())
                    continue;

                for (int j = 0; j < downstreamBreakends.size(); ++j)
                {
                    final SvBreakend downBreakend = downstreamBreakends.get(j);
                    final Transcript downTrans = downstreamTranscripts.get(j);

                    if(downBreakend.getSV().isNullBreakend())
                        continue;

                    GeneFusion possibleFusion = checkFusionLogic(upTrans, downTrans, false);

                    // form one any way but mark it as not meeting standard fusion rules
                    if(possibleFusion == null)
                        possibleFusion = new GeneFusion(upTrans, downTrans, false, false);

                    if (topCandidateFusion == null
                    || isCandidateBetter(topCandidateFusion, topUpBreakend, topDownBreakend, possibleFusion, upBreakend, downBreakend, rnaFusion))
                    {
                        topCandidateFusion = possibleFusion;
                        topUpBreakend = upBreakend;
                        topDownBreakend = downBreakend;
                    }
                }
            }

            if(topCandidateFusion != null)
            {
                rnaFusion.setTranscriptData(topCandidateFusion.upstreamTrans(), true, true, 0);
                rnaFusion.setTranscriptData(topCandidateFusion.downstreamTrans(), false, true, 0);
                rnaFusion.setViableFusion(topCandidateFusion.viable() && topCandidateFusion.phaseMatched());

                // add cluster and chain info
                setFusionClusterChainInfo(rnaFusion, topUpBreakend, topDownBreakend);
            }
        }
        else
        {
            // select the closest breakend's transcript
            for(int i = 0; i <= 1 ; ++i)
            {
                boolean isUpstream = (i == 0);
                long rnaPosition = isUpstream ? rnaFusion.PositionUp : rnaFusion.PositionDown;

                List<Transcript> transcriptList;
                List<SvBreakend> breakendList;
                boolean isViable = false;

                // use the viable transcripts if present, otherwise the nearest
                if(isUpstream)
                {
                    if(!upstreamTranscripts.isEmpty())
                    {
                        isViable = true;
                        transcriptList = upstreamTranscripts;
                        breakendList = upstreamBreakends;
                    }
                    else
                    {
                        transcriptList = nearUpstreamTranscripts;
                        breakendList = nearUpstreamBreakends;
                    }
                }
                else
                {
                    if(!downstreamTranscripts.isEmpty())
                    {
                        isViable = true;
                        transcriptList = downstreamTranscripts;
                        breakendList = downstreamBreakends;
                    }
                    else
                    {
                        transcriptList = nearDownstreamTranscripts;
                        breakendList = nearDownstreamBreakends;
                    }
                }

                Transcript closestTrans = null;
                SvBreakend closestBreakend = null;
                long closestDistance = 0;

                for (int j = 0; j < transcriptList.size(); ++j)
                {
                    final Transcript trans = transcriptList.get(j);
                    final SvBreakend breakend = breakendList.get(j);

                    long distance = abs(rnaPosition - trans.svPosition());
                    if(closestTrans == null || distance < closestDistance)
                    {
                        closestDistance = distance;
                        closestTrans = trans;
                        closestBreakend = breakend;
                    }
                }

                if(closestTrans != null)
                {
                    if(isViable)
                    {
                        rnaFusion.setTranscriptData(closestTrans, isUpstream, isViable, 0);
                    }
                    else
                    {
                        // for non-viable breakends, provide the exons skipped count
                        final String geneId = isUpstream ? rnaFusion.GeneUp : rnaFusion.GeneDown;
                        final int rnaExonData[] = mEnsemblDataCache.getExonRankings(geneId, rnaPosition);
                        final int svPosExonData[] = mEnsemblDataCache.getExonRankings(geneId, closestBreakend.position());

                        int exonsSkipped = abs(rnaExonData[EXON_RANK_MIN] - svPosExonData[EXON_RANK_MIN]);
                        rnaFusion.setTranscriptData(closestTrans, isUpstream, isViable, exonsSkipped);
                    }

                    setFusionClusterInfo(rnaFusion, closestBreakend, isUpstream);
                }
            }
        }
    }

    private void setFusionClusterChainInfo(final RnaFusionData rnaFusionData, final SvBreakend breakendUp, final SvBreakend breakendDown)
    {
        setFusionClusterInfo(rnaFusionData, breakendUp, true);
        setFusionClusterInfo(rnaFusionData, breakendDown, false);

        if(breakendUp.getSV().equals(breakendDown.getSV(), true))
            return;

        SvCluster clusterUp = breakendUp.getSV().getCluster();
        SvChain chainUp = clusterUp.findChain(breakendUp.getSV());

        SvCluster clusterDown = breakendDown.getSV().getCluster();
        SvChain chainDown = clusterDown.findChain(breakendDown.getSV());

        if(chainUp != null && chainUp == chainDown)
        {
            final int chainData[] = chainUp.breakendsAreChained(breakendUp.getSV(), breakendUp.usesStart(), breakendDown.getSV(), breakendUp.usesStart());
            final String chainInfo = String.format("%d;%d", chainData[CHAIN_LINK_COUNT], chainData[CHAIN_LENGTH]);

            // data: ChainLinks;ChainLength
            rnaFusionData.setChainInfo(chainInfo);
        }
    }

    private void setFusionClusterInfo(final RnaFusionData rnaFusionData, final SvBreakend breakend, boolean isUpstream)
    {
        // data: ClusterId;ClusterCount;ChainId;ChainCount
        SvCluster cluster = breakend.getSV().getCluster();
        SvChain chain = cluster.findChain(breakend.getSV());

        final String clusterInfo = String.format("%d;%d;%d;%d",
                cluster.id(), cluster.getSvCount(),
                chain != null ? chain.id() : cluster.getChainId(breakend.getSV()), chain != null ? chain.getSvCount() : 1);

        rnaFusionData.setClusterInfo(clusterInfo, isUpstream);
    }

    private boolean isCandidateBetter(final GeneFusion currentFusion, final SvBreakend beCurrentStart, final SvBreakend beCurrentEnd,
            final GeneFusion candidateFusion, final SvBreakend beCandidateStart, final SvBreakend beCandidateEnd, final RnaFusionData rnaFusion)
    {
        /*
            if(fusion.reportable() && fusion.phaseMatched())
            {
                return fusion;
            }
            else if(fusion.phaseMatched())
            {
                if(possibleFusion == null || !possibleFusion.phaseMatched())
                    possibleFusion = fusion;
            }
            else if(fusion.reportable())
            {
                if(possibleFusion == null || !possibleFusion.reportable())
                    possibleFusion = fusion;
            }
            else if(possibleFusion == null)
            {
                possibleFusion = fusion;
            }
         */

        // give priority to same SV
        boolean currentSameSV = beCurrentStart.getSV() == beCurrentEnd.getSV();
        boolean candidateSameSV = beCandidateStart.getSV() == beCandidateEnd.getSV();

        if(currentSameSV != candidateSameSV)
            return candidateSameSV;

        // then whether chained
        final SvChain currentChainStart = beCurrentStart.getSV().getCluster().findChain(beCurrentStart.getSV());
        boolean currentSameChain = currentChainStart != null && currentChainStart == beCurrentEnd.getSV().getCluster().findChain(beCurrentEnd.getSV());

        final SvChain candidateChainStart = beCandidateStart.getSV().getCluster().findChain(beCandidateStart.getSV());
        boolean candidateSameChain = candidateChainStart != null && candidateChainStart == beCandidateEnd.getSV().getCluster().findChain(beCandidateEnd.getSV());

        if(currentSameChain != candidateSameChain)
            return candidateSameChain;

        // then in same cluster
        boolean currentSameCluster = beCurrentStart.getSV().getCluster() == beCurrentEnd.getSV().getCluster();
        boolean candidateSameCluster = beCandidateStart.getSV().getCluster() == beCandidateEnd.getSV().getCluster();

        if(currentSameCluster != candidateSameCluster)
            return candidateSameCluster;

        // lastly the nearest to the RNA positions
        double currentPosDiff = (abs(rnaFusion.PositionUp - beCurrentStart.position()) + abs(rnaFusion.PositionDown - beCurrentEnd.position())) * 0.5;
        double candidatePosDiff = (abs(rnaFusion.PositionUp - beCandidateStart.position()) + abs(rnaFusion.PositionDown - beCandidateEnd.position())) * 0.5;

        return candidatePosDiff < currentPosDiff;
    }

    private boolean isViableBreakend(final SvBreakend breakend, long rnaPosition, byte geneStrand, boolean isUpstream)
    {
        boolean requireHigherBreakendPos = isUpstream ? (geneStrand == 1) : (geneStrand == -1);

        long position = breakend.position();

        if(requireHigherBreakendPos)
        {
            // factor in any uncertainty around the precise breakend, eg from homology
            position += breakend.usesStart() ? breakend.getSV().getSvData().startIntervalOffsetEnd() : breakend.getSV().getSvData().endIntervalOffsetEnd();

            return (position >= rnaPosition);
        }
        else
        {
            position += breakend.usesStart() ? breakend.getSV().getSvData().startIntervalOffsetStart() : breakend.getSV().getSvData().endIntervalOffsetStart();

            return (position <= rnaPosition);
        }
    }


    public void close()
    {
        if(mFusionFinder != null)
            mFusionFinder.onCompleted();
    }
}
