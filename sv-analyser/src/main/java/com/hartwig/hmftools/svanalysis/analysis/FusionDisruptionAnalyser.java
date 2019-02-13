package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.PRE_GENE_PROMOTOR_DISTANCE;

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
        mSvFusionList = Lists.newArrayList();
        mVisWriter = null;
        mChrBreakendMap = null;
    }

    public final SvGeneTranscriptCollection getGeneTranscriptCollection() { return mEnsemblDataCache; }

    public void loadFusionReferenceData(final CommandLine cmdLineArgs, final String outputDir, final String ensemblDataDir)
    {
        mOutputDir = outputDir;

        mFusionFinder = new SvFusionAnalyser(cmdLineArgs, mEnsemblDataCache, mOutputDir);

        List<HmfTranscriptRegion> transcriptRegions = HmfGenePanelSupplier.allGeneList37();
        mChromosomeTranscriptMap = Multimaps.fromRegions(transcriptRegions);

        if(!ensemblDataDir.isEmpty())
        {
            mEnsemblDataCache.setDataPath(ensemblDataDir);
            mEnsemblDataCache.loadEnsemblData();
        }
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

    private static String CHECK_VAR_ID = "";
    // private static String CHECK_VAR_ID = "4718976";
    private static int CHECK_CLUSTER_ID = -1;
    // private static int CHECK_CLUSTER_ID = 94;

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
        findFusions(svList, clusters); // skipped until chained fusion logic extended

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

            if(var.id().equals(CHECK_VAR_ID))
            {
                LOGGER.debug("specific var({})", var.posId());
            }

            checkFusions(var.getGenesList(true), var.getGenesList(false), var.getCluster());
        }

        boolean checkClusters = false;

        if(checkClusters)
        {
            int maxClusterSize = 50;

            // for now only consider simple SVs and resolved small clusters
            for (final SvCluster cluster : clusters)
            {
                if (cluster.id() == CHECK_CLUSTER_ID)
                {
                    LOGGER.debug("specific cluster");
                }

                if (cluster.getCount() == 1) // simple clusters already checked
                    continue;

                // if(cluster.hasReplicatedSVs() || !cluster.isFullyChained() || cluster.getTypeCount(SGL) > 0)
                //    continue;

                if (cluster.getUniqueSvCount() > maxClusterSize)
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

        mFusions.addAll(fusions);

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

        String clusterInfo = String.format("%d,%d,%s", cluster.id(), cluster.getUniqueSvCount(), cluster.getResolvedType());

        mFusionFinder.writeFusions(fusions, mSampleId, clusterInfo, true);

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

        for(final RnaFusionData rnaFusion : rnaFusionList)
        {
            mFusionFinder.setRnaFusionData(rnaFusion);

            // find all SVs with breakends in either gene, take all possible pairings
            // must be correctly positioned before/after rna exon position

            List<SvBreakend> upstreamBreakends = Lists.newArrayList();
            List<SvBreakend> downstreamBreakends = Lists.newArrayList();

            for(int i = 0; i <= 1 ; ++i)
            {
                boolean isUpstream = (i == 0);
                String chromosome = isUpstream ? rnaFusion.ChrUp : rnaFusion.ChrDown;
                long rnaPosition = isUpstream ? rnaFusion.PositionUp : rnaFusion.PositionDown;
                byte geneStrand = isUpstream ? rnaFusion.StrandUp : rnaFusion.StrandDown;
                List<SvBreakend> viableBreakends = isUpstream ? upstreamBreakends : downstreamBreakends;
                String geneName = isUpstream ? rnaFusion.GeneUp : rnaFusion.GeneDown;

                final List<SvBreakend> breakendList = mChrBreakendMap.get(chromosome);

                if(breakendList == null)
                    continue;

                for(final SvBreakend breakend : breakendList)
                {
                    if(breakend.getSV().isNoneSegment())
                        continue;

                    // check breakend falls in genic region
                    boolean inGene = breakend.getSV().getGenesList(breakend.usesStart()).stream()
                            .filter(x -> x.GeneName.equals(geneName))
                            .count() > 0;

                    if(!inGene)
                        continue;

                    if(isViableBreakend(breakend, rnaPosition, geneStrand, isUpstream))
                    {
                        viableBreakends.add(breakend);
                    }
                }
            }

            // run them through fusion logic (ie a pair of breakend lists), but don't require phase matching
            if(!upstreamBreakends.isEmpty() && !downstreamBreakends.isEmpty())
            {
                GeneFusion topCandidateFusion = null;
                SvBreakend topUpBreakend = null;
                SvBreakend topDownBreakend = null;

                for (final SvBreakend upBreakend : upstreamBreakends)
                {
                    List<GeneAnnotation> upGenesList = upBreakend.getSV().getGenesList(upBreakend.usesStart())
                            .stream()
                            .filter(x -> x.GeneName.equals(rnaFusion.GeneUp))
                            .collect(Collectors.toList());

                    for (final SvBreakend downBreakend : downstreamBreakends)
                    {
                        List<GeneAnnotation> downGenesList = downBreakend.getSV().getGenesList(downBreakend.usesStart())
                                .stream()
                                .filter(x -> x.GeneName.equals(rnaFusion.GeneDown))
                                .collect(Collectors.toList());

                        List<GeneFusion> possibleFusions = mFusionFinder.findFusions(upGenesList, downGenesList, false);

                        if (possibleFusions.isEmpty())
                            continue;

                        // find a reportable fusion if possible
                        GeneFusion possibleFusion = selectPossibleFusion(possibleFusions);

                        if (possibleFusion == null)
                            continue;

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
                    rnaFusion.setMatchedFusion(topCandidateFusion);
            }
            else
            {
                // set whatever breakend in can be found
                Transcript upTrans = null;
                Transcript downTrans = null;

                if(!upstreamBreakends.isEmpty())
                {
                    final List<GeneAnnotation> genes = upstreamBreakends.get(0).getSV().getGenesList(upstreamBreakends.get(0).usesStart());
                    upTrans = !genes.get(0).transcripts().isEmpty() ? genes.get(0).transcripts().get(0) : null;

                }

                if(!downstreamBreakends.isEmpty())
                {
                    final List<GeneAnnotation> genes = downstreamBreakends.get(0).getSV().getGenesList(downstreamBreakends.get(0).usesStart());
                    downTrans = !genes.get(0).transcripts().isEmpty() ? genes.get(0).transcripts().get(0) : null;
                }

                rnaFusion.setBreakends(upTrans, downTrans);
            }

            mFusionFinder.writeRnaMatchData(mSampleId, rnaFusion);
            // SvFusion svFusion = new SvFusion();
        }
    }

    private boolean isCandidateBetter(final GeneFusion currentFusion, final SvBreakend beCurrentStart, final SvBreakend beCurrentEnd,
            final GeneFusion candidateFusion, final SvBreakend beCandidateStart, final SvBreakend beCandidateEnd, final RnaFusionData rnaFusion)
    {
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

    private GeneFusion selectPossibleFusion(final List<GeneFusion> possibleFusions)
    {
        // find a reportable fusion if possible
        GeneFusion possibleFusion = null;
        for(final GeneFusion fusion : possibleFusions)
        {
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
        }

        return possibleFusion;
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
