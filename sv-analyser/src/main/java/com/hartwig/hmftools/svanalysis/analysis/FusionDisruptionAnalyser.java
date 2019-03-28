package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.EXON_RANK_MAX;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.EXON_RANK_MIN;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.svanalysis.types.RnaFusionData.RNA_SPLICE_TYPE_ONLY_REF;
import static com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser.checkFusionLogic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.svanalysis.annotators.VisualiserWriter;
import com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.types.RnaFusionData;
import com.hartwig.hmftools.svannotation.analysis.SvDisruptionAnalyser;
import com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class FusionDisruptionAnalyser
{
    private SvFusionAnalyser mFusionFinder;
    private SvDisruptionAnalyser mDisruptionFinder;

    private String mSampleId;
    private String mOutputDir;
    private SvGeneTranscriptCollection mGeneTransCollection;
    private Map<String, List<SvBreakend>> mChrBreakendMap;

    private boolean mSkipFusionCheck;
    private List<GeneFusion> mFusions;
    private Map<String, List<RnaFusionData>> mSampleRnaData;

    private VisualiserWriter mVisWriter;
    private BufferedWriter mRnaWriter;

    public static final String SAMPLE_RNA_FILE = "sample_rna_file";
    public static final String SKIP_FUSION_OUTPUT = "skip_fusion_output";

    private static final Logger LOGGER = LogManager.getLogger(FusionDisruptionAnalyser.class);

    public FusionDisruptionAnalyser()
    {
        mFusionFinder = null;
        mDisruptionFinder = null;
        mGeneTransCollection = new SvGeneTranscriptCollection();
        mOutputDir = "";
        mFusions = Lists.newArrayList();
        mSkipFusionCheck = false;
        mVisWriter = null;

        mSampleRnaData = Maps.newHashMap();
        mRnaWriter = null;

        mChrBreakendMap = null;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_RNA_FILE, true, "Sample RNA data to match");
        options.addOption(SKIP_FUSION_OUTPUT, false, "No fusion search or output");
    }

    public void loadFusionReferenceData(final CommandLine cmdLineArgs, final String outputDir, SvGeneTranscriptCollection ensemblDataCache)
    {
        mOutputDir = outputDir;

        mGeneTransCollection = ensemblDataCache;
        mFusionFinder = new SvFusionAnalyser(cmdLineArgs, ensemblDataCache, mOutputDir);

        mSkipFusionCheck = cmdLineArgs.hasOption(SKIP_FUSION_OUTPUT);

        if (cmdLineArgs.hasOption(SAMPLE_RNA_FILE))
        {
            loadSampleRnaData(cmdLineArgs.getOptionValue(SAMPLE_RNA_FILE));
        }
    }

    public final Set<String> getRnaSampleIds() { return mSampleRnaData.keySet(); }
    public final List<GeneFusion> getFusions() { return mFusions; }
    public void setVisWriter(VisualiserWriter writer) { mVisWriter = writer; }

    public static void setSvGeneData(final List<SvVarData> svList, SvGeneTranscriptCollection geneCollection,
            boolean applyPromotorDistance, boolean selectiveLoading)
    {
        int upstreamDistance = applyPromotorDistance ? PRE_GENE_PROMOTOR_DISTANCE : 0;

        if(selectiveLoading)
        {
            // only load transcript info for the genes covered
            Map<String,Boolean> restrictedGeneIds = Maps.newHashMap();

            for (final SvVarData var : svList)
            {
                // isSpecificSV(var);
                for (int be = SVI_START; be <= SVI_END; ++be)
                {
                    if (be == SVI_END && var.isNullBreakend())
                        continue;

                    boolean isStart = isStart(be);

                    geneCollection.populateGeneIdList(restrictedGeneIds, var.chromosome(isStart), var.position(isStart), upstreamDistance);
                }
            }

            geneCollection.loadEnsemblTranscriptData(restrictedGeneIds);
        }

        for(final SvVarData var : svList)
        {
            // isSpecificSV(var);

            List<GeneAnnotation> genesList = Lists.newArrayList();

            for (int be = SVI_START; be <= SVI_END; ++be)
            {
                if (be == SVI_END && var.isNullBreakend())
                    continue;

                boolean isStart = isStart(be);

                genesList.addAll(geneCollection.findGeneAnnotationsBySv(
                        var.dbId(), isStart, var.chromosome(isStart), var.position(isStart), upstreamDistance));
            }

            if (genesList.isEmpty())
                continue;

            List<GeneAnnotation> startGenes = Lists.newArrayList();
            List<GeneAnnotation> endGenes = Lists.newArrayList();

            for (GeneAnnotation gene : genesList)
            {
                gene.setSvData(var.getSvData());

                if (gene.isStart())
                    startGenes.add(gene);
                else
                    endGenes.add(gene);
            }

            // LOGGER.debug("SV({}) matched {} start genes, {} end genes", var.id(), startGenes.size(), endGenes.size());

            var.setGenesList(startGenes, true);
            var.setGenesList(endGenes, false);
        }
    }

    public void run(final String sampleId, final List<SvVarData> svList, final List<SvCluster> clusters, Map<String, List<SvBreakend>> chrBreakendMap)
    {
        mSampleId = sampleId;
        mChrBreakendMap = chrBreakendMap;

        if(!mSkipFusionCheck)
            findFusions(svList, clusters);

        if(!mSampleRnaData.isEmpty())
            assessRnaFusions();
    }

    private void findFusions(final List<SvVarData> svList, final List<SvCluster> clusters)
    {
        if(mSampleId.isEmpty() || mFusionFinder == null)
            return;

        mFusions.clear();

        boolean checkSoloSVs = false;
        boolean checkClusters = true;

        if(checkSoloSVs)
        {
            // always report SVs by themselves
            for (final SvVarData var : svList)
            {
                if (var.isReplicatedSv())
                    continue;

                if (var.isNullBreakend())
                    continue;

                checkFusions(var.getGenesList(true), var.getGenesList(false), var.getCluster(),
                        null, 0, 0);
            }
        }

        if(checkClusters)
        {
            // for now only consider simple SVs and resolved small clusters
            for (final SvCluster cluster : clusters)
            {
                if (cluster.getSvCount() == 1) // simple clusters already checked
                    continue;

                if(cluster.getChains().isEmpty())
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

        // the lower link (index-wise) uses the other end of the link
        // the upper link uses the other end of the link as well

        for (int lpIndex1 = 0; lpIndex1 < linkedPairs.size(); ++lpIndex1)
        {
            SvLinkedPair pair1 = linkedPairs.get(lpIndex1);

            SvVarData lowerSV = pair1.first();

            if(lowerSV.isNullBreakend())
                continue;

            SvBreakend lowerBreakend = pair1.first().getBreakend(!pair1.firstLinkOnStart());

            List<GeneAnnotation> genesListLower = lowerSV.getGenesList(lowerBreakend.usesStart());

            if (genesListLower.isEmpty())
                continue;

            if(lpIndex1 > 0)
            {
                SvLinkedPair prevPair = linkedPairs.get(lpIndex1 - 1);
                checkBreakendExtendsPastGene(lowerBreakend, genesListLower, prevPair.getFirstBreakend().position());

                if (genesListLower.isEmpty())
                    continue;
            }

            long totalLinkLength = 0;

            for (int lpIndex2 = lpIndex1; lpIndex2 < linkedPairs.size(); ++lpIndex2)
            {
                SvLinkedPair pair2 = linkedPairs.get(lpIndex2);

                SvVarData upperSV = pair2.second();

                if(upperSV.isNullBreakend())
                    continue;

                SvBreakend upperBreakend = upperSV.getBreakend(!pair1.secondLinkOnStart());
                List<GeneAnnotation> genesListUpper = upperSV.getGenesList(upperBreakend.usesStart());

                if (genesListUpper.isEmpty())
                    continue;

                if(lpIndex2 < linkedPairs.size() - 1)
                {
                    SvLinkedPair nextPair = linkedPairs.get(lpIndex2 + 1);
                    checkBreakendExtendsPastGene(upperBreakend, genesListUpper, nextPair.getSecondBreakend().position());

                    if (genesListUpper.isEmpty())
                        continue;
                }

                totalLinkLength += pair2.length();

                if (pairTraversesGene(pair2))
                {
                    // if the TI traverses any gene it invalidates any fusion

                    // move ahead to start again at the next index
                    lpIndex1 = lpIndex2;
                    break;
                }

                // test the fusion between these 2 breakends
                LOGGER.info("cluster({}) chain({}) potential fusion: be1({} {}) & be2({} {}) link indices({} -> {}) traversalLen({})",
                        cluster.id(), chain.id(), lowerBreakend.toString(), genesListLower.get(0).GeneName,
                        upperBreakend.toString(), genesListUpper.get(0).GeneName, lpIndex1, lpIndex2, totalLinkLength);

                checkFusions(genesListLower, genesListUpper, cluster, chain, lpIndex2 - lpIndex1 + 1, totalLinkLength);

                // if the upper breakend is in a gene, any further fusions will need to start from this location
                lpIndex1 = lpIndex2;
                break;
            }
        }
    }

    private void checkBreakendExtendsPastGene(SvBreakend breakend, List<GeneAnnotation> genesList, long otherPosition)
    {
        int index = 0;
        while(index < genesList.size())
        {
            GeneAnnotation geneAnnotation = genesList.get(index);
            EnsemblGeneData geneData = mGeneTransCollection.getGeneDataById(geneAnnotation.StableId);

            if(geneData == null)
            {
                genesList.remove(index);
                continue;
            }

            // check whether the other position extends beyond the gene bounds
            if((breakend.orientation() == 1 && otherPosition > geneData.GeneStart)
            || (breakend.orientation() == -1 && otherPosition < geneData.GeneEnd))
            {
                genesList.remove(index);
                continue;
            }

            ++index;
        }
    }

    private boolean pairTraversesGene(SvLinkedPair pair)
    {
        // check traversal of any gene by this section
        long lowerPos = pair.getBreakend(true).position();
        long upperPos = pair.getBreakend(false).position();

        List<EnsemblGeneData> geneDataList = mGeneTransCollection.getChrGeneDataMap().get(pair.chromosome());

        for(EnsemblGeneData geneData : geneDataList)
        {
            if(lowerPos > geneData.GeneEnd)
                continue;

            if(upperPos < geneData.GeneStart)
                break;

            return true;
        }

        return false;
    }

    private void checkFusions(List<GeneAnnotation> breakendGenes1, List<GeneAnnotation> breakendGenes2,
            final SvCluster cluster, final SvChain chain, int links, long chainLength)
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

        if(chain != null)
        {
            LOGGER.info("sample({}) chained fusion: cluster({}) chain({}) links({}) length({})",
                    mSampleId, cluster.id(), chain.id(), links, chainLength);
        }

        String clusterInfo = String.format("%d,%d,%s;%d;%d;%d",
                cluster.id(), cluster.getSvCount(), cluster.getResolvedType(),
                chain != null ? chain.id() : -1, links, chainLength);

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
        final List<RnaFusionData> rnaFusionList = mSampleRnaData.get(mSampleId);

        if (rnaFusionList == null || rnaFusionList.isEmpty())
            return;

        LOGGER.debug("assessing {} RNA fusions", rnaFusionList.size());

        for (final RnaFusionData rnaFusion : rnaFusionList)
        {
            setRnaFusionData(rnaFusion);

            annotateRnaFusions(rnaFusion);

            writeRnaMatchData(mSampleId, rnaFusion);
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
        List<SvBreakend> viableUpBreakends = Lists.newArrayList();
        List<SvBreakend> viableDownBreakends = Lists.newArrayList();
        List<Transcript> viableUpTranscripts = Lists.newArrayList();
        List<Transcript> viableDownTranscripts = Lists.newArrayList();

        // transcripts on the correct side and orientation of the RNA boundary
        List<Transcript> nearUpTranscripts = Lists.newArrayList();
        List<Transcript> nearDownTranscripts = Lists.newArrayList();
        List<SvBreakend> nearUpBreakends = Lists.newArrayList();
        List<SvBreakend> nearDownBreakends = Lists.newArrayList();

        // non-viable transcripts to be used if no others are found
        List<Transcript> genicUpTranscripts = Lists.newArrayList();
        List<Transcript> genicDownTranscripts = Lists.newArrayList();
        List<SvBreakend> genicUpBreakends = Lists.newArrayList();
        List<SvBreakend> genicDownBreakends = Lists.newArrayList();

        boolean isExactRnaExon = rnaFusion.SpliceType.equals(RNA_SPLICE_TYPE_ONLY_REF);

        for(int i = 0; i <= 1 ; ++i)
        {
            boolean isUpstream = (i == 0);
            String chromosome = isUpstream ? rnaFusion.ChrUp : rnaFusion.ChrDown;
            long rnaPosition = isUpstream ? rnaFusion.PositionUp : rnaFusion.PositionDown;
            byte geneStrand = isUpstream ? rnaFusion.StrandUp : rnaFusion.StrandDown;
            List<SvBreakend> viableBreakends = isUpstream ? viableUpBreakends : viableDownBreakends;
            List<SvBreakend> nearBreakends = isUpstream ? nearUpBreakends : nearDownBreakends;
            List<SvBreakend> genicBreakends = isUpstream ? genicUpBreakends : genicDownBreakends;
            List<Transcript> viableTranscripts = isUpstream ? viableUpTranscripts : viableDownTranscripts;
            List<Transcript> nearTranscripts = isUpstream ? nearUpTranscripts : nearDownTranscripts;
            List<Transcript> genicTranscripts = isUpstream ? genicUpTranscripts : genicDownTranscripts;
            String geneName = isUpstream ? rnaFusion.GeneUp : rnaFusion.GeneDown;

            final List<SvBreakend> breakendList = mChrBreakendMap.get(chromosome);

            if(breakendList == null)
                continue;

            for(final SvBreakend breakend : breakendList)
            {
                final SvVarData var = breakend.getSV();

                if(var.isNoneSegment())
                    continue;

                // isSpecificSV(var);

                // check whether breakend falls in genic region
                List<GeneAnnotation> genesList = var.getGenesList(breakend.usesStart())
                        .stream()
                        .filter(x -> x.GeneName.equals(geneName))
                        .collect(Collectors.toList());

                if(genesList.isEmpty())
                    continue;

                // check that breakend has correct orientation and position relative to RNA breakend
                boolean correctLocation = isViableBreakend(breakend, rnaPosition, geneStrand, isUpstream);

                // check whether any of the breakend's transcripts falls within the nearest exon of the RNA fusion breakpoint
                for(final Transcript trans : genesList.get(0).transcripts())
                {
                    if(trans.isCanonical())
                    {
                        if(correctLocation)
                        {
                            nearBreakends.add(breakend);
                            nearTranscripts.add(trans);
                        }
                        else
                        {
                            genicBreakends.add(breakend);
                            genicTranscripts.add(trans);
                        }
                    }

                    if(correctLocation && mFusionFinder.isTranscriptBreakendViableForRnaBoundary(
                            trans, isUpstream,  breakend.position(), rnaPosition, isExactRnaExon))
                    {
                        viableBreakends.add(breakend);
                        viableTranscripts.add(trans);
                        break;
                    }
                }
            }
        }

        LOGGER.debug("rna fusion({}) breakend matches: upstream(viable={} near={} genic={}) downstream(viable={} near={} genic={})",
                rnaFusion.Name, viableUpBreakends.size(), nearUpBreakends.size(), genicUpBreakends.size(),
                viableDownBreakends.size(), nearDownBreakends.size(), genicDownBreakends.size());

        // run them through fusion logic (ie a pair of breakend lists), but don't require phase matching
        if(!viableUpBreakends.isEmpty() && !viableDownBreakends.isEmpty())
        {
            GeneFusion topCandidateFusion = null;
            SvBreakend topUpBreakend = null;
            SvBreakend topDownBreakend = null;

            for (int i = 0; i < viableUpBreakends.size(); ++i)
            {
                final SvBreakend upBreakend = viableUpBreakends.get(i);
                final Transcript upTrans = viableUpTranscripts.get(i);

                if(upBreakend.getSV().isNullBreakend())
                    continue;

                for (int j = 0; j < viableDownBreakends.size(); ++j)
                {
                    final SvBreakend downBreakend = viableDownBreakends.get(j);
                    final Transcript downTrans = viableDownTranscripts.get(j);

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

                        LOGGER.debug("rnaFusion({}) first pair({} & {})", rnaFusion.Name, upBreakend.toString(), downBreakend.toString());
                    }
                }
            }

            if(topCandidateFusion != null)
            {
                rnaFusion.setTranscriptData(
                        true, topCandidateFusion.upstreamTrans(), topUpBreakend,
                        true, true,  0);

                rnaFusion.setTranscriptData(
                        false, topCandidateFusion.downstreamTrans(), topDownBreakend,
                        true, true,0);

                rnaFusion.setViableFusion(topCandidateFusion.viable() && topCandidateFusion.phaseMatched());
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
                boolean correctLocation = false;

                // use the viable transcripts if present, otherwise the nearest
                if(isUpstream)
                {
                    if(!viableUpTranscripts.isEmpty())
                    {
                        isViable = true;
                        correctLocation = true;
                        transcriptList = viableUpTranscripts;
                        breakendList = viableUpBreakends;
                    }
                    else if(!nearUpTranscripts.isEmpty())
                    {
                        correctLocation = true;
                        transcriptList = nearUpTranscripts;
                        breakendList = nearUpBreakends;
                    }
                    else
                    {
                        transcriptList = genicUpTranscripts;
                        breakendList = genicUpBreakends;
                    }
                }
                else
                {
                    if(!viableDownTranscripts.isEmpty())
                    {
                        isViable = true;
                        correctLocation = true;
                        transcriptList = viableDownTranscripts;
                        breakendList = viableDownBreakends;
                    }
                    else if(!nearDownTranscripts.isEmpty())
                    {
                        correctLocation = true;
                        transcriptList = nearDownTranscripts;
                        breakendList = nearDownBreakends;
                    }
                    else
                    {
                        transcriptList = genicDownTranscripts;
                        breakendList = genicDownBreakends;
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
                    int exonsSkipped = 0;

                    if(!isViable)
                    {
                        // for non-viable breakends, provide the exons skipped count
                        final String geneId = isUpstream ? rnaFusion.GeneUp : rnaFusion.GeneDown;
                        final int rnaExonData[] = mGeneTransCollection.getExonRankings(geneId, rnaPosition);
                        final int svPosExonData[] = mGeneTransCollection.getExonRankings(geneId, closestBreakend.position());

                        exonsSkipped = abs(rnaExonData[EXON_RANK_MIN] - svPosExonData[EXON_RANK_MIN]);
                    }

                    rnaFusion.setTranscriptData(isUpstream, closestTrans, closestBreakend, isViable, correctLocation, exonsSkipped);

                    LOGGER.debug("rnaFusion({}) {} closest breakend({}) distance({})",
                            rnaFusion.Name, isUpstream ? "up" :"down", closestBreakend.toString(), closestDistance);
                }
            }
        }

        rnaFusion.setFusionClusterChainInfo();
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

        SvVarData currentStartSV = beCurrentStart.getSV();
        SvVarData currentEndSV = beCurrentEnd.getSV();
        SvVarData candidateStartSV = beCandidateStart.getSV();
        SvVarData candidateEndSV = beCandidateEnd.getSV();

        // give priority to same SV
        boolean currentSameSV = currentStartSV == currentEndSV;
        boolean candidateSameSV = candidateStartSV == candidateEndSV;

        if(currentSameSV != candidateSameSV)
            return candidateSameSV;

        // then whether chained
        boolean currentSameCluster = currentStartSV.getCluster() == currentEndSV.getCluster();
        boolean candidateSameCluster = candidateStartSV.getCluster() == candidateEndSV.getCluster();

        if(currentSameCluster != candidateSameCluster)
        {
            LOGGER.debug("current pair({} & {}) clusters({} & {}), candidate pair({} & {}) clusters({} & {})",
                    currentStartSV.id(), currentEndSV.id(), currentStartSV.getCluster().id(), currentEndSV.getCluster().id(),
                    candidateStartSV.id(), candidateEndSV.id(), candidateStartSV.getCluster().id(), candidateEndSV.getCluster().id());

            return candidateSameCluster;
        }

        if(currentSameCluster && candidateSameCluster)
        {
            // check whether one pair is in the same chain and the other not
            SvChain currentMatchingChain = currentStartSV.getCluster().findSameChainForSVs(currentStartSV, currentEndSV);

            SvChain candidateMatchingChain = candidateStartSV.getCluster().findSameChainForSVs(candidateStartSV, candidateEndSV);

            LOGGER.debug("current pair({} & {}) clusters({} chain={}), candidate pair({} & {}) clusters({} chain={})",
                    currentStartSV.id(), currentEndSV.id(), currentStartSV.getCluster().id(),
                    currentMatchingChain != null ? currentMatchingChain.id() : "diff",
                    candidateStartSV.id(), candidateEndSV.id(), candidateStartSV.getCluster().id(),
                    candidateMatchingChain != null ? candidateMatchingChain.id() : "diff");

            if(currentMatchingChain != null && candidateMatchingChain == null)
                return false;
            if(currentMatchingChain == null && candidateMatchingChain != null)
                return true;
        }

        // otherwise revert to whichever positions are closest to the RNA breakends

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

    public void setRnaFusionData(final RnaFusionData rnaFusion)
    {
        int[] transUpExonData = mGeneTransCollection.getExonRankings(rnaFusion.GeneUp, rnaFusion.PositionUp);
        rnaFusion.setExonUpRank(transUpExonData[EXON_RANK_MIN], transUpExonData[EXON_RANK_MAX]);

        transUpExonData = mGeneTransCollection.getExonRankings(rnaFusion.GeneDown, rnaFusion.PositionDown);
        rnaFusion.setExonDownRank(transUpExonData[EXON_RANK_MIN], transUpExonData[EXON_RANK_MAX]);
    }


    public void writeRnaMatchData(final String sampleId, final RnaFusionData rnaFusion)
    {
        try
        {
            if(mRnaWriter == null)
            {
                String outputFilename = mOutputDir;

                outputFilename += "RNA_MATCH_DATA.csv";

                mRnaWriter = createBufferedWriter(outputFilename, false);

                mRnaWriter.write("SampleId,FusionName,GeneUp,GeneDown,ViableFusion");

                mRnaWriter.write(",SvIdUp,ChrUp,PosUp,RnaPosUp,OrientUp,StrandUp,TypeUp,ClusterInfoUp");
                mRnaWriter.write(",TransViableUp,TransValidLocUp,TransIdUp,ExonsSkippedUp,RegionTypeUp,CodingTypeUp,ExonUp,DisruptiveUp,DistancePrevUp");

                mRnaWriter.write(",SvIdDown,ChrDown,PosDown,RnaPosDown,OrientDown,StrandDown,TypeDown,ClusterInfoDown");
                mRnaWriter.write(",TransViableDown,TransValidLocDown,TransIdDown,ExonsSkippedDown,RegionTypeDown,CodingTypeDown,ExonDown,DisruptiveDown,DistancePrevDown");

                mRnaWriter.write(",ChainInfo,JunctionReadCount,SpanningFragCount,SpliceType");
                mRnaWriter.write(",ExonMinRankUp,ExonMaxRankUp,ExonMinRankDown,ExonMaxRankDown");

                mRnaWriter.newLine();
            }

            BufferedWriter writer = mRnaWriter;

            writer.write(String.format("%s,%s,%s,%s,%s",
                    sampleId, rnaFusion.Name, rnaFusion.GeneUp, rnaFusion.GeneDown, rnaFusion.isViableFusion()));

            final Transcript transUp = rnaFusion.getTrans(true);

            if(transUp != null)
            {
                writer.write(String.format(",%d,%s,%d,%d,%d,%d,%s,%s",
                        transUp.parent().id(), transUp.parent().chromosome(), transUp.parent().position(), rnaFusion.PositionUp,
                        transUp.parent().orientation(), transUp.parent().Strand, transUp.parent().type(),
                        rnaFusion.getClusterInfo(true)));

                writer.write(String.format(",%s,%s,%s,%d,%s,%s,%d,%s,%d",
                        rnaFusion.isTransViable(true), rnaFusion.isTransCorrectLocation(true),
                        transUp.StableId, rnaFusion.getExonsSkipped(true),
                        transUp.regionType(), transUp.codingType(),
                        transUp.exonUpstream(), transUp.isDisruptive(), transUp.exonDistanceUp()));
            }
            else
            {
                writer.write(String.format(",%s,%s,%d,%d,%d,%d,%s,%s",
                        "", rnaFusion.ChrUp, 0, rnaFusion.PositionUp,
                        0, rnaFusion.StrandUp, "", ""));

                writer.write(String.format(",%s,%s,,,,,,,",
                        rnaFusion.isTransViable(true), rnaFusion.isTransCorrectLocation(true)));
            }

            final Transcript transDown = rnaFusion.getTrans(false);

            if(transDown != null)
            {
                writer.write(
                        String.format(",%d,%s,%d,%d,%d,%d,%s,%s",
                                transDown.parent().id(), transDown.parent().chromosome(), transDown.parent().position(), rnaFusion.PositionDown,
                                transDown.parent().orientation(), transDown.parent().Strand, transDown.parent().type(),
                                rnaFusion.getClusterInfo(false)));

                writer.write(
                        String.format(",%s,%s,%s,%d,%s,%s,%d,%s,%d",
                                rnaFusion.isTransViable(false), rnaFusion.isTransCorrectLocation(false),
                                transDown.StableId, rnaFusion.getExonsSkipped(false),
                                transDown.regionType(), transDown.codingType(),
                                transDown.exonDownstream(), transDown.isDisruptive(), transDown.exonDistanceUp()));
            }
            else
            {
                writer.write(String.format(",%s,%s,%d,%d,%d,%d,%s,%s",
                        "", rnaFusion.ChrDown, 0, rnaFusion.PositionDown,
                        0, rnaFusion.StrandDown, "", ""));

                writer.write(String.format(",%s,%s,,,,,,,",
                        rnaFusion.isTransViable(false), rnaFusion.isTransCorrectLocation(false)));
            }

            writer.write(String.format(",%s,%d,%d,%s",
                    !rnaFusion.getChainInfo().isEmpty() ? rnaFusion.getChainInfo() : "0;0",
                    rnaFusion.JunctionReadCount, rnaFusion.SpanningFragCount, rnaFusion.SpliceType));

            writer.write(String.format(",%d,%d,%d,%d",
                    rnaFusion.exonMinRankUp(), rnaFusion.exonMaxRankUp(), rnaFusion.exonMinRankDown(), rnaFusion.exonMaxRankDown()));

            writer.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing RNA match data: {}", e.toString());
        }
    }

    public final Map<String, List<RnaFusionData>> getSampleRnaData() { return mSampleRnaData; }
    public final List<RnaFusionData> getSampleRnaData(final String sampleId) { return mSampleRnaData.get(sampleId); }

    private static int COL_SAMPLEID = 0;
    private static int COL_NAME = 1;
    private static int COL_JUNCT_RC = 2;
    private static int COL_SPAN_RC = 3;
    private static int COL_SPLICE = 4;
    private static int COL_GENE_UP = 5;
    private static int COL_CHR_UP = 7;
    private static int COL_POS_UP = 8;
    private static int COL_STRAND_UP = 9;
    private static int COL_GENE_DOWN = 10;
    private static int COL_CHR_DOWN = 12;
    private static int COL_POS_DOWN = 13;
    private static int COL_STRAND_DOWN = 14;

    private boolean loadSampleRnaData(final String filename)
    {
        if (filename.isEmpty() || !Files.exists(Paths.get(filename)))
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("empty RNA data file({})", filename);
                return false;
            }

            line = fileReader.readLine(); // skip header

            String currentSampleId = "";
            List<RnaFusionData> rnaDataList = Lists.newArrayList();

            while (line != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                // check if still on the same variant
                final String sampleId = items[COL_SAMPLEID];

                if(currentSampleId.isEmpty() || !currentSampleId.equals(sampleId))
                {
                    currentSampleId = sampleId;
                    rnaDataList = Lists.newArrayList();
                    mSampleRnaData.put(currentSampleId, rnaDataList);
                }

                // check that gene names match Ensembl
                String geneUp = items[COL_GENE_UP];
                String geneDown = items[COL_GENE_DOWN];

                geneUp = checkAlternateGeneName(geneUp);
                geneDown = checkAlternateGeneName(geneDown);

                RnaFusionData rnaData = new RnaFusionData(
                        items[COL_NAME], geneUp, geneDown, items[COL_CHR_UP], items[COL_CHR_DOWN],
                        Long.parseLong(items[COL_POS_UP]), Long.parseLong(items[COL_POS_DOWN]),
                        Byte.parseByte(items[COL_STRAND_UP]), Byte.parseByte(items[COL_STRAND_DOWN]),
                        Integer.parseInt(items[COL_JUNCT_RC]),Integer.parseInt(items[COL_SPAN_RC]), items[COL_SPLICE]);

                rnaDataList.add(rnaData);

                line = fileReader.readLine();
            }

        }
        catch(IOException e)
        {
            LOGGER.warn("failed to load sample RNA data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    private String checkAlternateGeneName(final String geneName)
    {
        if(geneName.equals("AC005152.2"))
            return "SOX9-AS1";

        if(geneName.equals("AC016683.6"))
            return "PAX8-AS1";

        if(geneName.equals("AC007092.1"))
            return "LINC01122";

        if(geneName.equals("C10ORF112"))
            return "MALRD1";

        if(geneName.equals("C5orf50"))
            return "SMIM23";

        if(geneName.equals("C10orf68") || geneName.equals("C10orf112"))
            return geneName.toUpperCase();

        if(geneName.equals("C17orf76-AS1"))
            return "FAM211A-AS1";

        if(geneName.equals("IGH-@"))
            return "";

        if(geneName.equals("MKLN1-AS1"))
            return "LINC-PINT";

        if(geneName.equals("PHF15"))
            return "JADE2";

        if(geneName.equals("PHF17"))
            return "JADE1";

        if(geneName.equals("RP11-134P9.1"))
            return "LINC01136";

        if(geneName.equals("RP11-973F15.1"))
            return "LINC01151";

        if(geneName.equals("RP11-115K3.2"))
            return "YWHAEP7";

        if(geneName.equals("RP11-3B12.1"))
            return "POT1-AS1";

        if(geneName.equals("RP11-199O14.1"))
            return "CASC20";

        if(geneName.equals("RP11-264F23.3"))
            return "CCND2-AS1";

        if(geneName.equals("RP11-93L9.1"))
            return "LINC01091";

        return geneName;
    }

    public void close()
    {
        if(mFusionFinder != null)
            mFusionFinder.onCompleted();

        closeBufferedWriter(mRnaWriter);
    }
}
