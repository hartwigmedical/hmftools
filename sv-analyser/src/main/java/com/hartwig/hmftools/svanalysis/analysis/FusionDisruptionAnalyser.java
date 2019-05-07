package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation.isUpstream;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_3P_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_5P_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_BOTH_PROM;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_KNOWN;
import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_NONE;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.isSpecificCluster;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isSpecificSV;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.EXON_RANK_MAX;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.EXON_RANK_MIN;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.svanalysis.types.RnaFusionData.RNA_SPLICE_TYPE_ONLY_REF;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.getExonRankings;
import static com.hartwig.hmftools.svannotation.analysis.SvDisruptionAnalyser.markNonDisruptiveTranscripts;
import static com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser.checkFusionLogic;
import static com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser.validFusionTranscript;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
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
import org.apache.commons.math3.util.Pair;
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
    private List<String> mKnownFusionGenes;

    private boolean mSkipFusionCheck;
    private boolean mLogReportableOnly;
    private List<GeneFusion> mFusions;
    private Map<String, List<RnaFusionData>> mSampleRnaData;
    private PerformanceCounter mPerfCounter;

    private VisualiserWriter mVisWriter;
    private BufferedWriter mRnaWriter;

    private List<String> mRestrictedGenes;
    private boolean mRequirePhaseMatching;
    private boolean mReportKnownFusionData;

    public static final String SAMPLE_RNA_FILE = "sample_rna_file";
    public static final String SKIP_FUSION_OUTPUT = "skip_fusion_output";
    public static final String PRE_GENE_BREAKEND_DISTANCE = "fusion_gene_distance";
    public static final String RESTRICTED_GENE_LIST = "restricted_fusion_genes";
    public static final String NO_PHASE_MATCH_REQD = "no_fusion_phase_match";
    public static final String LOG_REPORTABLE_ONLY = "log_reportable_fusion";
    public static final String LOG_KNOWN_FUSION_DATA = "log_known_fusion_data";

    private static final Logger LOGGER = LogManager.getLogger(FusionDisruptionAnalyser.class);

    public FusionDisruptionAnalyser()
    {
        mFusionFinder = null;
        mDisruptionFinder = null;
        mGeneTransCollection = new SvGeneTranscriptCollection();
        mOutputDir = "";
        mFusions = Lists.newArrayList();
        mSkipFusionCheck = false;
        mLogReportableOnly = false;
        mVisWriter = null;

        mSampleRnaData = Maps.newHashMap();
        mRnaWriter = null;

        mPerfCounter = new PerformanceCounter("Fusions");

        mChrBreakendMap = null;
        mKnownFusionGenes = Lists.newArrayList();
        mRestrictedGenes = Lists.newArrayList();
        mRequirePhaseMatching = true;
        mReportKnownFusionData = false;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_RNA_FILE, true, "Sample RNA data to match");
        options.addOption(SKIP_FUSION_OUTPUT, false, "No fusion search or output");
        options.addOption(PRE_GENE_BREAKEND_DISTANCE, true, "Distance after to a breakend to consider in a gene");
        options.addOption(RESTRICTED_GENE_LIST, true, "Restrict fusion search to specific genes");
        options.addOption(NO_PHASE_MATCH_REQD, false, "Check fusion without requiring phase-matching");
        options.addOption(LOG_REPORTABLE_ONLY, false, "Only write out reportable fusions");
        options.addOption(LOG_KNOWN_FUSION_DATA, false, "Only write out reportable fusions");
    }

    public void initialise(final CommandLine cmdLineArgs, final String outputDir, SvGeneTranscriptCollection ensemblDataCache)
    {
        mOutputDir = outputDir;

        mGeneTransCollection = ensemblDataCache;
        mFusionFinder = new SvFusionAnalyser(cmdLineArgs, ensemblDataCache, mOutputDir);
        populateKnownFusionGenes();

        if(cmdLineArgs != null)
        {
            mSkipFusionCheck = cmdLineArgs.hasOption(SKIP_FUSION_OUTPUT);

            if (!mSkipFusionCheck)
            {
                String annotationHeaders = "PhaseMatched,ClusterId,ClusterCount,ResolvedType,OverlapUp,OverlapDown,ChainInfo";
                String fusionFileName = "SVA_FUSIONS.csv";
                mFusionFinder.initialiseOutputFile(fusionFileName, annotationHeaders);
            }

            if (cmdLineArgs.hasOption(PRE_GENE_BREAKEND_DISTANCE))
            {
                int preGeneBreakendDistance = Integer.parseInt(cmdLineArgs.getOptionValue(PRE_GENE_BREAKEND_DISTANCE));
                PRE_GENE_PROMOTOR_DISTANCE = preGeneBreakendDistance;
            }

            if (cmdLineArgs.hasOption(SAMPLE_RNA_FILE))
            {
                loadSampleRnaData(cmdLineArgs.getOptionValue(SAMPLE_RNA_FILE));
            }

            if(cmdLineArgs.hasOption(RESTRICTED_GENE_LIST))
            {
                String restrictedGenesStr = cmdLineArgs.getOptionValue(RESTRICTED_GENE_LIST);
                mRestrictedGenes = Arrays.stream(restrictedGenesStr.split(";")).collect(Collectors.toList());

                LOGGER.info("restricting fusion genes to: {}", restrictedGenesStr);
            }

            mLogReportableOnly = cmdLineArgs.hasOption(LOG_REPORTABLE_ONLY);
            mReportKnownFusionData = cmdLineArgs.hasOption(LOG_KNOWN_FUSION_DATA);

            if(cmdLineArgs.hasOption(NO_PHASE_MATCH_REQD))
                mRequirePhaseMatching = false;
        }
    }

    public final Set<String> getRnaSampleIds() { return mSampleRnaData.keySet(); }
    public final List<GeneFusion> getFusions() { return mFusions; }
    public void setVisWriter(VisualiserWriter writer) { mVisWriter = writer; }

    // for testing
    public void setHasValidConfigData(boolean toggle) { mFusionFinder.setHasValidConfigData(toggle); }

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

            for (int be = SVI_START; be <= SVI_END; ++be)
            {
                if (be == SVI_END && var.isNullBreakend())
                    continue;

                boolean isStart = isStart(be);

                List<GeneAnnotation> genesList = geneCollection.findGeneAnnotationsBySv(
                        var.dbId(), isStart, var.chromosome(isStart), var.position(isStart), upstreamDistance);

                if (genesList.isEmpty())
                    continue;

                for (GeneAnnotation gene : genesList)
                {
                    gene.setSvData(var.getSvData());

                    int transIndex = 0;
                    while(transIndex < gene.transcripts().size())
                    {
                        Transcript transcript = gene.transcripts().get(transIndex);

                        // only retain transcript which are potential fusion candidates (with exception for canonical)
                        if(!validFusionTranscript(transcript) && !transcript.isCanonical())
                        {
                            gene.transcripts().remove(transIndex);
                        }
                        else
                        {
                            ++transIndex;
                        }
                    }
                }

                var.setGenesList(genesList, isStart);
            }

            // mark any transcripts as not disruptive prior to running any fusion logic
            markNonDisruptiveTranscripts(var.getGenesList(true), var.getGenesList(false));
        }
    }

    public void run(final String sampleId, final List<SvVarData> svList, final List<SvCluster> clusters, Map<String, List<SvBreakend>> chrBreakendMap)
    {
        mPerfCounter.start();

        mSampleId = sampleId;
        mChrBreakendMap = chrBreakendMap;

        if(!mSkipFusionCheck)
            findFusions(svList, clusters);

        if(!mSampleRnaData.isEmpty())
            assessRnaFusions();

        mPerfCounter.stop();
    }

    private void findFusions(final List<SvVarData> svList, final List<SvCluster> clusters)
    {
        if(mSampleId.isEmpty() || mFusionFinder == null)
            return;

        mFusions.clear();

        boolean checkSoloSVs = true;
        boolean checkClusters = true;
        boolean checkKnown = true;

        if(checkSoloSVs)
        {
            finalSingleSVFusions(svList);
        }

        if(checkClusters)
        {
            findChainedFusions(clusters);
        }

        if(checkKnown)
        {
            findUnclusteredKnownFusions(svList);
        }
    }

    private void finalSingleSVFusions(final List<SvVarData> svList)
    {
        // always report SVs by themselves
        for (final SvVarData var : svList)
        {
            if (var.isNullBreakend())
                continue;

            // skip SVs which have been chained or in larger
            if(var.getCluster().getSvCount() > 1 && var.getCluster().findChain(var) != null)
                continue;

            // isSpecificSV(var);
            mFusionFinder.setLogInvalidReasons(isSpecificSV(var));

            List<GeneAnnotation> genesListStart = Lists.newArrayList(var.getGenesList(true));
            List<GeneAnnotation> genesListEnd = Lists.newArrayList(var.getGenesList(false));

            applyGeneRestrictions(genesListStart);
            applyGeneRestrictions(genesListEnd);

            if(genesListStart.isEmpty() || genesListEnd.isEmpty())
                continue;

            List<GeneFusion> fusions = mFusionFinder.findFusions(genesListStart, genesListEnd, true, null);

            if (fusions.isEmpty())
                continue;

            if(mLogReportableOnly)
            {
                fusions = fusions.stream().filter(GeneFusion::reportable).collect(Collectors.toList());
            }

            final SvCluster cluster = var.getCluster();

            // check transcript disruptions
            for (final GeneFusion fusion : fusions)
            {
                String[] disruptedTranscriptsStr = {"", ""};

                for(int i = 0; i <=1; ++i)
                {
                    boolean isUpstream = (i == 0);

                    // look at each gene in turn
                    Transcript transcript = isUpstream ? fusion.upstreamTrans() : fusion.downstreamTrans();
                    GeneAnnotation gene = isUpstream ? fusion.upstreamTrans().parent() : fusion.downstreamTrans().parent();

                    SvBreakend breakend = var.getBreakend(gene.isStart());

                    disruptedTranscriptsStr[i] = checkTranscriptDisruptionInfo(breakend, transcript);
                }

                String clusterInfo = String.format("%s,%d,%d,%s,%s,%s,%s",
                        fusion.phaseMatched(), cluster.id(), cluster.getSvCount(), cluster.getResolvedType(),
                        disruptedTranscriptsStr[0], disruptedTranscriptsStr[1], NO_CHAIN_INFO);

                fusion.setAnnotations(clusterInfo);
                writeFusionData(fusion, cluster);
            }

            mFusions.addAll(fusions);
        }

        mFusionFinder.setLogInvalidReasons(false);
    }

    private void findChainedFusions(final List<SvCluster> clusters)
    {
        // for now only consider simple SVs and resolved small clusters
        for (final SvCluster cluster : clusters)
        {
            if (cluster.getSvCount() == 1) // simple clusters already checked
                continue;

            if (cluster.getChains().isEmpty())
                continue;

            isSpecificCluster(cluster);

            List<GeneFusion> chainFusions = Lists.newArrayList();
            List<GeneFusion> validFusions = Lists.newArrayList();

            for (final SvChain chain : cluster.getChains())
            {
                findChainedFusions(cluster, chain, chainFusions, validFusions);
            }

            if(chainFusions.isEmpty())
                continue;

            // now all fusions have been gathered from this chain, set the reportable one (if any)
            chainFusions.stream().forEach(x -> x.setReportable(false));

            LOGGER.debug("cluster({}) found chained({} valid={}) fusions", cluster.id(), chainFusions.size(), validFusions.size());

            // consider fusions from amongst unique gene-pairings
            List<String> genePairings = Lists.newArrayList();

            for(int i = 0; i < validFusions.size(); ++i)
            {
                GeneFusion fusion = validFusions.get(i);
                String genePair = fusion.upstreamTrans().parent().GeneName + "_" + fusion.downstreamTrans().parent().GeneName;

                if(genePairings.contains(genePair))
                    continue;

                genePairings.add(genePair);

                // gather up all matching fusions
                List<GeneFusion> genePairFusions = Lists.newArrayList();
                genePairFusions.add(fusion);

                for(int j = i+1; j < validFusions.size(); ++j)
                {
                    GeneFusion nextFusion = validFusions.get(j);
                    String nextGenePair = nextFusion.upstreamTrans().parent().GeneName + "_" + nextFusion.downstreamTrans().parent().GeneName;

                    if(nextGenePair.equals(genePair))
                    {
                        genePairFusions.add(nextFusion);
                    }
                }

                // only chained fusions with unterminated ends and valid traversal are considered as reportable
                mFusionFinder.setReportableGeneFusions(genePairFusions);
            }

            for(GeneFusion fusion : chainFusions)
            {
                if (mLogReportableOnly && !fusion.reportable())
                    continue;

                writeFusionData(fusion, cluster);
            }
        }
    }

    private static int FUSION_MAX_CHAIN_LENGTH = 100000;

    private void findChainedFusions(final SvCluster cluster, final SvChain chain,
            List<GeneFusion> chainFusions, List<GeneFusion> validFusions)
    {
        // look for fusions formed by breakends connected in a chain

        // given a chain sAe - sBe - sCe - sDe, where As and De are the open ends of the chain, the following fusions need to be tested:
        // a) each SV in isolation, ie as a single SV
        // b) the left-most / lower breakend of each SV with with right-most / upper breakend of SVs higher up in the chain

        // whenever a linked pair is traversed by a fusion, it cannot touch or traverse genic regions without disrupting the fusion
        final List<SvLinkedPair> linkedPairs = chain.getLinkedPairs();

        for (int lpIndex1 = 0; lpIndex1 <= linkedPairs.size(); ++lpIndex1)
        {
            SvVarData lowerSV = null;
            SvBreakend lowerBreakend = null;

            // the lower link takes the other breakend of the current linked pair's 'first' SV
            // and in order for it to also test the last SV in isolation, also takes the last pair's second (upper) breakend
            if(lpIndex1 < linkedPairs.size())
            {
                SvLinkedPair pair = linkedPairs.get(lpIndex1);
                lowerSV = pair.first();
                lowerBreakend = pair.first().getBreakend(!pair.firstLinkOnStart());
            }
            else
            {
                SvLinkedPair prevPair = linkedPairs.get(lpIndex1 - 1);
                lowerSV = prevPair.second();
                lowerBreakend = prevPair.getSecondBreakend();
            }

            // isSpecificSV(lowerSV);

            if (lowerSV.isNullBreakend())
                continue;

            if(lowerSV.isReplicatedSv())
                lowerSV = lowerSV.getOrigSV();

            List<GeneAnnotation> genesListLower = Lists.newArrayList(lowerSV.getGenesList(lowerBreakend.usesStart()));
            applyGeneRestrictions(genesListLower);

            if (genesListLower.isEmpty())
                continue;

            List<SvLinkedPair> traversedPairs = Lists.newArrayList();

            for (int lpIndex2 = lpIndex1; lpIndex2 <= linkedPairs.size(); ++lpIndex2)
            {
                SvVarData upperSV = null;
                SvBreakend upperBreakend = null;

                // the upper link takes the breakend of the current linked pair's 'second' SV
                // and beyond all the links, it must also test fusions with the chain's upper open breakend
                if(lpIndex2 < linkedPairs.size())
                {
                    SvLinkedPair pair = linkedPairs.get(lpIndex2);
                    upperSV = pair.first();
                    upperBreakend = pair.getFirstBreakend();
                }
                else
                {
                    // at the last link, take the open breakend of the chain
                    upperBreakend = chain.getOpenBreakend(false);

                    if(upperBreakend == null) // can be null for SGLs at end of chain
                        break;

                    upperSV = upperBreakend.getSV();
                }

                if(upperSV.isReplicatedSv())
                    upperSV = upperSV.getOrigSV();

                // isSpecificSV(upperSV);

                List<GeneAnnotation> genesListUpper = Lists.newArrayList(upperSV.getGenesList(upperBreakend.usesStart()));
                applyGeneRestrictions(genesListUpper);

                // if a new linked pair has been traversed to reach this upper breakend, record its length
                // and test whether it traverses any genic region, and if so invalidate the fusion
                if(lpIndex2 > lpIndex1)
                {
                    SvLinkedPair lastPair = linkedPairs.get(lpIndex2 - 1);
                    traversedPairs.add(lastPair);
                }

                if(genesListUpper.isEmpty())
                {
                    // skip past this link and breakend to the next one, keeping the possibility of a fusion with the lower breakend open
                    continue;
                }

                /*
                if(lpIndex2 > lpIndex1)
                {
                    LOGGER.debug("cluster({}) chain({}) testing chained fusion: be1({} {}) & be2({} {}) link indices({} -> {})",
                            cluster.id(), chain.id(), lowerBreakend.toString(), genesListLower.get(0).GeneName,
                            upperBreakend.toString(), genesListUpper.get(0).GeneName, lpIndex1, lpIndex2);
                }
                */

                // test the fusion between these 2 breakends
                if(isRepeatedBreakendPair(lowerBreakend, upperBreakend, chainFusions))
                    continue;

                 boolean logFusionReasons = isSpecificSV(lowerBreakend.getSV()) & isSpecificSV(upperBreakend.getSV());
                 mFusionFinder.setLogInvalidReasons(logFusionReasons);

                List<GeneFusion> fusions = mFusionFinder.findFusions(genesListLower, genesListUpper, true, null);

                if(fusions.isEmpty())
                    continue;

                if (lpIndex2 > lpIndex1)
                {
                    // a chain cannot be an exon-exon fusion, so cull any of these
                    fusions = fusions.stream().filter(x -> !x.isExonic()).collect(Collectors.toList());
                }

                if(mLogReportableOnly)
                {
                    fusions = fusions.stream().filter(x -> x.reportable()).collect(Collectors.toList());

                    if(fusions.isEmpty())
                        continue;
                }

                chainFusions.addAll(fusions);

                for (GeneFusion fusion : fusions)
                {
                    // if the fusion from the upstream gene is on the positive strand, then it will have a fusion direction of +1
                    // whenever it goes through a subsequent linked pair by joining to the first (lower) breakend in the pair
                    int upGeneStrand = fusion.upstreamTrans().parent().Strand;
                    boolean isPrecodingUpstream = fusion.upstreamTrans().preCoding();
                    boolean fusionLowerToUpper = fusion.upstreamTrans().parent().position() == lowerBreakend.position();

                    // check any traversed genes
                    long totalLinkLength = 0;
                    boolean validTraversal = true;
                    boolean allTraversalAssembled = true;

                    for(SvLinkedPair pair : traversedPairs)
                    {
                        totalLinkLength += pair.length();

                        if(pair.isInferred())
                            allTraversalAssembled = false;

                        // if going lower to upper, if the orientation of the first breakend in the pair is opposite to the strand of
                        // the upstream gene, then the fusion direction for that pair is the same as a the upstream gene
                        // otherwise it needs to be switched
                        int fusionDirection = 0;

                        if(fusionLowerToUpper)
                        {
                            fusionDirection = pair.getFirstBreakend().orientation() != upGeneStrand ? upGeneStrand : -upGeneStrand;
                        }
                        else
                        {
                            fusionDirection = pair.getSecondBreakend().orientation() != upGeneStrand ? upGeneStrand : -upGeneStrand;
                        }

                        if(validTraversal && pairTraversesGene(pair, fusionDirection, isPrecodingUpstream))
                            validTraversal = false;
                    }

                    String[] disruptedTranscriptsStr = { "", "" };

                    for (int i = 0; i <= 1; ++i)
                    {
                        boolean isUpstream = (i == 0);

                        // look at each gene in turn
                        Transcript transcript = isUpstream ? fusion.upstreamTrans() : fusion.downstreamTrans();
                        GeneAnnotation gene = isUpstream ? fusion.upstreamTrans().parent() : fusion.downstreamTrans().parent();

                        SvBreakend breakend = lowerBreakend.position() == gene.position() ? lowerBreakend : upperBreakend;

                        boolean isChainEnd = (breakend == lowerBreakend && lpIndex1 == 0)
                                || (breakend == upperBreakend && lpIndex2 == linkedPairs.size());

                        if (isChainEnd)
                        {
                            disruptedTranscriptsStr[i] = checkTranscriptDisruptionInfo(breakend, transcript);
                        }
                        else
                        {
                            // looks within and beyond the chain to check for disruptions
                            int linkIndex = breakend == lowerBreakend ? lpIndex1 - 1 : lpIndex2;
                            disruptedTranscriptsStr[i] = checkTranscriptDisruptionInfo(breakend, transcript, chain, linkIndex);
                        }
                    }

                    int linksCount = lpIndex2 - lpIndex1;

                    String chainInfo = String.format("%d;%d;%d;%s;%s",
                            chain.id(), linksCount, totalLinkLength, validTraversal, allTraversalAssembled);

                    String clusterChainInfo = String.format("%s,%d,%d,%s,%s,%s,%s",
                            fusion.phaseMatched(), cluster.id(), cluster.getSvCount(), cluster.getResolvedType(),
                            disruptedTranscriptsStr[0], disruptedTranscriptsStr[1], chainInfo);

                    fusion.setAnnotations(clusterChainInfo);

                    boolean chainLengthOk = fusion.getKnownFusionType() == REPORTABLE_TYPE_KNOWN || totalLinkLength <= FUSION_MAX_CHAIN_LENGTH;
                    boolean notDisrupted = !isDisrupted(disruptedTranscriptsStr[0]) && !isDisrupted(disruptedTranscriptsStr[1]);

                    if(validTraversal && chainLengthOk && notDisrupted)
                    {
                        validFusions.add(fusion);
                    }
                }
            }
        }

        mFusionFinder.setLogInvalidReasons(false);
    }

    private boolean isRepeatedBreakendPair(final SvBreakend be1, final SvBreakend be2, final List<GeneFusion> fusions)
    {
        for(GeneFusion fusion : fusions)
        {
            if(fusion.upstreamTrans().parent().id() == be1.getSV().dbId() && fusion.downstreamTrans().parent().id() == be2.getSV().dbId())
                return true;
            else if(fusion.upstreamTrans().parent().id() == be2.getSV().dbId() && fusion.downstreamTrans().parent().id() == be1.getSV().dbId())
                return true;
        }

        return false;
    }

    // ChainId,ChainLinks,ChainLength,ValidTraversal,TraversalAssembled
    private static String NO_CHAIN_INFO = "-1;0;0;true;false";
    public static int FCI_CHAIN_ID = 0;
    public static int FCI_CHAIN_LINKS = 1;
    public static int FCI_CHAIN_LENGTH = 2;
    public static int FCI_VALID_TRAVERSAL = 3;
    public static int FCI_TRAV_ASSEMBLY = 4;

    private void applyGeneRestrictions(List<GeneAnnotation> genesList)
    {
        if(mRestrictedGenes.isEmpty())
            return;

        int index = 0;
        while(index < genesList.size())
        {
            GeneAnnotation gene = genesList.get(index);

            if(mRestrictedGenes.contains(gene.GeneName))
                ++index;
            else
                genesList.remove(index);
        }
    }

    // FacingBreakends,HasAssembledLink,ClusterFacingBreakends,TotalBreakends,MinDistance,DisruptedExons,Terminated
    private static String NO_DISRUPTION_INFO = "0;false;0;0;0;0;0";
    public static int FDI_FACING_BES = 0;
    public static int FDI_HAS_ASSEMBLY = 1;
    public static int FDI_FACING_CL_BES = 2;
    public static int FDI_TOTLAL_BES = 3;
    public static int FDI_MIN_DISTANCE = 4;
    public static int FDI_DIS_EXONS = 5;
    public static int FDI_TERMINATED = 6;

    public static boolean isDisrupted(final String disruptionInfo)
    {
        String[] disruptionData = disruptionInfo.split(";");

        if(disruptionData.length != FDI_TERMINATED+1)
            return true;

        // test the exons disrupted and terminated fields
        return !disruptionData[FDI_DIS_EXONS].equals("0") || !disruptionData[FDI_TERMINATED].equals("0");
    }

    private String checkTranscriptDisruptionInfo(final SvBreakend breakend, final Transcript transcript)
    {
        // check all breakends which fall within the bounds of this transcript, including any which are exonic
        List<SvBreakend> breakendList = mChrBreakendMap.get(breakend.chromosome());

        List<TranscriptExonData> exonDataList = mGeneTransCollection.getTranscriptExons(transcript.parent().StableId, transcript.StableId);

        if(exonDataList == null || exonDataList.isEmpty())
            return "";

        int totalBreakends = 0;
        int clusterFacingBreakends = 0;
        int facingBreakends = 0;
        int disruptedExons = 0;
        long minDistance = -1;

        final SvCluster cluster = breakend.getSV().getCluster();

        int index = breakend.getChrPosIndex();

        while(true)
        {
            index += breakend.orientation() == -1 ? +1 : -1;

            if (index < 0 || index >= breakendList.size())
                break;

            SvBreakend nextBreakend = breakendList.get(index);

            // exit once the next breakend extends beyond the gene bounds
            if((breakend.orientation() == 1 && nextBreakend.position() < transcript.TranscriptStart)
            || (breakend.orientation() == -1 && nextBreakend.position() > transcript.TranscriptEnd))
            {
                break;
            }

            ++totalBreakends;

            if (nextBreakend.orientation() != breakend.orientation())
            {
                // skip breakends which cannot be chained by min TI length
                int minTiLength = getMinTemplatedInsertionLength(breakend, nextBreakend);
                long breakendDistance = abs(breakend.position() - nextBreakend.position());

                if(breakendDistance < minTiLength)
                    continue;

                if (minDistance == -1)
                    minDistance = breakendDistance;

                ++facingBreakends;

                if (nextBreakend.getSV().getCluster() == cluster)
                    ++clusterFacingBreakends;
            }
        }

        if(facingBreakends > 0)
        {
            // is any facing breakend assembled?
            final SvLinkedPair tiLink = breakend.getSV().getLinkedPair(breakend.usesStart());
            boolean hasAssembledLink = tiLink != null && tiLink.isAssembled();

            return String.format("%d;%s;%d;%d;%d;%d;%d",
                facingBreakends, hasAssembledLink, clusterFacingBreakends, totalBreakends, minDistance, disruptedExons, 0);
        }
        else
        {
            return NO_DISRUPTION_INFO;
        }
    }

    private String checkTranscriptDisruptionInfo(final SvBreakend breakend, final Transcript transcript, final SvChain chain, int linkIndex)
    {
        // starting with this breakend and working onwards from it in the chain, check for any disruptions to the transcript
        // this includes subsequent links within the same chain and transcript
        SvLinkedPair startPair = chain.getLinkedPairs().get(linkIndex);
        boolean traverseUp = startPair.getFirstBreakend() == breakend; // whether to search up or down the chain

        List<TranscriptExonData> exonDataList = mGeneTransCollection.getTranscriptExons(transcript.parent().StableId, transcript.StableId);

        if(exonDataList == null || exonDataList.isEmpty())
            return "";

        int totalBreakends = 0;
        int facingBreakends = 0;
        int disruptedExons = 0;
        boolean transcriptTerminated = false;
        long minDistance = startPair.length();
        boolean allLinksAssembled = true;

        boolean isUpstream = isUpstream(transcript.parent());

        while(linkIndex >= 0 && linkIndex <= chain.getLinkedPairs().size() - 1)
        {
            SvLinkedPair pair = chain.getLinkedPairs().get(linkIndex);

            if(pair.isInferred())
                allLinksAssembled = false;

            // identify next exon after this TI
            // the breakend's transcript info cannot be used because it faces the opposite way from the fusing breakend
            SvBreakend nextBreakend = traverseUp ? pair.getSecondBreakend() : pair.getFirstBreakend();

            // exit if the next breakend is now past the end of the transcript or if the breakend is down-stream of coding
            if(nextBreakend.orientation() == 1)
            {
                if(nextBreakend.position() > transcript.TranscriptEnd)
                    break;
                else if(!isUpstream && transcript.CodingEnd != null && nextBreakend.position() > transcript.CodingEnd)
                    break;
            }
            else
            {
                if(nextBreakend.position() < transcript.TranscriptStart)
                    break;
                else if(!isUpstream && transcript.CodingStart != null && nextBreakend.position() < transcript.CodingStart)
                    break;
            }

            transcriptTerminated = true;

            ++totalBreakends;
            ++facingBreakends;

            break;

            // no longer check that subsequent links within the same transcript are valid

            /*
            int[] nextBreakendExons = getExonRankings(exonDataList, nextBreakend.position());

            if(nextBreakendExons[EXON_RANK_MIN] == nextBreakendExons[EXON_RANK_MAX])
            {
                // TI ends in an exon which is not permitted
                disruptedExons += max(transcript.ExonMax - nextBreakendExons[EXON_RANK_MAX], 0);
                transcriptTerminated = true;
                break;
            }

            if(nextBreakend.getSV().isNullBreakend())
                break;

            // for now only allow simple SVs that aren't disruptive
            if(!nextBreakend.getSV().isSimpleType())
            {
                disruptedExons += max(transcript.ExonMax - nextBreakendExons[EXON_RANK_MAX], 0);
                transcriptTerminated = true;
                break;
            }

            // otherwise check where this next SV goes
            SvBreakend nextOtherBreakend = nextBreakend.getOtherBreakend();

            // it must remain within the same intronic section to be valid
            int[] nextOtherBreakendExons = getExonRankings(exonDataList, nextOtherBreakend.position());

            if(nextOtherBreakendExons[EXON_RANK_MIN] != nextBreakendExons[EXON_RANK_MIN]
            || nextOtherBreakendExons[EXON_RANK_MAX] != nextBreakendExons[EXON_RANK_MAX])
            {
                disruptedExons += max(nextOtherBreakendExons[EXON_RANK_MIN] - nextBreakendExons[EXON_RANK_MIN],
                    nextOtherBreakendExons[EXON_RANK_MAX] - nextBreakendExons[EXON_RANK_MAX]);

                transcriptTerminated = true;
                break;
            }

            ++totalBreakends;

            linkIndex += traverseUp ? 1 : -1;
            */
        }

        if(facingBreakends > 0)
        {
            return String.format("%d;%s;%d;%d;%d;%d;%d",
                    facingBreakends, allLinksAssembled, facingBreakends, totalBreakends,
                    minDistance, disruptedExons, transcriptTerminated ? 1 : 0);
        }
        else
        {
            return NO_DISRUPTION_INFO;
        }
    }

    private boolean pairTraversesGene(SvLinkedPair pair, int fusionDirection, boolean isPrecodingUpstream)
    {
        // for this pair to not affect the fusion, the section it traverses cannot cross any gene's splice acceptor
        // with the same strand direction unless that is part of a fully traversed non-coding 5' exon

        long lowerPos = pair.getBreakend(true).position();
        long upperPos = pair.getBreakend(false).position();

        List<EnsemblGeneData> geneDataList = mGeneTransCollection.getChrGeneDataMap().get(pair.chromosome());

        for(EnsemblGeneData geneData : geneDataList)
        {
            if(lowerPos > geneData.GeneEnd)
                continue;

            if(upperPos < geneData.GeneStart)
                break;

            if(geneData.Strand == fusionDirection)
            {
                // check whether a splice acceptor is encountered within this window
                List<TranscriptExonData> transExonDataList = mGeneTransCollection.getTransExonData(geneData.GeneId);

                if(transExonDataList == null)
                    continue;

                for(TranscriptExonData exonData : transExonDataList)
                {
                    if(exonData.ExonRank == 1)
                        continue;

                    if((geneData.Strand == 1 && lowerPos <= exonData.ExonStart && upperPos >= exonData.ExonStart)
                    || (geneData.Strand == -1 && lowerPos <= exonData.ExonEnd && upperPos >= exonData.ExonEnd))
                    {
                        // allow an exon to be fully traversed if the upstream transcript is pre-coding
                        if(isPrecodingUpstream && lowerPos <= exonData.ExonStart && upperPos >= exonData.ExonEnd)
                        {
                            if(geneData.Strand == 1 && (exonData.CodingStart == null || upperPos < exonData.CodingStart))
                                continue;
                            else if(geneData.Strand == -1 && (exonData.CodingEnd == null || lowerPos > exonData.CodingEnd))
                                continue;
                        }

                        /*
                        LOGGER.debug("pair({}) fusionDirection({}) traverses splice acceptor({} {}) exon(rank{} pos={})",
                                pair.toString(), fusionDirection, geneData.GeneName, exonData.TransName,
                                exonData.ExonRank, exonData.ExonStart, exonData.ExonEnd);
                        */

                        return true;
                    }
                }
            }
        }

        return false;
    }

    private void populateKnownFusionGenes()
    {
        if(mFusionFinder.getKnownFusionsModel() == null)
            return;

        for(Pair<String,String> genePair : mFusionFinder.getKnownFusionsModel().fusions().keySet())
        {
            if(!mKnownFusionGenes.contains(genePair.getFirst()))
                mKnownFusionGenes.add(genePair.getFirst());

            if(!mKnownFusionGenes.contains(genePair.getSecond()))
                mKnownFusionGenes.add(genePair.getSecond());
        }
    }

    public static final String INVALID_REASON_UNCLUSTERED = "Unclustered";
    public static final String INVALID_REASON_UNCHAINED = "Unchained";

    private void findUnclusteredKnownFusions(final List<SvVarData> svList)
    {
        // look for unclustered or unchained SVs with breakends in known fusion pairs
        Map<SvBreakend, GeneAnnotation> knownFusionBreakends = Maps.newHashMap();

        for(final SvVarData var :svList)
        {
            for (int be = SVI_START; be <= SVI_END; ++be)
            {
                boolean isStart = isStart(be);

                List<GeneAnnotation> genesList = Lists.newArrayList(var.getGenesList(isStart));

                for(GeneAnnotation gene : genesList)
                {
                    if(mKnownFusionGenes.contains(gene.GeneName))
                    {
                        knownFusionBreakends.put(var.getBreakend(isStart), gene);
                        break;
                    }
                }
            }
        }

        for (Map.Entry<SvBreakend, GeneAnnotation> entry1 : knownFusionBreakends.entrySet())
        {
            SvBreakend be1 = entry1.getKey();
            GeneAnnotation gene1 = entry1.getValue();
            List<GeneAnnotation> genes1 = Lists.newArrayList(gene1);

            for (Map.Entry<SvBreakend, GeneAnnotation> entry2 : knownFusionBreakends.entrySet())
            {
                SvBreakend be2 = entry2.getKey();
                GeneAnnotation gene2 = entry2.getValue();
                List<GeneAnnotation> genes2 = Lists.newArrayList(gene2);

                if(be1 == be2)
                    continue;

                if(!areKnownFusionGenes(gene1.GeneName, gene2.GeneName))
                    continue;

                if (hasFusionGenePair(gene1, gene2))
                    continue;

                if(be1.getSV().getCluster() == be2.getSV().getCluster())
                {
                    // skip if already in the same chain
                    SvCluster cluster = be1.getSV().getCluster();

                    if(cluster.findSameChainForSVs(be1.getSV(), be2.getSV()) != null)
                        continue;
                }

                isSpecificSV(be1.getSV());

                List<GeneFusion> fusions = mFusionFinder.findFusions(genes1, genes2, true, null);

                if (fusions.isEmpty())
                    continue;

                if(mLogReportableOnly)
                {
                    fusions = fusions.stream().filter(GeneFusion::reportable).collect(Collectors.toList());
                }

                // check transcript disruptions
                for (final GeneFusion fusion : fusions)
                {
                    String[] disruptedTranscriptsStr = {"", ""};

                    for(int i = 0; i <=1; ++i)
                    {
                        boolean isUpstream = (i == 0);

                        // look at each gene in turn
                        Transcript transcript = isUpstream ? fusion.upstreamTrans() : fusion.downstreamTrans();
                        GeneAnnotation gene = isUpstream ? fusion.upstreamTrans().parent() : fusion.downstreamTrans().parent();

                        SvBreakend breakend = gene.GeneName.equals(gene1.GeneName) ? be1 : be2;

                        disruptedTranscriptsStr[i] = checkTranscriptDisruptionInfo(breakend, transcript);
                    }

                    final SvCluster cluster1 = be1.getSV().getCluster();
                    final SvCluster cluster2 = be2.getSV().getCluster();

                    String combinedClusterInfo = String.format("Unclustered;%d;%s;%d;%s",
                            cluster1.getSvCount(), cluster1.getResolvedType(), cluster2.getSvCount(), cluster2.getResolvedType());

                    // PhaseMatched,ClusterId,ClusterCount,ResolvedType,OverlapUp,OverlapDown,ChainInfo
                    String clusterInfo = String.format("%s,%d,%d,%s,%s,%s,%s",
                            fusion.phaseMatched(), cluster1.id(), cluster2.id(), combinedClusterInfo,
                            disruptedTranscriptsStr[0], disruptedTranscriptsStr[1], NO_CHAIN_INFO);

                    fusion.setAnnotations(clusterInfo);
                    writeFusionData(fusion, cluster1);
                }

                mFusions.addAll(fusions);
            }
        }

        if(mReportKnownFusionData)
        {
            logKnownFusionData(knownFusionBreakends);
        }
    }

    private void logKnownFusionData(Map<SvBreakend, GeneAnnotation> knownFusionBreakends)
    {
        for (Map.Entry<SvBreakend, GeneAnnotation> entry1 : knownFusionBreakends.entrySet())
        {
            SvBreakend be1 = entry1.getKey();
            GeneAnnotation gene1 = entry1.getValue();
            List<GeneAnnotation> genes1 = Lists.newArrayList(gene1);

            // skip SVs where all transcripts are non-disruptive
            if(!gene1.hasAnyDisruptiveTranscript())
                continue;

            for (Map.Entry<SvBreakend, GeneAnnotation> entry2 : knownFusionBreakends.entrySet())
            {
                SvBreakend be2 = entry2.getKey();
                GeneAnnotation gene2 = entry2.getValue();
                List<GeneAnnotation> genes2 = Lists.newArrayList(gene2);

                if (be1 == be2)
                    continue;

                if(!gene2.hasAnyDisruptiveTranscript())
                    continue;

                if (!areKnownFusionGenes(gene1.GeneName, gene2.GeneName))
                    continue;

                if (hasFusionGenePair(gene1, gene2)) // already captured even if not valid
                    continue;

                // possible reasons for not forming a valid fusions are:
                // 1. Unclustered or unchained SVs
                // 2. Incorrect orientations
                // 3. Invalid transcript locations (eg post-coding)
                // 4. Unphased

                List<String> invalidReasons = Lists.newArrayList();

                if (be1.getSV().getCluster() != be2.getSV().getCluster())
                {
                    invalidReasons.add(INVALID_REASON_UNCLUSTERED);
                }
                else
                {
                    SvCluster cluster = be1.getSV().getCluster();

                    if(cluster.findSameChainForSVs(be1.getSV(), be2.getSV()) == null)
                    {
                        invalidReasons.add(INVALID_REASON_UNCHAINED);
                    }
                }

                List<GeneFusion> fusions = mFusionFinder.findFusions(genes1, genes2, true, invalidReasons);

                /*
                if (!fusions.isEmpty())
                {
                    LOGGER.warn("sample({}) genes({} & {}) SVs({} & {}) found fusion, should have previously",
                            mSampleId, gene1.GeneName, gene2.GeneName, be1.getSV().id(), be2.getSV().id());
                }
                */

                String invalidReasonStr = "";

                for(String reason : invalidReasons)
                {
                    invalidReasonStr = appendStr(invalidReasonStr, reason, ';');
                }

                // SampleId,GeneUp,GeneDown,SvId1,SvId2,InvalidReasons
                String fusionInfo = String.format("%s,%s,%s,%s,%s,%s",
                        mSampleId, gene1.GeneName, gene2.GeneName, be1.getSV().id(), be2.getSV().id(), invalidReasonStr);

                LOGGER.info("KNOWN_FUSION_DATA: {}", fusionInfo);
            }
        }
    }

    private boolean hasFusionGenePair(final GeneAnnotation geneUp, final GeneAnnotation geneDown)
    {
        for(final GeneFusion fusion : mFusions)
        {
            final GeneAnnotation upGene = fusion.upstreamTrans().parent();
            final GeneAnnotation downGene = fusion.downstreamTrans().parent();

            if(upGene.GeneName.equals(geneUp.GeneName) && downGene.GeneName.equals(geneDown.GeneName))
                return true;
        }

        return false;
    }

    private boolean areKnownFusionGenes(final String geneUp, final String geneDown)
    {
        if(mFusionFinder.getKnownFusionsModel() == null)
            return false;

        for(Pair<String,String> genePair : mFusionFinder.getKnownFusionsModel().fusions().keySet())
        {
            if(genePair.getFirst().equals(geneUp) && genePair.getSecond().equals(geneDown))
                return true;
        }

        return false;
    }

    private void writeFusionData(GeneFusion fusion, final SvCluster cluster)
    {
        mFusions.add(fusion);

        if (fusion.reportable() || !mLogReportableOnly)
        {
            mFusionFinder.writeFusionData(fusion, mSampleId);
        }

        if(fusion.reportable() && mVisWriter != null)
        {
            mVisWriter.addGeneExonData(cluster.id(),
                    fusion.upstreamTrans().parent().StableId, fusion.upstreamTrans().parent().GeneName,
                    fusion.upstreamTrans().StableId, fusion.upstreamTrans().parent().chromosome(), "FUSION");

            mVisWriter.addGeneExonData(cluster.id(),
                    fusion.downstreamTrans().parent().StableId, fusion.downstreamTrans().parent().GeneName,
                    fusion.downstreamTrans().StableId, fusion.downstreamTrans().parent().chromosome(), "FUSION");
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

        // move from consideration to de-link RNA data from SV types
        mSampleRnaData.remove(mSampleId);
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

                    GeneFusion possibleFusion = checkFusionLogic(upTrans, downTrans, false, null);

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

                rnaFusion.setViableFusion(topCandidateFusion.viable(), topCandidateFusion.phaseMatched());
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
                        String geneId = closestTrans.parent().StableId;
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
        if(currentFusion.viable() != candidateFusion.viable())
        {
            return candidateFusion.viable();
        }

        if(currentFusion.phaseMatched() != candidateFusion.phaseMatched())
        {
            return candidateFusion.phaseMatched();
        }

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
        if(areKnownFusionGenes(rnaFusion.GeneUp, rnaFusion.GeneDown))
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_KNOWN);
            return;
        }

        KnownFusionsModel refFusionData = mFusionFinder.getKnownFusionsModel();

        if(refFusionData == null)
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_NONE);
            return;
        }

        boolean fivePrimeProm = refFusionData.fivePrimePromiscuousMatch(Lists.newArrayList(rnaFusion.GeneUp));
        boolean threePrimeProm = refFusionData.threePrimePromiscuousMatch(Lists.newArrayList(rnaFusion.GeneDown));

        if(fivePrimeProm && threePrimeProm)
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_BOTH_PROM);
        }
        else if(fivePrimeProm)
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_5P_PROM);
        }
        else if(threePrimeProm)
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_3P_PROM);
        }
        else
        {
            rnaFusion.setKnownType(REPORTABLE_TYPE_NONE);
        }
    }

    public void writeRnaMatchData(final String sampleId, final RnaFusionData rnaFusion)
    {
        try
        {
            if(mRnaWriter == null)
            {
                String outputFilename = mOutputDir;

                outputFilename += "SVA_RNA_DATA.csv";

                mRnaWriter = createBufferedWriter(outputFilename, false);

                mRnaWriter.write("SampleId,FusionName,GeneUp,GeneDown,ViableFusion,PhaseMatched,KnownType");

                mRnaWriter.write(",SvIdUp,ChrUp,PosUp,RnaPosUp,OrientUp,StrandUp,TypeUp,ClusterInfoUp");
                mRnaWriter.write(",TransViableUp,TransValidLocUp,TransIdUp,ExonsSkippedUp,RegionTypeUp,CodingTypeUp,ExonUp,DisruptiveUp,DistancePrevUp");

                mRnaWriter.write(",SvIdDown,ChrDown,PosDown,RnaPosDown,OrientDown,StrandDown,TypeDown,ClusterInfoDown");
                mRnaWriter.write(",TransViableDown,TransValidLocDown,TransIdDown,ExonsSkippedDown,RegionTypeDown,CodingTypeDown,ExonDown,DisruptiveDown,DistancePrevDown");

                mRnaWriter.write(",ChainInfo,JunctionReadCount,SpanningFragCount,SpliceType");
                mRnaWriter.write(",ExonMinRankUp,ExonMaxRankUp,ExonMinRankDown,ExonMaxRankDown");

                mRnaWriter.newLine();
            }

            BufferedWriter writer = mRnaWriter;

            writer.write(String.format("%s,%s,%s,%s,%s,%s,%s",
                    sampleId, rnaFusion.Name, rnaFusion.GeneUp, rnaFusion.GeneDown,
                    rnaFusion.isViableFusion(), rnaFusion.isPhaseMatchedFusion(), rnaFusion.getKnownType()));

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
                        transUp.ExonUpstream, transUp.isDisruptive(), transUp.exonDistanceUp()));
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
                                transDown.ExonDownstream, transDown.isDisruptive(), transDown.exonDistanceUp()));
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
            return "IGHJ6";

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
        mPerfCounter.logStats();

        if(mFusionFinder != null)
            mFusionFinder.onCompleted();

        closeBufferedWriter(mRnaWriter);
    }
}
