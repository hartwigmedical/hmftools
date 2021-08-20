package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.typeAsInt;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxOutput.ITEM_DELIM_CHR;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.DOUBLE_MINUTES;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.runAnnotation;
import static com.hartwig.hmftools.linx.analysis.DoubleMinuteData.getMajorAlleleJcnRatio;
import static com.hartwig.hmftools.linx.analysis.DoubleMinuteData.variantExceedsBothAdjacentJcn;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getSvTypesStr;
import static com.hartwig.hmftools.linx.chaining.ChainFinder.LR_METHOD_DM_CLOSE;
import static com.hartwig.hmftools.linx.types.LinxConstants.ADJACENT_JCN_RATIO;
import static com.hartwig.hmftools.linx.types.ResolvedType.DOUBLE_MINUTE;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_DM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.linx.CohortDataWriter;
import com.hartwig.hmftools.linx.CohortFileInterface;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.chaining.ChainFinder;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class DoubleMinuteFinder implements CohortFileInterface
{
    private CnDataLoader mCnDataLoader;
    private EnsemblDataCache mGeneTransCache;
    private final Map<String, List<SvBreakend>> mChrBreakendMap;
    private final ChainFinder mChainFinder;

    private final List<Integer> mProcessedClusters;
    private final Map<Integer, DoubleMinuteData> mDoubleMinutes;

    private final CohortDataWriter mCohortDataWriter;
    private final BufferedWriter mSampleWriter;
    private final boolean mLogCandidates;
    private final boolean mShowCandidates;
    private final LinxConfig mConfig;

    protected static final double JCN_UPPER_THRESHOLD = 8;
    private static final double JCN_THRESHOLD = 5;
    private static final double LOWER_ADJACENT_JCN_RATIO = 2;
    private static final double MIN_PERC_OF_MAX_JCN = 0.25;
    protected static final int MIN_SEGMENT_DEPTH_WINDOW_COUNT = 6;

    private static final String COHORT_WRITER_DM = "DoubleMinute";

    public DoubleMinuteFinder(
            final LinxConfig config, final CohortDataWriter cohortDataWriter, final Map<String,List<SvBreakend>> chrBreakendMap)
    {
        mChrBreakendMap = chrBreakendMap;
        mChainFinder = new ChainFinder();
        mCnDataLoader = null;
        mGeneTransCache = null;

        mProcessedClusters = Lists.newArrayList();
        mDoubleMinutes = Maps.newHashMap();
        mConfig = config;

        mLogCandidates = runAnnotation(config.RequiredAnnotations, DOUBLE_MINUTES);
        mShowCandidates = runAnnotation(config.RequiredAnnotations, "SHOW_DM");

        mCohortDataWriter = cohortDataWriter;

        if(mLogCandidates && mConfig.isSingleSample())
        {
            mSampleWriter = initialiseWriter(mConfig.OutputDataPath, config.getSampleIds().get(0));
        }
        else
        {
            mSampleWriter = null;
        }
    }

    public void setGeneTransCache(final EnsemblDataCache geneDataCache) { mGeneTransCache = geneDataCache; }
    public void setCopyNumberAnalyser(CnDataLoader cnAnalyser) { mCnDataLoader = cnAnalyser; }

    public final Map<Integer, DoubleMinuteData> getDoubleMinutes() { return mDoubleMinutes; }

    public void clear()
    {
        mProcessedClusters.clear();
        mDoubleMinutes.clear();
    }

    public void analyseCluster(SvCluster cluster)
    {
        analyseCluster(cluster, false);
    }

    public void analyseCluster(SvCluster cluster, boolean reassess)
    {
        // don't re-analyse and cluster unless it has changed constituent SVs
        if(mProcessedClusters.contains(cluster.id()))
        {
            if(!reassess)
                return;

            cluster.clearDoubleMinuteData();
        }
        else
        {
            mProcessedClusters.add(cluster.id());
        }

        boolean hasValidSv = false;
        final List<SvVarData> candidateDmSVs = Lists.newArrayList();
        final List<SvVarData> candidateFlankedSVs = Lists.newArrayList();

        for (SvVarData var : cluster.getSVs())
        {
            if (var.jcn() < JCN_THRESHOLD)
                continue;

            double svAdjMAPRatio = getAdjacentMajorAPRatio(var);

            if(svAdjMAPRatio >= ADJACENT_JCN_RATIO)
                hasValidSv = true;

            if(svAdjMAPRatio >= LOWER_ADJACENT_JCN_RATIO)
                candidateDmSVs.add(var);

            candidateFlankedSVs.add(var);
        }

        if(!hasValidSv)
            return;

        // take any candidate which are at high enough JCN relative to the max for the group
        double maxDmJcn = candidateDmSVs.stream().mapToDouble(x -> x.jcn()).max().orElse(0);
        final double minJcn = max(JCN_THRESHOLD, MIN_PERC_OF_MAX_JCN * maxDmJcn);
        final List<SvVarData> dmSVs = Lists.newArrayList();

        candidateDmSVs.stream().filter(x -> x.jcn() >= minJcn).forEach(x -> dmSVs.add(x));

        // dismiss any single SV which isn't a DUP
        if(dmSVs.size() == 1)
        {
            final SvVarData var = dmSVs.get(0);

            if(var.type() != DUP)
                return;

            if(!variantExceedsBothAdjacentJcn(var))
                return;
        }

        final List<SvChain> dmChains = createDMChains(cluster, dmSVs, false);

        if(dmSVs.size() > 2 && cluster.requiresReplication())
        {
            final List<SvChain> variedJcnChains = createDMChains(cluster, dmSVs, true);

            // take the more effective of the 2 approaches
            int uniformChainedCount = countChainedSVs(dmSVs, dmChains);
            int uniformClosedCount = countCloseableChains(dmChains);

            int variedChainedCount = countChainedSVs(dmSVs, variedJcnChains);
            int variedClosedCount = countCloseableChains(variedJcnChains);

            if(variedClosedCount == variedJcnChains.size() && uniformClosedCount < dmChains.size()
            && variedChainedCount >= uniformChainedCount)
            {
                LNX_LOGGER.debug("cluster({}) varied JCN chains({} closed={} chained={}) better than uniform({} closed={} chained={}",
                        cluster.id(), variedJcnChains.size(), variedClosedCount, variedChainedCount,
                        dmChains.size(), uniformClosedCount, uniformChainedCount);

                dmChains.clear();
                dmChains.addAll(variedJcnChains);
            }
        }

        int chainedCount = countChainedSVs(dmSVs, dmChains);
        int unchainedCount = dmSVs.size() - chainedCount;

        // single DUPs are closed in a loop, as is any other chain involving all the DM SVs
        for(SvChain dmChain : dmChains)
        {
            if(!dmChain.isClosedLoop() && dmChain.couldCloseChain())
            {
                int maxIndex = dmChain.getLinkedPairs().stream().mapToInt(x -> x.getLinkIndex()).max().orElse(0);
                dmChain.closeChain(LR_METHOD_DM_CLOSE, maxIndex + 1);
            }
        }

        boolean fullyChained = !dmChains.isEmpty();

        if(unchainedCount > 0)
        {
            fullyChained = false;
        }

        for(SvChain chain : dmChains)
        {
            if(chain.isClosedLoop())
                continue;

            if(!chain.getChainEndSV(true).isSglBreakend() && !chain.getChainEndSV(false).isSglBreakend())
            {
                fullyChained = false;
                break;
            }
        }

        DoubleMinuteData dmData = new DoubleMinuteData(cluster, dmSVs);

        dmData.Chains.addAll(dmChains);
        dmData.FullyChained = fullyChained;

        dmData.annotate(mChrBreakendMap);

        boolean isDM = dmData.isDoubleMinute();

        if(!mLogCandidates && !mShowCandidates && !isDM)
            return;

        LNX_LOGGER.debug("cluster({}) dmSVs({}) chains({}) unchainedSVs({}) {}",
                cluster.id(), dmSVs.size(), dmChains.size(), unchainedCount,
                fullyChained ? "fully chained" : "invalid chain");

        mDoubleMinutes.put(cluster.id(), dmData);

        // mark each chains individually as meeting the DM criteria or not, and subsequently dismiss any SVs not part of those valid chains
        if(isDM)
        {
            // cache DM data against the cluster since it used in the chaining routine amongst other things
            cluster.setDoubleMinuteData(dmData.ValidSVs, dmData.ValidChains);

            for(SvChain dmChain : dmData.ValidChains)
            {
                // cache single DUP chains now since the cluster may not go through the chaining routine
                if(dmChain.getSvCount() == 1 && dmChain.getSvList().get(0).type() == DUP)
                {
                    final SvVarData dup = dmChain.getSvList().get(0);
                    if(!cluster.getChains().stream().anyMatch(x -> x.getSvList().size() == 1 && x.getSvList().contains(dup)))
                    {
                        cluster.addChain(dmChain, false);
                    }
                }

                // all DMs, even if partial or not fully closed will be marked for the visualiser
                dmChain.setDoubleMinute(true);
            }

            // only resolve clusters of size 1 and 2 as DMs, otherwise just annotate the cluster as containing or being a DM
            if(cluster.getSvCount() <= 2 && cluster.getSvCount() == dmSVs.size())
                cluster.setResolved(false, DOUBLE_MINUTE);

            cluster.addAnnotation(CLUSTER_ANNOT_DM);
        }
        else if(mShowCandidates)
        {
            cluster.setDoubleMinuteData(dmData.SVs, dmData.Chains);

            for(SvChain dmChain : dmData.Chains)
            {
                // cache single DUP chains now since the cluster may not go through the chaining routine
                if(dmChain.getSvCount() == 1 && dmChain.getSvList().get(0).type() == DUP)
                {
                    final SvVarData dup = dmChain.getSvList().get(0);
                    if(!cluster.getChains().stream().anyMatch(x -> x.getSvList().size() == 1 && x.getSvList().contains(dup)))
                    {
                        cluster.addChain(dmChain, false);
                    }
                }

                // all DMs, even if partial or not fully closed will be marked for the visualiser
                dmChain.setDoubleMinute(true);
            }
        }
    }

    protected static double getAdjacentMajorAPRatio(final SvVarData var)
    {
        // get the largest ratio of JCN to the adjacent major AP
        double maxRatio = getMajorAlleleJcnRatio(var.getBreakend(true));

        if(!var.isSglBreakend())
            maxRatio = max(getMajorAlleleJcnRatio(var.getBreakend(false)), maxRatio);

        return maxRatio;
    }

    private List<SvChain> createDMChains(final SvCluster cluster, final List<SvVarData> dmSvList, boolean applyReplication)
    {
        // first extract stand-alone DUPs to avoid them chaining in just because they can
        final List<SvChain> dmChains = Lists.newArrayList();
        int chainId = 0;

        final List<SvVarData> chainSvList = Lists.newArrayList();

        for(SvVarData var : dmSvList)
        {
            if(var.type() != DUP)
            {
                chainSvList.add(var);
            }
            else
            {
                boolean hasOverlap = dmSvList.stream().filter(x -> x != var)
                        .anyMatch(x -> positionWithin(x.position(true), var.position(true), var.position(false))
                                || positionWithin(x.position(false), var.position(true), var.position(false)));

                if(!hasOverlap && variantExceedsBothAdjacentJcn(var))
                {
                    // special case creating a chain out of a DUP
                    SvChain chain = new SvChain(chainId++);
                    LinkedPair pair = LinkedPair.from(var, var, true, false);
                    chain.addLink(pair, true);
                    chain.setJcnData(var.jcn(), var.jcnUncertainty());
                    chain.closeChain();
                    chain.setDoubleMinute(true);
                    dmChains.add(chain);
                }
                else
                {
                    chainSvList.add(var);
                }
            }
        }

        mChainFinder.initialise(cluster, chainSvList, applyReplication);
        mChainFinder.formChains(false);

        for(SvChain chain : mChainFinder.getUniqueChains())
        {
            chain.setId(chainId++);
            dmChains.add(chain);
        }

        mChainFinder.clear();
        return dmChains;
    }

    private static int countChainedSVs(final List<SvVarData> svList, final List<SvChain> chains)
    {
        return (int)svList.stream().filter(x -> chains.stream().anyMatch(y -> y.getSvList().contains(x))).count();
    }

    private static int countCloseableChains(final List<SvChain> chains)
    {
        return (int)chains.stream().filter(x -> x.couldCloseChain()).count();
    }

    @Override
    public String fileType() { return COHORT_WRITER_DM; }

    @Override
    public BufferedWriter createWriter(final String outputDir)
    {
        return initialiseWriter(outputDir, "");
    }

    private static BufferedWriter initialiseWriter(final String outputDir, final String sampleId)
    {
        if(outputDir == null || outputDir.isEmpty())
            return null;

        try
        {
            String outputFileName = !sampleId.isEmpty() ?
                    outputDir + sampleId + ".linx.ecdna.csv" : outputDir + "LNX_ECDNA.csv";

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            if(sampleId.isEmpty())
                writer.write("SampleId,");

            writer.write("ClusterId,ClusterDesc,ResolvedType,ClusterCount");
            writer.write(",SamplePurity,SamplePloidy,IsDM,DMSvCount,DMSvTypes,SvIds,Chromosomes");
            writer.write(",Chains,FullyChained,ClosedChains,ClosedSegLength,ChainedSVs,Replication");
            writer.write(",ClosedBreakends,ClosedJcnTotal,OpenBreakends,OpenJcnTotal,OpenJcnMax");
            writer.write(",NonSegFoldbacks,NonSegFoldbackJcnTotal,SimpleDels");
            writer.write(",IntExtCount,IntExtJcnTotal,IntExtMaxJcn,FbIntCount,FbIntJcnTotal,FbIntJcnMax");
            writer.write(",SglbIntCount,SglIntJcnTotal,SglIntJcnMax,InfIntCount,InfIntJcnTotal,InfIntJcnMax");
            writer.write(",MaxCopyNumber,MinJcn,MaxJcn,AmpGenes,CrossCentro,MinOfMaxAdjMajRatio,MinOfMinAdjMajRatio");
            writer.newLine();
            return writer;
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error initialising DM file: {}", e.toString());
            return null;
        }
    }

    public void reportCluster(final String sampleId, final SvCluster cluster)
    {
        if(!mLogCandidates)
            return;

        final DoubleMinuteData dmData = mDoubleMinutes.get(cluster.id());

        if(dmData == null)
            return;

        String svIds = "";
        final int[] typeCounts = new int[StructuralVariantType.values().length];
        final List<String> chromosomes = Lists.newArrayList();

        double minDMJcn = 0;
        double maxDMCopyNumber = 0;
        double minDMCopyNumber = 0;

        for(final SvVarData var : dmData.ValidSVs)
        {
            ++typeCounts[typeAsInt(var.type())];

            if(minDMJcn == 0)
                minDMJcn = var.jcn();
            else
                minDMJcn = min(minDMJcn, var.jcn());

            maxDMCopyNumber = max(maxDMCopyNumber, max(var.copyNumber(true), var.copyNumber(false)));

            double minCopyNumber = var.isSglBreakend() ? var.copyNumber(true) : min(var.copyNumber(true), var.copyNumber(false));
            if(minDMCopyNumber == 0 || minCopyNumber < minDMCopyNumber)
                minDMCopyNumber = minCopyNumber;

            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(var.isSglBreakend() && se== SE_END)
                    continue;

                final SvBreakend breakend = var.getBreakend(se);
                final String chromosome = breakend.chromosome();

                if(!chromosomes.contains(chromosome))
                    chromosomes.add(chromosome);
            }

            svIds = appendStr(svIds, var.idStr(), ITEM_DELIM_CHR);
        }

        final String dmTypesStr = getSvTypesStr(typeCounts);

        final PurityContext purityContext = mCnDataLoader.getPurityContext();
        double samplePurity = purityContext != null ? purityContext.bestFit().purity() : 0;
        double samplePloidy = purityContext != null ? purityContext.bestFit().ploidy() : 0;

        // get amplified genes list by looking at all section traversed by this chain or breakends with genes in them?
        String amplifiedGenesStr = "";
        for(SvChain chain : dmData.ValidChains)
        {
            String amplifiedGenes = getAmplifiedGenesList(chain);
            if(!amplifiedGenes.isEmpty())
                amplifiedGenesStr = appendStr(amplifiedGenesStr, amplifiedGenes, ITEM_DELIM_CHR);
        }

        final String chromosomeStr = appendStrList(chromosomes, ITEM_DELIM_CHR);

        BufferedWriter writer = null;

        if(mConfig.isSingleSample())
        {
            writer = mSampleWriter;
        }
        else
        {
            writer = mCohortDataWriter.getWriter(this);
        }

        try
        {
            if(!mConfig.isSingleSample())
            {
                writer.write(String.format("%s,", sampleId));
            }

            writer.write(String.format("%d,%s,%s,%d",
                    cluster.id(), cluster.getDesc(), cluster.getResolvedType(), cluster.getSvCount()));

            writer.write(String.format(",%.1f,%.1f,%s,%d,%s,%s,%s",
                    samplePurity, samplePloidy, dmData.isDoubleMinute(), dmData.ValidSVs.size(), dmTypesStr, svIds, chromosomeStr));

            writer.write(String.format(",%d,%s,%d,%d,%d,%s",
                    dmData.ValidChains.size(), dmData.FullyChained, dmData.ValidChains.stream().filter(x -> x.isClosedLoop()).count(),
                    dmData.ClosedSegmentLength, dmData.SVs.size() - dmData.UnchainedSVs.size(),
                    dmData.ValidChains.stream().anyMatch(x -> x.hasRepeatedSV())));

            writer.write(String.format(",%d,%.1f,%d,%.1f,%.1f,%d,%.1f,%d",
                    dmData.ClosedBreakends, dmData.ClosedJcnTotal, dmData.OpenBreakends, dmData.OpenJcnTotal, dmData.OpenJcnMax,
                    dmData.NonSegmentFoldbacks, dmData.NonSegmentFoldbackJcnTotal, dmData.SimpleDels));

            writer.write(String.format(",%s", dmData.internalTypeCountsAsStr()));

            double minMinAmr = 0;
            if(dmData.SVs.size() == 1)
            {
                final SvVarData var = dmData.SVs.get(0);

                minMinAmr = min(
                        getMajorAlleleJcnRatio(var.getBreakend(true)),
                        getMajorAlleleJcnRatio(var.getBreakend(false)));
            }

            writer.write(String.format(",%.1f,%.1f,%.1f,%s,%s,%.1f,%.1f",
                    maxDMCopyNumber, minDMJcn, dmData.MaxJcn, amplifiedGenesStr, dmData.ChainsCentromere,
                    dmData.MinAdjMAJcnRatio, minMinAmr));

            writer.newLine();
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing DM data: {}", e.toString());
        }
    }

    private String getAmplifiedGenesList(final SvChain chain)
    {
        if(mGeneTransCache == null)
            return "";

        String genesStr = "";
        for(LinkedPair pair : chain.getLinkedPairs())
        {
            String chromosome = pair.chromosome();

            List<GeneData> genesList = mGeneTransCache.findGenesByRegion(
                    chromosome, pair.getBreakend(true).position(), pair.getBreakend(false).position());

            if(genesList.isEmpty())
                continue;

            for(final GeneData geneData : genesList)
            {
                genesStr = appendStr(genesStr, geneData.GeneName, ';');
            }
        }

        return genesStr;
    }

    public void close()
    {
        closeBufferedWriter(mSampleWriter);
    }
}
