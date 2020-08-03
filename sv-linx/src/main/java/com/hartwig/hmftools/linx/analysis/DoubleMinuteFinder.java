package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.typeAsInt;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxOutput.SUBSET_DELIM;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getSvTypesStr;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.isOverlapping;
import static com.hartwig.hmftools.linx.chaining.ChainFinder.LR_METHOD_DM_CLOSE;
import static com.hartwig.hmftools.linx.types.LinxConstants.ADJACENT_JCN_RATIO;
import static com.hartwig.hmftools.linx.types.ResolvedType.DOUBLE_MINUTE;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_DM;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.linx.chaining.ChainFinder;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

public class DoubleMinuteFinder
{
    private CnDataLoader mCnDataLoader;
    private EnsemblDataCache mGeneTransCache;
    private final Map<String, List<SvBreakend>> mChrBreakendMap;
    private final ChainFinder mChainFinder;

    private final List<Integer> mProcessedClusters;
    private final Map<Integer, DoubleMinuteData> mDoubleMinutes;

    private BufferedWriter mFileWriter;
    private boolean mLogCandidates;

    private static final double JCN_THRESHOLD = 8;
    private static final double JCN_LOWER_THRESHOLD = 6;
    private static final double LOWER_ADJACENT_JCN_RATIO = 2;
    private static final double MIN_PERC_OF_MAX_JCN = 0.25;

    public DoubleMinuteFinder(final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        mChrBreakendMap = chrBreakendMap;
        mChainFinder = new ChainFinder();
        mCnDataLoader = null;
        mGeneTransCache = null;

        mProcessedClusters = Lists.newArrayList();
        mDoubleMinutes = Maps.newHashMap();

        mFileWriter = null;
    }

    public void setGeneTransCache(final EnsemblDataCache geneDataCache) { mGeneTransCache = geneDataCache; }
    public void setCopyNumberAnalyser(CnDataLoader cnAnalyser) { mCnDataLoader = cnAnalyser; }

    public void setOutputDir(final String outputDir, boolean logCandidates)
    {
        initialiseWriter(outputDir);
        mLogCandidates = logCandidates;
    }

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
        final List<SvVarData> candidateDMSVs = Lists.newArrayList();

        for (SvVarData var : cluster.getSVs())
        {
            if (var.jcn() < JCN_LOWER_THRESHOLD)
                continue;

            double jcn = var.jcn();
            double svAdjMAPRatio = getAdjacentMajorAPRatio(var);

            if(!hasValidSv && jcn >= JCN_THRESHOLD && svAdjMAPRatio >= ADJACENT_JCN_RATIO)
                hasValidSv = true;

            if(jcn >= JCN_LOWER_THRESHOLD && svAdjMAPRatio >= LOWER_ADJACENT_JCN_RATIO)
                candidateDMSVs.add(var);
        }

        if(!hasValidSv)
            return;

        // take any candidate which are at high enough JCN relative to the max for the group
        final double maxDmJcn = candidateDMSVs.stream().mapToDouble(x -> x.jcn()).max().orElse(0);
        final List<SvVarData> dmSVs = Lists.newArrayList();

        candidateDMSVs.stream().filter(x -> x.jcn() >= MIN_PERC_OF_MAX_JCN * maxDmJcn).forEach(x -> dmSVs.add(x));

        // dismiss any single SV which isn't a DUP
        if(dmSVs.size() == 1 && dmSVs.get(0).type() != DUP)
            return;

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

            // all DMs, even if partial or not fully closed will be marked for the visualiser
            dmChain.setDoubleMinute(true);
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

        if(!mLogCandidates && !isDM)
            return;

        LNX_LOGGER.debug("cluster({}) dmSVs({}) chains({}) unchainedSVs({}) {}",
                cluster.id(), dmSVs.size(), dmChains.size(), unchainedCount,
                fullyChained ? "fully chained" : "invalid chain");

        mDoubleMinutes.put(cluster.id(), dmData);

        // mark each chains individually as meeting the DM criteria or not, and subsequently dismiss any SVs not part of those valid chains
        if(isDM)
        {
            // cache DM data against the cluster since it used in the chaining routine amongst other things
            cluster.setDoubleMinuteData(dmSVs, dmChains);

            for(SvChain dmChain : dmChains)
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
            }

            // only resolve clusters of size 1 and 2 as DMs, otherwise just annotate the cluster as containing or being a DM
            if(cluster.getSvCount() <= 2 && cluster.getSvCount() == dmSVs.size())
                cluster.setResolved(false, DOUBLE_MINUTE);

            cluster.addAnnotation(CLUSTER_ANNOT_DM);
        }
    }

    public static double getAdjacentMajorAPRatio(final SvVarData var)
    {
        // get the largest ratio of JCN to the adjacent major AP
        double maxRatio = 0;

        for (int se = SE_START; se <= SE_END; ++se)
        {
            if(se == SE_END && var.isSglBreakend())
                continue;

            final SvBreakend breakend = var.getBreakend(se);

            // gets the major allele JCN on the lower side of the breakend (ie opposite to orientation)
            double adjacentMaJcn = breakend.majorAlleleJcn(breakend.orientation() == -1);

            // if against 0, then just ensure it will pass
            maxRatio = max(var.jcn() / max(adjacentMaJcn, 0.01), maxRatio);
        }

        return maxRatio;
    }

    private final List<SvChain> createDMChains(final SvCluster cluster, final List<SvVarData> dmSvList, boolean applyReplication)
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
                if(dmSvList.stream().filter(x -> x != var).anyMatch(x -> isOverlapping(x, var)))
                {
                    chainSvList.add(var);
                }
                else
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

    private void initialiseWriter(final String outputDir)
    {
        if(outputDir == null || outputDir.isEmpty())
            return;
        
        try
        {
            String outputFileName = outputDir + "LNX_DOUBLE_MINUTES.csv";
            mFileWriter = createBufferedWriter(outputFileName, false);

            mFileWriter.write("SampleId,ClusterId,ClusterDesc,ResolvedType,ClusterCount");
            mFileWriter.write(",SamplePurity,SamplePloidy,IsDM,DMSvCount,DMSvTypes,SvIds,Chromosomes");
            mFileWriter.write(",Chains,FullyChained,ClosedChains,ClosedSegLength,ChainedSVs,Replication");
            mFileWriter.write(",ClosedBreakends,ClosedJcnTotal,OpenBreakends,OpenJcnTotal,OpenJcnMax");
            mFileWriter.write(",NonSegFoldbacks,NonSegFoldbackJcnTotal,SimpleDels");
            mFileWriter.write(",IntExtCount,IntExtJcnTotal,IntExtMaxJcn,FbIntCount,FbIntJcnTotal,FbIntJcnMax");
            mFileWriter.write(",SglbIntCount,SglIntJcnTotal,SglIntJcnMax,InfIntCount,InfIntJcnTotal,InfIntJcnMax");
            mFileWriter.write(",MaxCopyNumber,MinJcn,MaxJcn,AmpGenes,CrossCentro,MinAdjMAJcnRatio");
            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error initialising DM file: {}", e.toString());
        }
    }

    public void reportCluster(final String sampleId, final SvCluster cluster)
    {
        if(mFileWriter == null)
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

        for(final SvVarData var : dmData.SVs)
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

            svIds = appendStr(svIds, var.idStr(), SUBSET_DELIM);
        }

        final String dmTypesStr = getSvTypesStr(typeCounts);

        final PurityContext purityContext = mCnDataLoader.getPurityContext();
        double samplePurity = purityContext != null ? purityContext.bestFit().purity() : 0;
        double samplePloidy = purityContext != null ? purityContext.bestFit().ploidy() : 0;

        // get amplified genes list by looking at all section traversed by this chain or breakends with genes in them?
        String amplifiedGenesStr = "";
        for(SvChain chain : dmData.Chains)
        {
            String amplifiedGenes = getAmplifiedGenesList(chain);
            if(!amplifiedGenes.isEmpty())
                amplifiedGenesStr = appendStr(amplifiedGenesStr, amplifiedGenes, SUBSET_DELIM);
        }

        final String chromosomeStr = appendStrList(chromosomes, SUBSET_DELIM);

        try
        {
            mFileWriter.write(String.format("%s,%d,%s,%s,%d",
                    sampleId, cluster.id(), cluster.getDesc(), cluster.getResolvedType(), cluster.getSvCount()));

            mFileWriter.write(String.format(",%.1f,%.1f,%s,%d,%s,%s,%s",
                    samplePurity, samplePloidy, dmData.isDoubleMinute(), dmData.SVs.size(), dmTypesStr, svIds, chromosomeStr));

            mFileWriter.write(String.format(",%d,%s,%d,%d,%d,%s",
                    dmData.Chains.size(), dmData.FullyChained, dmData.Chains.stream().filter(x -> x.isClosedLoop()).count(),
                    dmData.ClosedSegmentLength, dmData.SVs.size() - dmData.UnchainedSVs.size(),
                    dmData.Chains.stream().anyMatch(x -> x.hasRepeatedSV())));

            mFileWriter.write(String.format(",%d,%.1f,%d,%.1f,%.1f,%d,%.1f,%d",
                    dmData.ClosedBreakends, dmData.ClosedJcnTotal, dmData.OpenBreakends, dmData.OpenJcnTotal, dmData.OpenJcnMax,
                    dmData.NonSegmentFoldbacks, dmData.NonSegmentFoldbackJcnTotal, dmData.SimpleDels));

            mFileWriter.write(String.format(",%s", dmData.internalTypeCountsAsStr()));

            mFileWriter.write(String.format(",%.1f,%.1f,%.1f,%s,%s,%.1f",
                    maxDMCopyNumber, minDMJcn, dmData.MaxJcn, amplifiedGenesStr, dmData.ChainsCentromere, dmData.MinAdjMAJcnRatio));

            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing DM data: {}", e.toString());
        }
    }

    private final String getAmplifiedGenesList(final SvChain chain)
    {
        if(mGeneTransCache == null)
            return "";

        String genesStr = "";
        for(LinkedPair pair : chain.getLinkedPairs())
        {
            String chromosome = pair.chromosome();

            List<EnsemblGeneData> genesList = mGeneTransCache.findGenesByRegion(
                    chromosome, pair.getBreakend(true).position(), pair.getBreakend(false).position());

            if(genesList.isEmpty())
                continue;

            for(final EnsemblGeneData geneData : genesList)
            {
                genesStr = appendStr(genesStr, geneData.GeneName, ';');
            }
        }

        return genesStr;
    }

    public void close()
    {
        closeBufferedWriter(mFileWriter);
    }
}
