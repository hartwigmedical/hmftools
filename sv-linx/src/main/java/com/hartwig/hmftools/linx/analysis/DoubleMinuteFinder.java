package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.typeAsInt;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxOutput.SUBSET_DELIM;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getSvTypesStr;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.isOverlapping;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.positionWithin;
import static com.hartwig.hmftools.linx.chaining.ChainFinder.LR_METHOD_DM_CLOSE;
import static com.hartwig.hmftools.linx.types.DoubleMinuteData.INT_SEG_COUNT;
import static com.hartwig.hmftools.linx.types.ResolvedType.DOUBLE_MINUTE;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_DM;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LINK_TYPE_TI;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.linx.chaining.ChainFinder;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.types.DoubleMinuteData;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

public class DoubleMinuteFinder
{
    private CnDataLoader mCnDataLoader;
    private EnsemblDataCache mGeneTransCache;
    private final Map<String, List<SvBreakend>> mChrBreakendMap;
    private final ChainFinder mChainFinder;

    private final List<Integer> mProcessedClusters;
    private final Map<Integer, DoubleMinuteData> mDoubleMinutes;

    private String mOutputDir;
    private BufferedWriter mFileWriter;

    private static final double JCN_THRESHOLD = 8;
    private static final double ADJACENT_JCN_RATIO = 2.3;

    public DoubleMinuteFinder(final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        mChrBreakendMap = chrBreakendMap;
        mChainFinder = new ChainFinder();
        mCnDataLoader = null;
        mGeneTransCache = null;

        mProcessedClusters = Lists.newArrayList();
        mDoubleMinutes = Maps.newHashMap();

        mOutputDir = null;
        mFileWriter = null;
    }

    public void setGeneTransCache(final EnsemblDataCache geneDataCache) { mGeneTransCache = geneDataCache; }
    public void setCopyNumberAnalyser(CnDataLoader cnAnalyser) { mCnDataLoader = cnAnalyser; }
    public void setOutputDir(final String outputDir)
    {
        mOutputDir = outputDir;
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

        final List<SvVarData> candidateDMSVs = Lists.newArrayList();

        for (SvVarData var : cluster.getSVs())
        {
            if (var.jcnMax() < JCN_THRESHOLD)
                continue;

            double svAdjMAPRatio = getAdjacentMajorAPRatio(var);

            if(svAdjMAPRatio >= ADJACENT_JCN_RATIO)
                candidateDMSVs.add(var);
        }

        if(candidateDMSVs.isEmpty())
            return;

        final List<SvChain> dmChains = createDMChains(cluster, candidateDMSVs, false);

        // every SV must be in a chain even if they're not complete, and every chain must either be closed or have a SGL or INF on the end

        if(candidateDMSVs.size() > 2 && cluster.requiresReplication())
        {
            final List<SvChain> variedJcnChains = createDMChains(cluster, candidateDMSVs, true);

            // take the more effective of the 2 approaches
            int uniformChainedCount = countChainedSVs(candidateDMSVs, dmChains);
            int uniformClosedCount = countCloseableChains(dmChains);

            int variedChainedCount = countChainedSVs(candidateDMSVs, variedJcnChains);
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

        int chainedCount = countChainedSVs(candidateDMSVs, dmChains);
        int unchainedCount = candidateDMSVs.size() - chainedCount;

        // single DUPs are closed in a loop, as is any other chain involving all the DM SVs
        for(SvChain dmChain : dmChains)
        {
            if(!dmChain.isClosedLoop() && dmChain.couldCloseChain())
            {
                int maxIndex = dmChain.getLinkedPairs().stream().mapToInt(x -> x.getLinkIndex()).max().orElse(0);
                dmChain.closeChain(LR_METHOD_DM_CLOSE, maxIndex + 1);
                dmChain.setDoubleMinute(true);
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

        LNX_LOGGER.debug("cluster({}) dmSVs({}) chains({}) unchainedSVs({}) {}",
                cluster.id(), candidateDMSVs.size(), dmChains.size(), unchainedCount,
                fullyChained ? "fully chained" : "invalid chain");

        DoubleMinuteData dmData = new DoubleMinuteData(cluster, candidateDMSVs);

        dmData.Chains.addAll(dmChains);
        dmData.FullyChained = fullyChained;

        dmData.setChainCharacteristics(mChrBreakendMap);

        mDoubleMinutes.put(cluster.id(), dmData);

        final List<SvChain> closedChains = dmChains.stream().filter(x -> x.isClosedLoop()).collect(Collectors.toList());
        cluster.setDoubleMinuteData(candidateDMSVs, closedChains); // only the one for now

        cluster.addAnnotation(CLUSTER_ANNOT_DM);

        // cache DUP chains now since the cluster may not go through the chaining routine
        for(SvChain dmChain : dmChains)
        {
            if(dmChain.getSvCount() == 1 && dmChain.getSvList().get(0).type() == DUP)
            {
                final SvVarData dup = dmChain.getSvList().get(0);
                if(!cluster.getChains().stream().anyMatch(x -> x.getSvList().size() == 1 && x.getSvList().contains(dup)))
                {
                    cluster.addChain(dmChain, false);
                }
            }
        }


        if(candidateDMSVs.size() == cluster.getSvCount())
        {
            cluster.setResolved(false, DOUBLE_MINUTE);
        }
    }

    private static double getAdjacentMajorAPRatio(final SvVarData var)
    {
        // get the largest ratio of JCN to the adjacent major AP
        double maxRatio = 0;

        for (int se = SE_START; se <= SE_END; ++se)
        {
            if(se == SE_END && var.isSglBreakend())
                continue;

            final SvBreakend breakend = var.getBreakend(se);

            // gets the CN segment data on the lower side of the breakend (ie opposite to orientation)
            final SvCNData cnData = breakend.getSV().getCopyNumberData(breakend.usesStart(), breakend.orientation() == -1);

            if(cnData == null)
                continue;

            double adjacentMap = cnData.majorAlleleJcn();

            // if against 0, then just ensure it will pass
            maxRatio = adjacentMap > 0 ? max(var.jcn() / adjacentMap, maxRatio) : ADJACENT_JCN_RATIO * 2;
        }

        return maxRatio;
    }

    private final List<SvChain> createDMChains(final SvCluster cluster, final List<SvVarData> dmSvList, boolean applyReplication)
    {
        // first extract stand-alone DUPs to avoid them chaining in just because they can
        final List<SvChain> dmChains = Lists.newArrayList();

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
                    SvChain chain = new SvChain(0);
                    SvLinkedPair pair = new SvLinkedPair(var, var, LINK_TYPE_TI, true, false);
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

        dmChains.addAll(mChainFinder.getUniqueChains());
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

    public void reportCluster(final String sampleId, final SvCluster cluster)
    {
        if(mOutputDir == null || mOutputDir.isEmpty())
            return;

        if(!cluster.hasAnnotation(CLUSTER_ANNOT_DM))
            return;

        final DoubleMinuteData dmData = mDoubleMinutes.get(cluster.id());

        if(dmData == null)
            return;

        // a single DUP or a chain involving all high-JCN SVs which can be made into a loop
        final List<SvVarData> unchainedSVs = Lists.newArrayList();
        final List<SvLinkedPair> linkedPairs = Lists.newArrayList();

        dmData.Chains.forEach(x -> linkedPairs.addAll(x.getLinkedPairs()));

        for(SvVarData var : dmData.SVs)
        {
            if(!dmData.Chains.stream().anyMatch(x -> x.getSvList().contains(var)))
            {
                unchainedSVs.add(var);
            }
        }

        int intExtCount = 0; // count of SV that connect a closed segment to a non closed segment
        double intExtJcnTotal = 0; // Max JCN of SV that connect a closed segment to a non closed segment
        double intExtMaxJcn = 0; // SUM JCN of SV that connect a closed segment to a non closed segment
        double minAdjMAJcnRatio = 0;

        for(SvVarData var : cluster.getSVs())
        {
            if(dmData.SVs.contains(var))
            {
                if(minAdjMAJcnRatio == 0)
                    minAdjMAJcnRatio = getAdjacentMajorAPRatio(var);
                else
                    minAdjMAJcnRatio = min(getAdjacentMajorAPRatio(var), minAdjMAJcnRatio);

                continue;
            }

            if(var.isSglBreakend())
                continue;

            boolean[] breakendInLink = {false, false};

            for(int se = SE_START; se <= SE_END; ++se)
            {
                SvBreakend breakend = var.getBreakend(se);

                breakendInLink[se] = linkedPairs.stream()
                        .anyMatch(x -> x.chromosome().equals(breakend.chromosome())
                                && positionWithin(
                                breakend.position(), x.getBreakend(SE_START).position(), x.getBreakend(SE_END).position()));
            }

            if(breakendInLink[SE_START] != breakendInLink[SE_END])
            {
                ++intExtCount;
                intExtJcnTotal += var.jcn();
                intExtMaxJcn = max(intExtMaxJcn, var.jcn());
            }
        }

        String svIds = "";
        final int[] typeCounts = new int[StructuralVariantType.values().length];
        final List<String> chromosomes = Lists.newArrayList();

        double minDMJcn = 0;
        double maxDMJcn = 0;
        double maxDMCopyNumber = 0;
        double minDMCopyNumber = 0;

        for(final SvVarData var : dmData.SVs)
        {
            ++typeCounts[typeAsInt(var.type())];

            minDMJcn = max(var.jcnMin(), minDMJcn);
            maxDMJcn = max(var.jcn(), maxDMJcn);

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
            if (mFileWriter == null)
            {
                String outputFileName = mOutputDir;

                outputFileName += "LNX_DOUBLE_MINUTES.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,ClusterId,ClusterDesc,ResolvedType,ClusterCount");
                mFileWriter.write(",SamplePurity,SamplePloidy,DMSvCount,DMSvTypes,SvIds,Chromosomes");
                mFileWriter.write(",Chains,FullyChained,ClosedChains,ClosedSegLength,ChainedSVs,Replication");
                mFileWriter.write(",ClosedBreakends,ClosedJcnTotal,OpenBreakends,OpenJcnTotal,OpenJcnMax");
                mFileWriter.write(",IntExtCount,IntExtJcnTotal,IntExtMaxJcn,TotalSegmentCnChange");
                mFileWriter.write(",FbIntCount,FbIntJcnTotal,FbIntJcnMax,SglbIntCount,SglIntJcnTotal,SglIntJcnMax,InfIntCount,InfIntJcnTotal,InfIntJcnMax");
                mFileWriter.write(",MaxCopyNumber,MinJcn,MaxJcn,AmpGenes,CrossCentro,MinAdjMAJcnRatio");
                mFileWriter.newLine();
            }

            mFileWriter.write(String.format("%s,%d,%s,%s,%d",
                    sampleId, cluster.id(), cluster.getDesc(), cluster.getResolvedType(), cluster.getSvCount()));

            mFileWriter.write(String.format(",%.1f,%.1f,%d,%s,%s,%s",
                    samplePurity, samplePloidy, dmData.SVs.size(), dmTypesStr, svIds, chromosomeStr));

            mFileWriter.write(String.format(",%d,%s,%d,%d,%d,%s",
                    dmData.Chains.size(), dmData.FullyChained, dmData.Chains.stream().filter(x -> x.isClosedLoop()).count(),
                    dmData.ClosedSegmentLength, dmData.SVs.size() - unchainedSVs.size(),
                    dmData.Chains.stream().anyMatch(x -> x.hasRepeatedSV())));

            mFileWriter.write(String.format(",%d,%.1f,%d,%.1f,%.1f,%d,%.1f,%.1f,%.1f",
                    dmData.ClosedBreakends, dmData.ClosedJcnTotal, dmData.OpenBreakends, dmData.OpenJcnTotal, dmData.OpenJcnMax,
                    intExtCount, intExtJcnTotal, intExtMaxJcn, dmData.TotalSegmentCnChange));

            mFileWriter.write(String.format(",%s", dmData.internalTypeCountsAsStr()));

            mFileWriter.write(String.format(",%.1f,%.1f,%.1f,%s,%s,%.1f",
                    maxDMCopyNumber, minDMJcn, maxDMJcn, amplifiedGenesStr, dmData.ChainsCentromere, minAdjMAJcnRatio));

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
        for(SvLinkedPair pair : chain.getLinkedPairs())
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
