package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.typeAsInt;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getSvTypesStr;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNONTATION_DM;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LINK_TYPE_TI;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isSpecificSV;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvChain;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class DoubleMinuteFinder
{
    private CnDataLoader mCnAnalyser;
    private SvGeneTranscriptCollection mGeneTransCache;
    private String mOutputDir;

    private static double DM_PLOIDY_MIN_RATIO = 2.3;
    private static double DM_MIN_PLOIDY = 3;

    private static double DM_PLOIDY_INCOMPLETE_MIN_RATIO = 4;
    private static double DM_INCOMPLETE_MIN_PLOIDY = 10;

    private static int DM_MAX_SV_COUNT = 16;

    private static final Logger LOGGER = LogManager.getLogger(DoubleMinuteFinder.class);

    public DoubleMinuteFinder()
    {
        mCnAnalyser = null;
        mGeneTransCache = null;
        mOutputDir = null;
    }

    public void setGeneTransCache(final SvGeneTranscriptCollection geneTransCache) { mGeneTransCache = geneTransCache; }
    public void setCopyNumberAnalyser(CnDataLoader cnAnalyser) { mCnAnalyser = cnAnalyser; }
    public void setOuptutDir(final String outputDir) { mOutputDir = outputDir; }

    public void findPotentialDoubleMinuteClusters(final String sampleId, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        /* Identify potential DM clusters if:
            - each of their breakends is X times the major allele ploidy of the near CN segment
            - they have offsetting facing breakends on each arm
            - any SVs they skip over have ploidy X times less than the DM
            - skip any groups containing SGLs/NONEs

           In addition, check any loose SVs which meet the above criteria but weren't complete:
           - eg if they contain SGLs or NONES
           - they aren't all offset on an arm
        */

        List<SvBreakend> incompleteDMBreakends = Lists.newArrayList();

        final Map<String,List<SvCNData>> chrCopyNumberDataMap = mCnAnalyser.getChrCnDataMap();

        if(chrCopyNumberDataMap.isEmpty())
            return;

        for (final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            final String chromosome = entry.getKey();
            List<SvBreakend> breakendList = entry.getValue();

            double prevPloidy = 0;
            double minDMPloidy = 0;
            double minDMCopyNumber = 0;
            double maxDMCopyNumber = 0;
            boolean inPotentialDM = false;
            boolean isDmGroupResolved = false;

            final List<SvCNData> cnDataList = chrCopyNumberDataMap.get(chromosome);

            if(cnDataList == null || cnDataList.isEmpty())
            {
                LOGGER.warn("sample({}) missing CN data for chromosome({})", sampleId, chromosome);
                continue;
            }

            double telomereMAP = cnDataList.get(0).majorAllelePloidy();

            List<SvVarData> dmSVList = Lists.newArrayList();
            List<SvBreakend> dmBreakendList = Lists.newArrayList();
            int overlappedCount = 0;

            for(final SvBreakend breakend : breakendList)
            {
                double minPloidy = breakend.getSV().ploidyMin();
                double maxPloidy = breakend.getSV().ploidyMax();

                isSpecificSV(breakend.getSV());

                // first check the min ploidy vs the adjacent CN segment
                double adjacentMap = getAdjacentMajorAllelePloidy(breakend);

                if(Double.isNaN(adjacentMap))
                    continue;

                if(minPloidy >= DM_INCOMPLETE_MIN_PLOIDY && minPloidy >= DM_PLOIDY_INCOMPLETE_MIN_RATIO * adjacentMap)
                {
                    // keep track of all high-ploidy breakends for later analysis of incomplete DM groups
                    incompleteDMBreakends.add(breakend);
                }

                if(!inPotentialDM)
                {
                    if(breakend.orientation() != -1)
                    {
                        prevPloidy = maxPloidy;
                        continue;
                    }

                    if(minPloidy < DM_MIN_PLOIDY || minPloidy < telomereMAP * DM_PLOIDY_MIN_RATIO
                            || (prevPloidy > 0 && minPloidy < prevPloidy * DM_PLOIDY_MIN_RATIO))
                    {
                        prevPloidy = maxPloidy;
                        continue;
                    }

                    if(!isValidDMBreakend(breakend))
                        continue;

                    final SvVarData var = breakend.getSV();

                    // check the major allele ploidy outside this breakend
                    final SvCNData prevCnData = var.getCopyNumberData(breakend.usesStart(), true);
                    double prevMap = prevCnData.majorAllelePloidy();

                    if(minPloidy < prevMap * DM_PLOIDY_MIN_RATIO)
                        continue;

                    // satisfies the conditions to start a potential DM
                    minDMPloidy = minPloidy;
                    minDMCopyNumber = breakend.copyNumber();
                    maxDMCopyNumber = minDMCopyNumber;

                    inPotentialDM = true;

                    dmBreakendList.clear();
                    dmSVList.clear();
                    overlappedCount = 0;

                    dmBreakendList.add(breakend);
                    dmSVList.add(var);

                    isDmGroupResolved = false;
                }
                else
                {
                    // skip any low ploidy SVs
                    if(maxPloidy * DM_PLOIDY_MIN_RATIO < minDMPloidy)
                    {
                        ++overlappedCount;
                        continue;
                    }

                    final SvVarData var = breakend.getSV();

                    // cancel a potential group if it encounters a high-ploidy SGL
                    if(!isValidDMBreakend(breakend))
                    {
                        inPotentialDM = false;
                        LOGGER.debug("cancelling potential DM group(count={} first={}) at invalid SV", dmSVList.size(), dmSVList.get(0).posId());
                        continue;
                    }

                    minDMPloidy = min(minPloidy, minDMPloidy);
                    minDMCopyNumber = min(breakend.copyNumber(), minDMCopyNumber);
                    maxDMCopyNumber = max(breakend.copyNumber(), maxDMCopyNumber);

                    dmBreakendList.add(breakend);

                    if(!dmSVList.contains(var))
                        dmSVList.add(var);

                    if(breakend.orientation() == 1)
                    {
                        isDmGroupResolved = isResolvedDMGroup(dmBreakendList, dmSVList);

                        if (isDmGroupResolved)
                        {
                            // check that the next breakend drops back below the required threshold
                            if(minDMPloidy < adjacentMap * DM_PLOIDY_MIN_RATIO)
                            {
                                LOGGER.debug(String.format("cancelling potential DM group(count=%d first=%s): minDmPloidy(%.2f) vs nextMap(%.2f)",
                                        dmSVList.size(), dmSVList.get(0).posId(), minDMPloidy, adjacentMap));

                                inPotentialDM = false;
                                continue;
                            }

                            LOGGER.debug("DM group identified & complete");

                            // remove from consideration of incomplete DM groups
                            dmBreakendList.stream().forEach(x -> incompleteDMBreakends.remove(x));

                            reportPotentialGroup(sampleId, dmSVList, true, overlappedCount);
                            inPotentialDM = false;
                            continue;
                        }
                    }

                    if(dmSVList.size() >= DM_MAX_SV_COUNT && !isDmGroupResolved)
                    {
                        LOGGER.debug("cancelling potential DM group: count({}) firstCluster({}) firstSV({}) at max SV count",
                                dmSVList.size(), dmSVList.get(0).getCluster().id(), dmSVList.get(0).posId());
                        inPotentialDM = false;
                        continue;
                    }
                }
            }
        }

        if(incompleteDMBreakends.isEmpty())
            return;

        // now follow-up on any incomplete potential DM groups

        // first check that both breakends are present for all non-SGL SVs
        List<SvVarData> incompleteDMSvList = Lists.newArrayList();

        for(int i = 0; i < incompleteDMBreakends.size(); ++i)
        {
            SvBreakend breakend = incompleteDMBreakends.get(i);

            if(breakend.getSV().isNullBreakend())
            {
                incompleteDMSvList.add(breakend.getSV());
            }
            else
            {
                for (int j = i+1; j < incompleteDMBreakends.size(); ++j)
                {
                    SvBreakend otherBreakend = incompleteDMBreakends.get(j);
                    if (breakend.getSV() == otherBreakend.getSV())
                    {
                        incompleteDMSvList.add(breakend.getSV());
                        break;
                    }
                }
            }
        }

        if(incompleteDMSvList.isEmpty())
            return;

        // follow a chain through any linked chromosomes to create a set of related SVs
        while(!incompleteDMSvList.isEmpty())
        {
            List<String> chromosomeList = Lists.newArrayList();

            SvVarData var = incompleteDMSvList.get(0);

            if(!chromosomeList.contains(var.chromosome(true)))
                chromosomeList.add(var.chromosome(true));

            if(!var.isNullBreakend() && !chromosomeList.contains(var.chromosome(false)))
                chromosomeList.add(var.chromosome(false));

            List<SvVarData> dmSvList = Lists.newArrayList();
            dmSvList.add(var);

            boolean addedMoreChromosomes = true;

            while(addedMoreChromosomes)
            {
                addedMoreChromosomes = false;

                for (int j = 1; j < incompleteDMSvList.size(); ++j)
                {
                    SvVarData otherVar = incompleteDMSvList.get(j);

                    if (!chromosomeList.contains(otherVar.chromosome(true)) && !chromosomeList.contains(otherVar.chromosome(false)))
                        continue;

                    if(!dmSvList.contains(otherVar))
                        dmSvList.add(otherVar);

                    if(!chromosomeList.contains(otherVar.chromosome(true)))
                    {
                        addedMoreChromosomes = true;
                        chromosomeList.add(otherVar.chromosome(true));
                        break;
                    }

                    if(!var.isNullBreakend() && !chromosomeList.contains(otherVar.chromosome(false)))
                    {
                        addedMoreChromosomes = true;
                        chromosomeList.add(otherVar.chromosome(false));
                        break;
                    }
                }
            }

            if(dmSvList.isEmpty())
                break;

            dmSvList.stream().forEach(x -> incompleteDMSvList.remove(x));

            reportPotentialGroup(sampleId, dmSvList, false, 0);
        }
    }

    private static double getAdjacentMajorAllelePloidy(final SvBreakend breakend)
    {
        // gets the CN segment data on the lower side of the breakend (ie opposite to orientation)
        final SvCNData cnData = breakend.getSV().getCopyNumberData(breakend.usesStart(), breakend.orientation() == -1);

        if(cnData == null)
            return Double.NaN;

        return cnData.majorAllelePloidy();
    }

    private static boolean isValidDMBreakend(final SvBreakend breakend)
    {
        final SvVarData var = breakend.getSV();

        if (var.isNullBreakend())
            return false;

        if (var.isLocal())
            return true;

        // the other end must be a in a TI
        final SvLinkedPair remoteTI = var.getLinkedPair(!breakend.usesStart());

        if (remoteTI == null)
            return false;

        // check ploidy context of remote breakends as well
        for(int be = SE_START; be <= SE_END; ++be)
        {
            boolean isStart = isStart(be);
            final SvBreakend remoteBreakend = remoteTI.getBreakend(isStart);
            final SvVarData remoteSV = remoteBreakend.getSV();

            // the other end of this remote TI doesn't link back to the original chromosome
            final SvBreakend remoteOtherBreakend = remoteSV.getBreakend(!remoteBreakend.usesStart());

            if(remoteOtherBreakend == null)
                return false;

            if(!remoteOtherBreakend.getChrArm().equals(breakend.getChrArm()))
                return false;

            // check that the next breakend drops back below the required threshold
            final SvCNData applicableCNData = remoteBreakend.getSV().getCopyNumberData(remoteBreakend.usesStart(), isStart);

            if (applicableCNData == null)
            {
                LOGGER.warn("missing CN data for DM SV({})", remoteSV.id());
                return false;
            }

            double outsideMap = applicableCNData.majorAllelePloidy();

            if (var.ploidyMin() < outsideMap * DM_PLOIDY_MIN_RATIO)
                return false;
        }

        return true;
    }

    private static boolean isResolvedDMGroup(List<SvBreakend> breakendList, List<SvVarData> svList)
    {
        // check that every local SVs has recorded both its breakends and every remote SV
        // forms a TI with another remote TI
        int orientationTotal = 0;

        for(final SvVarData var : svList)
        {
            final SvBreakend startBreakend = var.getBreakend(true);
            final SvBreakend endBreakend = var.getBreakend(false);

            orientationTotal += calcConsistency(var);

            if(startBreakend.getChrArm().equals(endBreakend.getChrArm()))
            {
                if(!breakendList.contains(startBreakend) || !breakendList.contains(endBreakend))
                    return false;
            }
            else
            {
                // does this SV form a remote TI with another SV in this group?
                final SvBreakend remoteBreakend = breakendList.contains(startBreakend) ? endBreakend : startBreakend;
                final SvLinkedPair remoteTI = var.getLinkedPair(remoteBreakend.usesStart());

                if(remoteTI == null)
                    return false;

                final SvBreakend otherRemoteBreakend = remoteTI.getOtherBreakend(remoteBreakend);

                if(!svList.contains(otherRemoteBreakend.getSV()))
                    return false;
            }
        }

        if(orientationTotal != 0)
            return false;

        return true;
    }

    private void reportPotentialGroup(final String sampleId, List<SvVarData> dmSVList, boolean isComplete,
            int overlappedCount)
    {
        String clusterInfo = "";
        String svInfo = "";
        List<SvCluster> clusters = Lists.newArrayList();

        int[] typeCounts = new int[StructuralVariantType.values().length];
        List<String> chromosomes = Lists.newArrayList();

        double minDMPloidy = 0;
        double maxDMCopyNumber = 0;

        for(final SvVarData var : dmSVList)
        {
            ++typeCounts[typeAsInt(var.type())];

            if(minDMPloidy == 0 || var.ploidyMin() < minDMPloidy)
                minDMPloidy = var.ploidyMin();

            maxDMCopyNumber = max(maxDMCopyNumber, max(var.copyNumber(true), var.copyNumber(false)));

            if(!chromosomes.contains(var.chromosome(true)))
                chromosomes.add(var.chromosome(true));

            if(!var.isNullBreakend() && !chromosomes.contains(var.chromosome(false)))
                chromosomes.add(var.chromosome(false));

            svInfo = appendStr(svInfo, var.id(), ';');

            final SvCluster cluster = var.getCluster();
            if (!clusters.contains(cluster))
            {
                clusters.add(cluster);

                cluster.addAnnotation(CLUSTER_ANNONTATION_DM);
                clusterInfo = appendStr(clusterInfo, String.format("%d-%s-%s", cluster.id(), cluster.getDesc(), cluster.getResolvedType()), ';');
            }
        }

        final String dmTypesStr = getSvTypesStr(typeCounts);

        final PurityContext purityContext = mCnAnalyser.getPurityContext();
        double samplePurity = purityContext != null ? purityContext.bestFit().purity() : 0;
        double samplePloidy = purityContext != null ? purityContext.bestFit().ploidy() : 0;

        final SvChain chain = isComplete ? createDMChain(dmSVList) : null;

        long dmChainLength = chain != null ? chain.getLength(true) : 0;

        // get amplified genes list by looking at all section tranversed by this chain or breakends with genes in them?
        String amplifiedGenesStr = chain != null ? getAmplifiedGenesList(chain) : "";

        long posStart = 0;
        long posEnd = 0;

        if(isComplete && dmSVList.size() == 1)
        {
            SvVarData dup = dmSVList.get(0);
            posStart = dup.position(true);
            posEnd = dup.position(false);
        }

        String chromosomeStr = "";
        for(String chr : chromosomes)
        {
            chromosomeStr = appendStr(chromosomeStr, chr, ';');
        }

        // SampleId,SamplePurity,SamplePloidy,IsComplete,GroupCount,ClusterInfo,SvTypes,SvInfo,Chromosomes,DMPosStart,DMPosEnd,
        // MaxDMCopyNumber,MinDMPloidy,SVOverlapCount,DMLength,AmplifiedGenes
        String infoStr = String.format("%s,%.2f,%.2f,%s,%d,%s,%s,%s,%s,%d,%d",
                sampleId, samplePurity, samplePloidy, isComplete, dmSVList.size(), clusterInfo, dmTypesStr, svInfo,
                chromosomeStr, posStart, posEnd);

        infoStr += String.format(",%.2f,%.2f,%d,%d,%s",
                maxDMCopyNumber, minDMPloidy, overlappedCount, dmChainLength, amplifiedGenesStr);

        LOGGER.info("POTENTIAL_DM_DATA: {}", infoStr);
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

    private static final SvChain createDMChain(List<SvVarData> dmSVList)
    {
        if(dmSVList.size() == 1)
        {
            // special case creating a chain out of a DUP
            SvChain chain = new SvChain(0);
            final SvVarData var = dmSVList.get(0);
            if(var.type() != DUP)
                return null;

            SvLinkedPair pair = new SvLinkedPair(var, var, LINK_TYPE_TI, true, false);
            chain.addLink(pair, true);
            return chain;
        }

        // create a temporary cluster and try to chain it
        SvCluster dmCluster = new SvCluster(0);

        for(SvVarData var : dmSVList)
        {
            SvVarData copySV = new SvVarData(var, false);
            dmCluster.addVariant(copySV);
        }

        LinkFinder linkFinder = new LinkFinder();
        dmCluster.setAssemblyLinkedPairs(linkFinder.createAssemblyLinkedPairs(dmCluster));

        ChainFinder chainFinder = new ChainFinder();
        chainFinder.initialise(dmCluster);
        chainFinder.formClusterChains(false);

        if(dmCluster.getChains().size() != 1)
            return null;

        return dmCluster.getChains().get(0);
    }

    private static double DOUBLE_MINUTE_PLOIDY_THRESHOLD = 8;
    private static double DOUBLE_MINUTE_PLOIDY_GAP_RATIO = 3;

    public static void reportDoubleMinutes(final SvCluster cluster, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        // order SVs in descending ploidy order
        List<Double> ploidyList = Lists.newArrayList();
        List<SvVarData> indexSvList = Lists.newArrayList();
        boolean hasHighPloidy = false;

        for(final SvVarData var : cluster.getSVs())
        {
            double ploidy = var.ploidy();
            int i = 0;
            for(; i < ploidyList.size(); ++i)
            {
                Double otherPloidy = ploidyList.get(i);
                if(ploidy > otherPloidy)
                    break;
            }

            ploidyList.add(i, ploidy);
            indexSvList.add(i, var);

            if(ploidy >= DOUBLE_MINUTE_PLOIDY_THRESHOLD)
                hasHighPloidy = true;
        }

        if(!hasHighPloidy)
            return;

        boolean isPotentialDM = false;

        // keep track of any chromosome with high ploidy as an indication of false DMs and/or under-clustering
        List<String> highPloidyChromosomes = Lists.newArrayList();

        int svsAboveThreshold = 0;
        double minPloidyAboveThreshold = 0;
        for(int i = 0; i < ploidyList.size(); ++i)
        {
            Double ploidy = ploidyList.get(i);

            if(ploidy < DOUBLE_MINUTE_PLOIDY_THRESHOLD)
                break;

            ++svsAboveThreshold;
            minPloidyAboveThreshold = ploidy;

            final SvVarData var = indexSvList.get(i);

            for(int be = SE_START; be <= SE_END; ++be)
            {
                boolean useStart = isStart(be);

                if(var.isNullBreakend() && !useStart)
                    continue;

                if(!highPloidyChromosomes.contains(var.chromosome(useStart)))
                    highPloidyChromosomes.add(var.chromosome(useStart));
            }

            // check vs next
            if(i == ploidyList.size() - 1)
            {
                LOGGER.debug(String.format("cluster(%s count=%d) DM highPloidyCount(%d chr=%d) currentSV(%s) ploidy(%.2f) with no others",
                        cluster.id(), cluster.getSvCount(), svsAboveThreshold, highPloidyChromosomes.size(), var.posId(), ploidy));
                isPotentialDM = true;
                break;
            }
            else
            {
                double nextPloidy = ploidyList.get(i+1);

                if(nextPloidy * DOUBLE_MINUTE_PLOIDY_GAP_RATIO < ploidy)
                {
                    LOGGER.debug(String.format("cluster(%s count=%d) DM highPloidyCount(%d chr=%d) currentSV(%s) ploidy(%.2f) vs next(%.3f)",
                            cluster.id(), cluster.getSvCount(), svsAboveThreshold, highPloidyChromosomes.size(),
                            indexSvList.get(i).posId(), ploidy, nextPloidy));
                    isPotentialDM = true;
                    break;
                }
            }
        }

        if(isPotentialDM)
        {
            // check for high ploidy in other variants on the relevant chromosomes
            boolean otherClustersHaveHighPloidy = false;

            for(final String chromsome : highPloidyChromosomes)
            {
                final List<SvBreakend> breakendList = chrBreakendMap.get(chromsome);

                if(breakendList == null)
                    continue;

                for(final SvBreakend breakend : breakendList)
                {
                    if(breakend.getSV().getCluster() == cluster)
                        continue;

                    if(breakend.ploidy() * DOUBLE_MINUTE_PLOIDY_GAP_RATIO >= minPloidyAboveThreshold)
                    {
                        otherClustersHaveHighPloidy = true;
                        break;
                    }
                }

                if(otherClustersHaveHighPloidy)
                    break;
            }

            final String dmAnnotation = otherClustersHaveHighPloidy ? CLUSTER_ANNONTATION_DM + "_Unclear" : CLUSTER_ANNONTATION_DM;
            cluster.addAnnotation(dmAnnotation);
        }
    }

}
