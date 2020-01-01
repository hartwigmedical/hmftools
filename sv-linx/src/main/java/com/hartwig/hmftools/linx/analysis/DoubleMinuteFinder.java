package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.typeAsInt;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStrList;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatPloidy;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getSvTypesStr;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.linx.types.ResolvedType.DOUBLE_MINUTE;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_BFB_AMP;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_DM;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LINK_TYPE_TI;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.linx.chaining.ChainFinder;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.cn.TelomereCentromereCnData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.types.SvArmGroup;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.internal.$guava$.io.$ByteSink;

public class DoubleMinuteFinder
{
    private CnDataLoader mCnAnalyser;
    private SvGeneTranscriptCollection mGeneTransCache;
    private final ChainFinder mChainFinder;

    private final List<Integer> mProcessedClusters;
    private final Map<Integer, SvChain> mClusterChains;
    private final Map<Integer, List<SvVarData>> mClusterSVs;

    private String mOutputDir;
    private BufferedWriter mFileWriter;

    private static final double PLOIDY_THRESHOLD = 8;
    private static final double ADJACENT_PLOIDY_RATIO = 2.3;
    private static final double DOMINANT_FOLDBACK_RATIO_THRESHOLD = 0.75;
    private static final double DOMINANT_FOLDBACK_INF_THRESHOLD = 0.5;
    private static final double DOMINANT_FOLDBACK_PLOIDY_THRESHOLD = 4;
    private static final int INF_PAIR_MIN_DISTANCE = 50000;

    private static final Logger LOGGER = LogManager.getLogger(DoubleMinuteFinder.class);

    public DoubleMinuteFinder()
    {
        mChainFinder = new ChainFinder();
        mCnAnalyser = null;
        mGeneTransCache = null;

        mProcessedClusters = Lists.newArrayList();
        mClusterChains = Maps.newHashMap();
        mClusterSVs = Maps.newHashMap();

        mOutputDir = null;
        mFileWriter = null;
    }

    public void setGeneTransCache(final SvGeneTranscriptCollection geneTransCache) { mGeneTransCache = geneTransCache; }
    public void setCopyNumberAnalyser(CnDataLoader cnAnalyser) { mCnAnalyser = cnAnalyser; }
    public void setOutputDir(final String outputDir)
    {
        mOutputDir = outputDir;
    }
    public final Map<Integer, List<SvVarData>> getClusterSVs() { return mClusterSVs; }

    public void clear()
    {
        mProcessedClusters.clear();
        mClusterSVs.clear();
        mClusterChains.clear();
    }

    public void analyseCluster(SvCluster cluster)
    {
        analyseCluster(cluster, false);
    }

    public void analyseCluster(SvCluster cluster, boolean reassess)
    {
        if(mProcessedClusters.contains(cluster.id()))
        {
            if(!reassess)
                return;

            cluster.setDoubleMinuteData(Lists.newArrayList(), null); // clear any previous DM data
        }
        else
        {
            mProcessedClusters.add(cluster.id());
        }

        if(cluster.getSvCount() == 1 && cluster.getSV(0).type() != DUP)
            return;

        double maxClusterPloidy = 0;
        double maxInfPloidy = 0;

        for(SvVarData var : cluster.getSVs())
        {
            maxClusterPloidy = max(maxClusterPloidy, var.ploidyMin());

            if(var.isInferredSgl())
                maxInfPloidy = max(maxInfPloidy, var.ploidy());
        }

        double sumFbPloidy = 0;
        double maxFbPloidy = 0;

        if(!cluster.getFoldbacks().isEmpty())
        {
            for (SvVarData var : cluster.getFoldbacks())
            {
                sumFbPloidy += var.isChainedFoldback() ? var.ploidy() * 0.5 : var.ploidy();
                maxFbPloidy = max(maxFbPloidy, var.ploidy());
            }
        }

        // check plausibility of BFB explaining max observed ploidy allowing for INFs:
        // - BFB plausible if 2*sum(foldback ploidy) + maxINFPloidy > maxPloidy and maxFBPloidy > 0.15 * maxPloidy
        // - If foldbacks exist but one foldback is a ‘dominant foldback’  (>75% of the total foldback and 0.5 * maxINFPloidy) then BFB is still implausible

        boolean hasDominantFoldback = (maxFbPloidy > DOMINANT_FOLDBACK_RATIO_THRESHOLD * sumFbPloidy)
                && (maxFbPloidy >= DOMINANT_FOLDBACK_PLOIDY_THRESHOLD)
                && (maxFbPloidy > DOMINANT_FOLDBACK_INF_THRESHOLD * maxInfPloidy);

        if(!hasDominantFoldback && sumFbPloidy > 0 && 2 * sumFbPloidy + maxInfPloidy > maxClusterPloidy)
        {
            LOGGER.debug(String.format("cluster(%s) maxPloidy(%.1f) foldbacks(%d ploidy(max=%.1f sum=%.1f) explains AMP",
                    cluster.id(), maxClusterPloidy, cluster.getFoldbacks().size(), maxFbPloidy, sumFbPloidy));

            cluster.addAnnotation(CLUSTER_ANNOT_BFB_AMP);
            return;
        }

        if(maxClusterPloidy < PLOIDY_THRESHOLD)
            return;

        // other the criteria to be a DM are:
        // - at least 1 NON SIMPLE DEL variant has ploidy > 2.3x neigbouring major allele ploidy on at least one side
        // - minPloidy > 8
        // - NOT (INF=2) < 50k bases
        // - NOT amplifying a centromere or telomere

        // cluster satisfies the ploidy requirements - now attempt to find its boundaries, ie the SVs which
        // formed the DM by looking at ploidy and its ratio to adjacent major AP
        final List<SvVarData> highPloidySVs = Lists.newArrayList();
        boolean hasHighMinPloidy = false;

        if(cluster.getSvCount() == 1)
        {
            final SvVarData var = cluster.getSV(0);
            double svAdjMAPRatio = getAdjacentMajorAPRatio(var);
            if(svAdjMAPRatio < ADJACENT_PLOIDY_RATIO)
                return;

            highPloidySVs.add(var);
            hasHighMinPloidy = true;
        }
        else
        {
            for (final SvVarData var : cluster.getSVs())
            {
                if (var.ploidyMax() < maxClusterPloidy)
                    continue;

                if(var.type() != DEL && var.ploidyMin() >= PLOIDY_THRESHOLD)
                    hasHighMinPloidy = true;

                // at least one high-ploidy breakend must have a high ploidy relative to the
                // adjacent CN segment's major allele ploidy and a min ploidy above the threshold, not including DELs
                double svAdjMAPRatio = getAdjacentMajorAPRatio(var);
                boolean hasHighPloidyVsMAP = svAdjMAPRatio >= ADJACENT_PLOIDY_RATIO;

                if(hasHighPloidyVsMAP)
                {
                    highPloidySVs.add(var);
                }
            }
        }

        if(highPloidySVs.isEmpty() || !hasHighMinPloidy)
            return;

        LOGGER.debug(String.format("cluster(%s) possible DM: maxPloidy(%.1f) dmSvCount(%d)",
                cluster.id(), maxClusterPloidy, highPloidySVs.size()));

        if(highPloidySVs.size() == 2 && highPloidySVs.get(0).isInferredSgl() && highPloidySVs.get(1).isInferredSgl())
        {
            long distance = abs(highPloidySVs.get(0).position(true) - highPloidySVs.get(1).position(true));

            if(distance < INF_PAIR_MIN_DISTANCE)
            {
                LOGGER.debug(String.format("cluster(%s) possible DM: inf-pair distance(%d) too small",
                        cluster.id(), distance));
                return;
            }
        }

        final SvChain dmChain = createDMChain(cluster, highPloidySVs);

        // a single DUP or a chain involving all high-ploidy SVs which can be made into a loop
        boolean fullyChained = false;

        if(dmChain != null)
        {
            if(highPloidySVs.size() == 1)
            {
                fullyChained = true; // the single DUP case
            }
            else
            {
                fullyChained = dmChain.getSvCount() == highPloidySVs.size() && dmChain.isClosedLoop();

                if(fullyChained)
                {
                    dmChain.logLinks();
                }
            }
        }
        // check that no arm has a high telomere or centromere relative to the DM's ploidy - make an exception for fully chained DMs
        if(!amplifiedVsSamplePloidy(cluster, maxClusterPloidy) && !fullyChained)
        {
            for (SvArmGroup armGroup : cluster.getArmGroups())
            {
                if (!amplifiedVsArm(cluster, armGroup.chromosome(), armGroup.arm(), maxClusterPloidy))
                    return;
            }
        }


        mClusterChains.put(cluster.id(), dmChain);
        mClusterSVs.put(cluster.id(), highPloidySVs);

        cluster.setDoubleMinuteData(highPloidySVs, fullyChained ? dmChain : null);

        cluster.addAnnotation(CLUSTER_ANNOT_DM);

        if(highPloidySVs.size() == cluster.getSvCount())
        {
            cluster.setResolved(false, DOUBLE_MINUTE);
        }

        if(highPloidySVs.size() == 1 && fullyChained)
        {
            // single DUPs won't go through the chaining routine so cache this chain here
            final SvVarData var = highPloidySVs.get(0);
            dmChain.setPloidyData(var.ploidy(), var.ploidyUncertainty());
            cluster.addChain(dmChain, false);
        }

        LOGGER.debug(String.format("cluster(%s) identified DM: maxPloidy(%.1f) dmSvCount(%d) fullyChained(%s)",
                cluster.id(), maxClusterPloidy, highPloidySVs.size(), fullyChained));
    }

    private static double getAdjacentMajorAPRatio(final SvVarData var)
    {
        // get the largest ratio of ploidy to the adjacent major AP
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

            double adjacentMap = cnData.majorAllelePloidy();

            maxRatio = adjacentMap > 0 ? max(var.ploidy() / adjacentMap, maxRatio) : cnData.CopyNumber;
        }

        return maxRatio;
    }

    private boolean amplifiedVsArm(final SvCluster cluster, final String chromosome, final String arm, double dmPloidy)
    {
        final TelomereCentromereCnData tcData = mCnAnalyser.getChrTeleCentroData().get(chromosome);

        if(tcData == null)
            return true; // ignore if not available

        double centromereCopyNumber = arm == CHROMOSOME_ARM_P ? tcData.CentromerePArm : tcData.CentromereQArm;
        double telomereCopyNumber = arm == CHROMOSOME_ARM_P ? tcData.TelomerePArm : tcData.TelomereQArm;

        if(dmPloidy > centromereCopyNumber * ADJACENT_PLOIDY_RATIO && dmPloidy > telomereCopyNumber * ADJACENT_PLOIDY_RATIO)
        {
            return true;
        }

        LOGGER.debug(String.format("cluster(%s) possible DM: not amplified vs {}:{} telo(%s) centro(%s)",
                cluster.id(), chromosome, arm, formatPloidy(telomereCopyNumber), formatPloidy(telomereCopyNumber)));

        return false;
    }

    private boolean amplifiedVsSamplePloidy(final SvCluster cluster, double dmPloidy)
    {
        final PurityContext samplePurity = mCnAnalyser.getPurityContext();
        if(samplePurity == null)
            return false;

        double samplePloidy = samplePurity.score().maxPloidy();
        if(dmPloidy > 50 * samplePloidy) // effectively disabled for now
        {
            LOGGER.debug(String.format("cluster(%s) DM amplified vs sample(%s)",
                    cluster.id(), formatPloidy(samplePloidy)));

            return true;
        }

        return false;
    }

    private final SvChain createDMChain(SvCluster cluster, List<SvVarData> dmSVList)
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
            chain.setPloidyData(var.ploidy(), var.ploidyUncertainty());
            return chain;
        }

        mChainFinder.initialise(cluster, dmSVList, false);
        mChainFinder.formChains(false);

        if(mChainFinder.getUniqueChains().size() != 1)
            return null;

        SvChain chain = mChainFinder.getUniqueChains().get(0);

        mChainFinder.clear();
        return chain;
    }

    public void reportCluster(final String sampleId, final SvCluster cluster)
    {
        if(mOutputDir == null || mOutputDir.isEmpty())
            return;

        if(!cluster.hasAnnotation(CLUSTER_ANNOT_DM))
            return;

        final List<SvVarData> highPloidySVs = mClusterSVs.get(cluster.id());

        if(highPloidySVs.isEmpty())
            return;

        final SvChain chain = mClusterChains.get(cluster.id());

        // a single DUP or a chain involving all high-ploidy SVs which can be made into a loop
        boolean fullyChained = false;
        boolean chainsCentromere = false;

        if(chain != null)
        {
            if(highPloidySVs.size() == 1)
            {
                fullyChained = true; // the single DUP case
                chainsCentromere = highPloidySVs.get(0).isCrossArm();
            }
            else
            {
                fullyChained = chain.getSvCount() == highPloidySVs.size() && chain.isClosedLoop();

                if(fullyChained)
                {
                    chainsCentromere = chain.getLinkedPairs().stream().anyMatch(x -> x.firstBreakend().arm() != x.secondBreakend().arm());
                }
            }
        }

        String svIds = "";
        int[] typeCounts = new int[StructuralVariantType.values().length];
        final List<String> chromosomes = Lists.newArrayList();

        double minDMPloidy = 0;
        double maxDMPloidy = 0;
        double maxDMCopyNumber = 0;
        double minAdjMAPRatio = 0;
        long minPosition = 0;
        long maxPosition = 0;

        for(final SvVarData var : highPloidySVs)
        {
            ++typeCounts[typeAsInt(var.type())];

            minDMPloidy = max(var.ploidyMin(), minDMPloidy);
            maxDMPloidy = max(var.ploidy(), maxDMPloidy);

            maxDMCopyNumber = max(maxDMCopyNumber, max(var.copyNumber(true), var.copyNumber(false)));

            double svAdjMAPRatio = getAdjacentMajorAPRatio(var);

            if(svAdjMAPRatio > 0 && (minAdjMAPRatio == 0 || svAdjMAPRatio < minAdjMAPRatio))
                minAdjMAPRatio = svAdjMAPRatio;

            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(var.isSglBreakend() && se== SE_END)
                    continue;

                final SvBreakend breakend = var.getBreakend(se);
                final String chromosome = breakend.chromosome();

                if(!chromosomes.contains(chromosome))
                    chromosomes.add(chromosome);

                minPosition = minPosition == 0 ? breakend.position() : min(breakend.position(), minPosition);
                maxPosition = max(breakend.position(), maxPosition);
            }

            svIds = appendStr(svIds, var.idStr(), ';');
        }

        double sumFbPloidy = 0;
        double maxFbPloidy = 0;
        double sumSglPloidy = 0;
        double maxSglPloidy = 0;
        int nonDmSvsFullPloidy = 0;
        int nonDmSvsHalfPloidy = 0;

        for(final SvVarData var : cluster.getSVs())
        {
            if(var.isFoldback())
            {
                sumFbPloidy += var.isChainedFoldback() ? var.ploidy() * 0.5 : var.ploidy();
                maxFbPloidy = max(var.ploidy(), maxFbPloidy);
            }

            if(var.isSglBreakend())
            {
                sumSglPloidy += var.ploidy();
                maxSglPloidy = max(var.ploidy(), maxSglPloidy);
            }

            if(!highPloidySVs.contains(var))
            {
                if(var.ploidyMax() >= minDMPloidy)
                    ++nonDmSvsFullPloidy;

                if(var.ploidyMax() >= minDMPloidy * 0.5)
                    ++nonDmSvsHalfPloidy;
            }
        }

        double maxCentroTeloCopyNumber = 0;

        for(final SvArmGroup armGroup : cluster.getArmGroups())
        {
            final TelomereCentromereCnData tcData = mCnAnalyser.getChrTeleCentroData().get(armGroup.chromosome());

            if (tcData != null)
            {
                if (armGroup.arm() == CHROMOSOME_ARM_P)
                {
                    maxCentroTeloCopyNumber = max(maxCentroTeloCopyNumber, max(tcData.TelomerePArm, tcData.CentromerePArm));
                }
                else
                {
                    maxCentroTeloCopyNumber = max(maxCentroTeloCopyNumber, max(tcData.TelomereQArm, tcData.CentromereQArm));
                }
            }
        }

        final String dmTypesStr = getSvTypesStr(typeCounts);

        final PurityContext purityContext = mCnAnalyser.getPurityContext();
        double samplePurity = purityContext != null ? purityContext.bestFit().purity() : 0;
        double samplePloidy = purityContext != null ? purityContext.bestFit().ploidy() : 0;

        long dmChainLength = chain != null ? chain.getLength(false) : 0;
        int chainSvCount = chain != null ? chain.getSvCount() : 0;

        // get amplified genes list by looking at all section traversed by this chain or breakends with genes in them?
        String amplifiedGenesStr = chain != null ? getAmplifiedGenesList(chain) : "";

        final double chainData[] = chain != null ? getChainCharacteristics(cluster, chain, maxDMPloidy) : null;

        final String chromosomeStr = appendStrList(chromosomes, ';');

        try
        {
            if (mFileWriter == null)
            {
                String outputFileName = mOutputDir;

                outputFileName += "LNX_DOUBLE_MINUTES.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,ClusterId,ClusterDesc,ResolvedType,ClusterCount");
                mFileWriter.write(",SamplePurity,SamplePloidy,DMSvCount,DMSvTypes");
                mFileWriter.write(",FullyChained,ChainLength,ChainCount,SvIds,Chromosomes");
                mFileWriter.write(",MaxCopyNumber,MinPloidy,MaxPloidy,AmpGenes,ChainMinCnPercent,ChainDiffPloidies,ChainNonDMSVs");
                mFileWriter.write(",FbCount,FbSumPloidy,FbMaxPloidy,SglCount,SglSumPloidy,SglMaxPloidy");
                mFileWriter.write(",MinPosition,MaxPosition,MaxTeloCentroCn,NonDmSvsGtPloidy,NonDmSvsGtHalfPloidy,CrossCentro");

                for(Integer cbr : CN_SEGMENT_BUCKETS)
                {
                    mFileWriter.write(String.format(",CNR_%d", cbr));
                }

                mFileWriter.newLine();
            }

            mFileWriter.write(String.format("%s,%d,%s,%s,%d",
                    sampleId, cluster.id(), cluster.getDesc(), cluster.getResolvedType(), cluster.getSvCount()));

            mFileWriter.write(String.format(",%.1f,%.1f,%d,%s",
                    samplePurity, samplePloidy, highPloidySVs.size(), dmTypesStr));

            mFileWriter.write(String.format(",%s,%d,%d,%s,%s",
                    fullyChained, dmChainLength, chainSvCount, svIds, chromosomeStr));

            mFileWriter.write(String.format(",%.1f,%.1f,%.1f,%s",
                    maxDMCopyNumber, minDMPloidy, maxDMPloidy, amplifiedGenesStr));

            if(chainData != null)
            {
                mFileWriter.write(String.format(",%.2f,%.0f,%.0f",
                        chainData[CHAIN_DATA_MIN_CN], chainData[CHAIN_DATA_DIFF_PLOIDIES], chainData[CHAIN_DATA_NON_DM_SVS]));
            }
            else
            {
                mFileWriter.write(",0,0,0");
            }

            mFileWriter.write(String.format(",%d,%.1f,%.1f,%d,%.1f,%.1f",
                    cluster.getFoldbacks().size(), sumFbPloidy, maxFbPloidy,
                    cluster.getSglBreakendCount(), sumSglPloidy, maxSglPloidy));

            mFileWriter.write(String.format(",%d,%d,%.2f,%d,%d,%s",
                    chromosomes.size() == 1 ? minPosition : 0, chromosomes.size() == 1 ? maxPosition : 0,
                    maxCentroTeloCopyNumber, nonDmSvsFullPloidy, nonDmSvsHalfPloidy, chainsCentromere));

            mFileWriter.write(String.format("%s",getCopyNumberSegmentData(cluster, samplePloidy)));

            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing DM data: {}", e.toString());
        }
    }

    private static int CHAIN_DATA_MIN_CN = 0;
    private static int CHAIN_DATA_NON_DM_SVS = 1;
    private static int CHAIN_DATA_DIFF_PLOIDIES = 2;

    private double[] getChainCharacteristics(final SvCluster cluster, final SvChain chain, double maxDMPloidy)
    {
        final double[] chainData = {1.0, 0, 0};

        if(cluster.getSvCount() == 1)
            return chainData;

        double minCopyNumber = maxDMPloidy;

        List<Double> diffPloidies = Lists.newArrayList();
        List<SvVarData> nonDmSVs = Lists.newArrayList();

        for(SvLinkedPair pair : chain.getLinkedPairs())
        {
            final SvBreakend lowerBreakend = pair.getBreakend(true);
            final SvBreakend upperBreakend = pair.getBreakend(false);

            final List<SvBreakend> breakendList = cluster.getChrBreakendMap().get(pair.chromosome());

            for(int i = lowerBreakend.getClusterChrPosIndex() + 1; i <= upperBreakend.getClusterChrPosIndex() - 1; ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvVarData var = breakend.getSV();

                if(!chain.getSvList().contains(var))
                {
                    if (!nonDmSVs.contains(var))
                        nonDmSVs.add(var);

                    if (!diffPloidies.stream().anyMatch(x -> copyNumbersEqual(x, var.ploidy())))
                        diffPloidies.add(var.ploidy());
                }

                // ignore simple consecutive DELs
                if(i < upperBreakend.getClusterChrPosIndex() - 1)
                {
                    final SvBreakend nextBreakend = breakendList.get(i+1);

                    if(breakend.getSV().type() == DEL && breakend.getSV() == nextBreakend.getSV())
                    {
                        ++i;
                        continue;
                    }
                }

                minCopyNumber = min(max(breakend.copyNumberLowSide(), 0), minCopyNumber);
            }
        }

        chainData[CHAIN_DATA_MIN_CN] = minCopyNumber / maxDMPloidy;
        chainData[CHAIN_DATA_DIFF_PLOIDIES] = diffPloidies.size();
        chainData[CHAIN_DATA_NON_DM_SVS] = nonDmSVs.size();

        return chainData;
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

    private static List<Integer> CN_SEGMENT_BUCKETS = Lists.newArrayList(3, 6, 10, 20, 50, 100);

    private static String getCopyNumberSegmentData(final SvCluster cluster, double samplePloidy)
    {
        final long[] cnSegmentLengths = new long[CN_SEGMENT_BUCKETS.size()];

        for(List<SvBreakend> breakendList : cluster.getChrBreakendMap().values())
        {
            long prevPosition = 0;
            double prevCN = 0;
            boolean inSegment = false;
            double netPloidy = 0;

            for(SvBreakend breakend : breakendList)
            {
                if(!inSegment)
                {
                    if(breakend.orientation() == 1)
                    {
                        if(prevPosition > 0)
                        {
                            double avgCopyCN = (prevCN + breakend.copyNumber()) * 0.5;
                            long segmentLength = breakend.position() - prevPosition;
                            addCnSegmentData(cnSegmentLengths, segmentLength, samplePloidy, avgCopyCN);
                        }

                        continue;
                    }

                    inSegment = true;
                    prevPosition = breakend.position();
                    prevCN = breakend.copyNumber();
                    netPloidy = breakend.ploidy();
                }
                else
                {
                    long segmentLength = breakend.position() - prevPosition;

                    if(breakend.orientation() == -1)
                    {
                        // another breakend increasing CN - record segment CN up to this point
                        addCnSegmentData(cnSegmentLengths, segmentLength, samplePloidy, prevCN);

                        // move position on
                        prevPosition = breakend.position();
                        prevCN = breakend.copyNumber();
                        netPloidy += breakend.ploidy();
                    }
                    else
                    {
                        double avgCopyCN = (prevCN + breakend.copyNumber()) * 0.5;
                        addCnSegmentData(cnSegmentLengths, segmentLength, samplePloidy, avgCopyCN);

                        netPloidy -= breakend.ploidy();
                        prevPosition = breakend.position();

                        if(copyNumbersEqual(netPloidy, 0))
                            inSegment = false;
                    }
                }
            }
        }

        String cnSegmentData = "";

        for(int i = 0; i < cnSegmentLengths.length; ++i)
        {
            cnSegmentData += String.format(",%d", cnSegmentLengths[i]);
        }

        return cnSegmentData;
    }

    private static void addCnSegmentData(final long[] cnSegmentLengths, long length, double samplePloidy, double copyNumber)
    {
        // bucket into maximum 5 multiples of sample ploidy
        double cnRatio = max(copyNumber/samplePloidy, 1);

        if(cnRatio < CN_SEGMENT_BUCKETS.get(0))
            return;

        int index = 0;

        for(; index < cnSegmentLengths.length; ++index)
        {
            if(cnRatio < CN_SEGMENT_BUCKETS.get(index))
                break;
        }

        --index;

        cnSegmentLengths[index] += length;
    }

    public void close()
    {
        closeBufferedWriter(mFileWriter);
    }
}
