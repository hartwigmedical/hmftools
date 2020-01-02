package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;

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
import com.hartwig.hmftools.linx.types.DoubleMinuteData;
import com.hartwig.hmftools.linx.types.SvArmGroup;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class DoubleMinuteFinder
{
    private CnDataLoader mCnAnalyser;
    private SvGeneTranscriptCollection mGeneTransCache;
    private final ChainFinder mChainFinder;

    private final List<Integer> mProcessedClusters;
    private final Map<Integer, DoubleMinuteData> mDoubleMinutes;

    private String mOutputDir;
    private BufferedWriter mFileWriter;

    private static final double PLOIDY_THRESHOLD = 8;
    private static final double ADJACENT_PLOIDY_RATIO = 2.5;
    private static final int INF_PAIR_MIN_DISTANCE = 50000;

    /*
    private static final double DOMINANT_FOLDBACK_RATIO_THRESHOLD = 0.75;
    private static final double DOMINANT_FOLDBACK_INF_THRESHOLD = 0.5;
    private static final double DOMINANT_FOLDBACK_PLOIDY_THRESHOLD = 4;
   */

    private static final Logger LOGGER = LogManager.getLogger(DoubleMinuteFinder.class);

    public DoubleMinuteFinder()
    {
        mChainFinder = new ChainFinder();
        mCnAnalyser = null;
        mGeneTransCache = null;

        mProcessedClusters = Lists.newArrayList();
        mDoubleMinutes = Maps.newHashMap();

        mOutputDir = null;
        mFileWriter = null;
    }

    public void setGeneTransCache(final SvGeneTranscriptCollection geneTransCache) { mGeneTransCache = geneTransCache; }
    public void setCopyNumberAnalyser(CnDataLoader cnAnalyser) { mCnAnalyser = cnAnalyser; }
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

        if(cluster.getSvCount() == 1 && cluster.getSV(0).type() != DUP) // early exit - only consider cluster-1 if a DUP
            return;

        double maxClusterPloidy = 0;
        double maxInfPloidy = 0;

        for(SvVarData var : cluster.getSVs())
        {
            if(var.isSglBreakend())
                maxClusterPloidy = max(maxClusterPloidy, var.ploidyMin() * 0.5); // in case the SGL is a disguised foldback
            else
                maxClusterPloidy = max(maxClusterPloidy, var.ploidyMin());

            if(var.isInferredSgl())
                maxInfPloidy = max(maxInfPloidy, var.ploidy());
        }

        // cluster satisfies the ploidy requirements - now attempt to find its boundaries, ie the SVs which
        // formed the DM by looking at ploidy and its ratio to adjacent major AP
        final List<SvVarData> candidateDMSVs = Lists.newArrayList();

        if(cluster.getSvCount() == 1)
        {
            final SvVarData var = cluster.getSV(0);

            if(var.ploidyMin() < PLOIDY_THRESHOLD)
                return;

            double svAdjMAPRatio = getAdjacentMajorAPRatio(var);
            if(svAdjMAPRatio < ADJACENT_PLOIDY_RATIO)
                return;

            candidateDMSVs.add(var);
        }
        else
        {
            boolean hasHighMinPloidy = false;

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
                    candidateDMSVs.add(var);
                }
            }

            if(!hasHighMinPloidy)
                return;
        }

        if(candidateDMSVs.isEmpty())
            return;

        // determine whether foldbacks and/or SGLs could explain the amplification
        double sumFbPloidy = 0;
        double maxFbPloidy = 0;
        int foldbackCount = 0;
        final List<Double> maxSglPloidies = Lists.newArrayList();

        for (final SvVarData var : cluster.getSVs())
        {
            if(candidateDMSVs.contains(var))
                continue;

            if(var.isFoldback())
            {
                sumFbPloidy += var.isChainedFoldback() ? var.ploidy() * 0.5 : var.ploidy();
                foldbackCount += var.isChainedFoldback() ? 0.5 : 1;
                maxFbPloidy = max(maxFbPloidy, var.ploidy());
            }
            else if(var.isSglBreakend())
            {
                // limit the SGLs to the top 2 by ploidy
                captureMaxTwoSglPloidies(var.ploidy(), maxSglPloidies);
            }
        }

        // maximum plausible ploidy from BFB = min(2*Sum(FB ploidy) + (2 maximum INF/SGL ploidies), max(4*max(FB ploidy), 2*max(SGL /INF ploidy),HighestTelo/CentromereCN*2^(FBCount+min(2,INF/SGLCount))
        double sumSglPloidy = maxSglPloidies.stream().mapToDouble(x -> x).sum();
        double maxFbOrSglPloidy = max(6 * maxFbPloidy, maxSglPloidies.isEmpty() ? 0 : 3 * maxSglPloidies.get(0));

        double maxArmEndCopyNumber = getMaxArmEndCopyNumber(cluster);
        double armBfbCopyNumber = maxArmEndCopyNumber * pow(2, (foldbackCount + maxSglPloidies.size()));

        double maxBFBPloidy = min(2 * sumFbPloidy + sumSglPloidy, min(maxFbOrSglPloidy, armBfbCopyNumber));

        LOGGER.debug(String.format("cluster(%s) possible DM: maxPloidy(%.1f) dmSvCount(%d) maxBFBPloidy(%.1f fb=%.1f sgl=%.1f arm=%.1f)",
                cluster.id(), maxClusterPloidy, candidateDMSVs.size(), maxBFBPloidy, sumFbPloidy, sumSglPloidy, armBfbCopyNumber));

        if(maxBFBPloidy > maxClusterPloidy)
        {
            cluster.addAnnotation(CLUSTER_ANNOT_BFB_AMP);
            return;
        }

        /*
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

        */

        // other the criteria to be a DM are:
        // - NOT (INF=2) < 50k bases
        // - NOT amplifying a centromere or telomere

        if(candidateDMSVs.size() == 2 && candidateDMSVs.get(0).isInferredSgl() && candidateDMSVs.get(1).isInferredSgl())
        {
            long distance = abs(candidateDMSVs.get(0).position(true) - candidateDMSVs.get(1).position(true));

            if(distance < INF_PAIR_MIN_DISTANCE)
            {
                LOGGER.debug(String.format("cluster(%s) possible DM: inf-pair distance(%d) too small",
                        cluster.id(), distance));
                return;
            }
        }

        final SvChain dmChain = createDMChain(cluster, candidateDMSVs);

        // a single DUP or a chain involving all high-ploidy SVs which can be made into a loop
        boolean fullyChained = false;

        if(dmChain != null)
        {
            if(candidateDMSVs.size() == 1)
            {
                fullyChained = true; // the single DUP case
            }
            else
            {
                if(candidateDMSVs.size() == 2 && candidateDMSVs.get(0).isSglBreakend() && candidateDMSVs.get(1).isSglBreakend())
                {
                    // pair of facing SGLs - likely a DUP
                    fullyChained = true;
                }
                else
                {
                    fullyChained = dmChain.getSvCount() == candidateDMSVs.size() && dmChain.isClosedLoop();
                    dmChain.logLinks();
                }
            }
        }
        // check that no arm has a high telomere or centromere relative to the DM's ploidy - make an exception for fully chained DMs
        if(!amplifiedVsSamplePloidy(cluster, maxClusterPloidy) && !fullyChained)
        {
            if(maxClusterPloidy < maxArmEndCopyNumber * ADJACENT_PLOIDY_RATIO)
            {
                LOGGER.debug("cluster({}}) possible DM: not amplified vs max armEndCopyNumber({}}) centro({}})",
                        cluster.id(), formatPloidy(maxArmEndCopyNumber));
                return;
            }
        }

        DoubleMinuteData dmData = new DoubleMinuteData(cluster, candidateDMSVs);
        dmData.MaxBFBPloidy = maxBFBPloidy;

        if(dmChain != null)
        {
            dmData.Chains.add(dmChain);
            dmData.FullyChained = fullyChained;
        }

        // collect up other possible DM SVs for logging only for now

        mDoubleMinutes.put(cluster.id(), dmData);

        cluster.setDoubleMinuteData(candidateDMSVs, fullyChained ? dmChain : null);

        cluster.addAnnotation(CLUSTER_ANNOT_DM);

        if(candidateDMSVs.size() == cluster.getSvCount())
        {
            cluster.setResolved(false, DOUBLE_MINUTE);
        }

        if(candidateDMSVs.size() == 1 && fullyChained)
        {
            // single DUPs won't go through the chaining routine so cache this chain here
            final SvVarData var = candidateDMSVs.get(0);
            dmChain.setPloidyData(var.ploidy(), var.ploidyUncertainty());
            cluster.addChain(dmChain, false);
        }

        LOGGER.debug(String.format("cluster(%s) identified DM: maxPloidy(%.1f) dmSvCount(%d) fullyChained(%s)",
                cluster.id(), maxClusterPloidy, candidateDMSVs.size(), fullyChained));
    }

    private static void captureMaxTwoSglPloidies(double newPloidy, final List<Double> maxSglPloidies)
    {
        if(maxSglPloidies.isEmpty())
        {
            maxSglPloidies.add(newPloidy);
        }
        else if(maxSglPloidies.size() == 1)
        {
            if(maxSglPloidies.get(0) > newPloidy)
                maxSglPloidies.add(newPloidy);
            else
                maxSglPloidies.add(0, newPloidy);
        }
        else if(newPloidy > maxSglPloidies.get(1))
        {
            maxSglPloidies.remove(1);

            if(maxSglPloidies.get(0) > newPloidy)
                maxSglPloidies.add(newPloidy);
            else
                maxSglPloidies.add(0, newPloidy);
        }
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

    private double getMaxArmEndCopyNumber(final SvCluster cluster)
    {
        double maxCopyNumber = 0;

        for (SvArmGroup armGroup : cluster.getArmGroups())
        {
            final TelomereCentromereCnData tcData = mCnAnalyser.getChrTeleCentroData().get(armGroup.chromosome());

            if(tcData != null)
            {
                double centromereCopyNumber = armGroup.arm() == CHROMOSOME_ARM_P ? tcData.CentromerePArm : tcData.CentromereQArm;
                double telomereCopyNumber = armGroup.arm() == CHROMOSOME_ARM_P ? tcData.TelomerePArm : tcData.TelomereQArm;
                maxCopyNumber = max(max(telomereCopyNumber, centromereCopyNumber), maxCopyNumber);
            }
        }

        return maxCopyNumber;
    }

    private boolean amplifiedVsSamplePloidy(final SvCluster cluster, double dmPloidy)
    {
        final PurityContext samplePurity = mCnAnalyser.getPurityContext();
        if(samplePurity == null)
            return false;

        double samplePloidy = samplePurity.score().maxPloidy();
        if(dmPloidy > 50 * samplePloidy) // effectively disabled for now
        {
            LOGGER.debug("cluster({}}) DM amplified vs sample({}})",
                    cluster.id(), formatPloidy(samplePloidy));

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

        final DoubleMinuteData dmData = mDoubleMinutes.get(cluster.id());

        if(dmData == null)
            return;

        final SvChain chain = dmData.Chains.isEmpty() ? null : dmData.Chains.get(0);

        // a single DUP or a chain involving all high-ploidy SVs which can be made into a loop
        boolean chainsCentromere = false;

        if(chain != null)
        {
            if(dmData.SVs.size() == 1)
            {
                chainsCentromere = dmData.SVs.get(0).isCrossArm();
            }
            else
            {
                if(dmData.FullyChained)
                {
                    chainsCentromere = chain.getLinkedPairs().stream().anyMatch(x -> x.firstBreakend().arm() != x.secondBreakend().arm());
                }
            }
        }

        String svIds = "";
        final int[] typeCounts = new int[StructuralVariantType.values().length];
        final List<String> chromosomes = Lists.newArrayList();

        double minDMPloidy = 0;
        double maxDMPloidy = 0;
        double maxDMCopyNumber = 0;
        long minPosition = 0;
        long maxPosition = 0;

        for(final SvVarData var : dmData.SVs)
        {
            ++typeCounts[typeAsInt(var.type())];

            minDMPloidy = max(var.ploidyMin(), minDMPloidy);
            maxDMPloidy = max(var.ploidy(), maxDMPloidy);

            maxDMCopyNumber = max(maxDMCopyNumber, max(var.copyNumber(true), var.copyNumber(false)));

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

            if(!dmData.SVs.contains(var))
            {
                if(var.ploidyMax() >= minDMPloidy)
                    ++nonDmSvsFullPloidy;

                if(var.ploidyMax() >= minDMPloidy * 0.5)
                    ++nonDmSvsHalfPloidy;

                if((var.ploidy() >= dmData.MaxBFBPloidy || var.ploidy() > PLOIDY_THRESHOLD)
                && getAdjacentMajorAPRatio(var) >= ADJACENT_PLOIDY_RATIO)
                {
                    dmData.CandidateSVs.add(var);
                }
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
        final String amplifiedGenesStr = chain != null ? getAmplifiedGenesList(chain) : "";

        final double chainData[] = chain != null ? getChainCharacteristics(cluster, chain, maxDMPloidy) : null;

        final String chromosomeStr = appendStrList(chromosomes, ';');

        String possibleSvTypes = "";

        if(!dmData.CandidateSVs.isEmpty())
        {
            final int[] possibleTypeCounts = new int[StructuralVariantType.values().length];

            for (final SvVarData var : dmData.CandidateSVs)
            {
                ++possibleTypeCounts[typeAsInt(var.type())];
            }

            possibleSvTypes = getSvTypesStr(possibleTypeCounts);
        }

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
                mFileWriter.write(",MaxBFBPloidy,FbCount,FbSumPloidy,FbMaxPloidy,SglCount,SglSumPloidy,SglMaxPloidy");
                mFileWriter.write(",MinPosition,MaxPosition,MaxTeloCentroCn,CrossCentro");
                mFileWriter.write(",PossibleSVs,NonDmSvsGtPloidy,NonDmSvsGtHalfPloidy");

                for(Integer cbr : CN_SEGMENT_BUCKETS)
                {
                    mFileWriter.write(String.format(",CNR_%d", cbr));
                }

                mFileWriter.newLine();
            }

            mFileWriter.write(String.format("%s,%d,%s,%s,%d",
                    sampleId, cluster.id(), cluster.getDesc(), cluster.getResolvedType(), cluster.getSvCount()));

            mFileWriter.write(String.format(",%.1f,%.1f,%d,%s",
                    samplePurity, samplePloidy, dmData.SVs.size(), dmTypesStr));

            mFileWriter.write(String.format(",%s,%d,%d,%s,%s",
                    dmData.FullyChained, dmChainLength, chainSvCount, svIds, chromosomeStr));

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

            mFileWriter.write(String.format(",%.1f,%d,%.1f,%.1f,%d,%.1f,%.1f",
                    dmData.MaxBFBPloidy, cluster.getFoldbacks().size(), sumFbPloidy, maxFbPloidy,
                    cluster.getSglBreakendCount(), sumSglPloidy, maxSglPloidy));

            mFileWriter.write(String.format(",%d,%d,%.1f,%s",
                    chromosomes.size() == 1 ? minPosition : 0, chromosomes.size() == 1 ? maxPosition : 0,
                    maxCentroTeloCopyNumber, chainsCentromere));

            mFileWriter.write(String.format(",%s,%d,%d",
                    possibleSvTypes, nonDmSvsFullPloidy, nonDmSvsHalfPloidy));

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
