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
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatPloidy;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getSvTypesStr;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.linx.types.ResolvedType.DOUBLE_MINUTE;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_BFB_AMP;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_DM;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LINK_TYPE_TI;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

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
            maxClusterPloidy = max(maxClusterPloidy, var.ploidy());

            if(var.isInferredSgl())
                maxInfPloidy = max(maxInfPloidy, var.ploidy());
        }

        double sumFoldbackPloidy = 0;
        double maxFoldbackPloidy = 0;

        if(!cluster.getFoldbacks().isEmpty())
        {
            for (SvVarData var : cluster.getFoldbacks())
            {
                maxFoldbackPloidy = max(maxFoldbackPloidy, var.ploidy());
                sumFoldbackPloidy += var.ploidy();
            }
        }

        // check plausibility of BFB explaining max observed ploidy allowing for INFs:
        // - BFB plausible if 2*sum(foldback ploidy) + maxINFPloidy > maxPloidy and maxFBPloidy > 0.15 * maxPloidy
        // - If foldbacks exist but one foldback is a ‘dominant foldback’  (>75% of the total foldback and 0.5 * maxINFPloidy) then BFB is still implausible

        boolean hasDominantFoldback = (maxFoldbackPloidy > DOMINANT_FOLDBACK_RATIO_THRESHOLD * sumFoldbackPloidy)
                && (maxFoldbackPloidy >= DOMINANT_FOLDBACK_PLOIDY_THRESHOLD)
                && (maxFoldbackPloidy > DOMINANT_FOLDBACK_INF_THRESHOLD * maxInfPloidy);

        if(!hasDominantFoldback && sumFoldbackPloidy > 0 && 2 * sumFoldbackPloidy + maxInfPloidy > maxClusterPloidy)
        {
            LOGGER.debug(String.format("cluster(%s) maxPloidy(%.1f) foldbacks(%d ploidy(max=%.1f sum=%.1f) explains AMP",
                    cluster.id(), maxClusterPloidy, cluster.getFoldbacks().size(), maxFoldbackPloidy, sumFoldbackPloidy));

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
            highPloidySVs.add(cluster.getSV(0));
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

        // check that no arm has a high telomere or centromere relative to the DM's ploidy
        if(!amplifiedVsSamplePloidy(cluster, maxClusterPloidy))
        {
            for (SvArmGroup armGroup : cluster.getArmGroups())
            {
                if (!amplifiedVsArm(cluster, armGroup.chromosome(), armGroup.arm(), maxClusterPloidy))
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

            maxRatio = adjacentMap > 0 ? max(var.ploidyMin() / adjacentMap, maxRatio) : cnData.CopyNumber;
        }

        return maxRatio;
    }

    private boolean amplifiedVsArm(final SvCluster cluster, final String chromosome, final String arm, double dmPloidy)
    {
        final TelomereCentromereCnData tcData = mCnAnalyser.getChrTeleCentroData().get(chromosome);

        if(tcData == null)
            return false;

        double centromereCopyNumber = arm == CHROMOSOME_ARM_P ? tcData.CentromerePArm : tcData.CentromereQArm;
        double telomereCopyNumber = arm == CHROMOSOME_ARM_P ? tcData.TelomerePArm : tcData.TelomereQArm;

        if(dmPloidy > centromereCopyNumber * ADJACENT_PLOIDY_RATIO && dmPloidy > telomereCopyNumber * ADJACENT_PLOIDY_RATIO)
        {
            return true;
        }

        LOGGER.debug(String.format("cluster(%s) possible DM: not amplified vs telo(%s) centro(%s)",
                cluster.id(), formatPloidy(telomereCopyNumber),formatPloidy(telomereCopyNumber)));

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
            return true;
        }

        LOGGER.debug(String.format("cluster(%s) possible DM: not amplified vs sample(%s)",
                cluster.id(), formatPloidy(samplePloidy)));

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

        if(chain != null)
        {
            if(highPloidySVs.size() == 1)
                fullyChained = true; // the single DUP case
            else
                fullyChained = chain.getSvCount() == highPloidySVs.size() && chain.isClosedLoop();
        }

        String svIds = "";
        int[] typeCounts = new int[StructuralVariantType.values().length];
        List<String> chromosomes = Lists.newArrayList();

        double minDMPloidy = 0;
        double maxDMPloidy = 0;
        double maxDMCopyNumber = 0;
        double minAdjMAPRatio = 0;
        double maxCentroTeloCopyNumber = 0;
        long minPosition = 0;
        long maxPosition = 0;

        for(final SvVarData var : highPloidySVs)
        {
            ++typeCounts[typeAsInt(var.type())];

            if(minDMPloidy == 0 || var.ploidyMin() < minDMPloidy)
                minDMPloidy = var.ploidyMin();

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

                final TelomereCentromereCnData tcData = mCnAnalyser.getChrTeleCentroData().get(chromosome);

                if(tcData != null)
                {
                    if(var.arm(isStart(se)) == CHROMOSOME_ARM_P)
                    {
                        maxCentroTeloCopyNumber = max(maxCentroTeloCopyNumber, max(tcData.TelomerePArm, tcData.CentromerePArm));
                    }
                    else
                    {
                        maxCentroTeloCopyNumber = max(maxCentroTeloCopyNumber, max(tcData.TelomereQArm, tcData.CentromereQArm));
                    }
                }
            }

            svIds = appendStr(svIds, var.idStr(), ';');
        }

        /*
        // check for evidence of additional SVs which ought to be part of the DM, or cast it into doubt but exhibiting a spectrum of ploidy drops
        int nearHighPloidySvs = 0;
        double maxNonDmSvAdjMAPRatio = 0;
        double reducedRatio = 0.8 * ADJACENT_PLOIDY_RATIO;

        for(final SvVarData var : cluster.getSVs())
        {
            if(var.type() == DEL || var.ploidyMin() < PLOIDY_THRESHOLD || highPloidySVs.contains(var))
                continue;

            double svAdjMAPRatio = getAdjacentMajorAPRatio(var);
            maxNonDmSvAdjMAPRatio = max(svAdjMAPRatio, maxNonDmSvAdjMAPRatio);

            if(svAdjMAPRatio >= reducedRatio)
                ++nearHighPloidySvs;
        }
        */

        final String dmTypesStr = getSvTypesStr(typeCounts);

        final PurityContext purityContext = mCnAnalyser.getPurityContext();
        double samplePurity = purityContext != null ? purityContext.bestFit().purity() : 0;
        double samplePloidy = purityContext != null ? purityContext.bestFit().ploidy() : 0;

        long dmChainLength = chain != null ? chain.getLength(false) : 0;
        int chainSvCount = chain != null ? chain.getSvCount() : 0;

        // get amplified genes list by looking at all section traversed by this chain or breakends with genes in them?
        String amplifiedGenesStr = chain != null ? getAmplifiedGenesList(chain) : "";

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
                mFileWriter.write(",MaxCopyNumber,MinPloidy,MaxPloidy");
                mFileWriter.write(",AmpGenes,Foldbacks,MinPosition,MaxPosition,MaxTeloCentroCn");

                for(Integer cbr : CN_SEGMENT_BUCKETS)
                {
                    mFileWriter.write(String.format(",CNBR_%d", cbr));
                }

                mFileWriter.newLine();
            }

            mFileWriter.write(String.format("%s,%d,%s,%s,%d",
                    sampleId, cluster.id(), cluster.getDesc(), cluster.getResolvedType(), cluster.getSvCount()));

            mFileWriter.write(String.format(",%.2f,%.2f,%d,%s",
                    samplePurity, samplePloidy, highPloidySVs.size(), dmTypesStr));

            mFileWriter.write(String.format(",%s,%d,%d,%s,%s",
                    fullyChained, dmChainLength, chainSvCount, svIds, chromosomeStr));

            mFileWriter.write(String.format(",%.2f,%.2f,%.2f",
                    maxDMCopyNumber, minDMPloidy, maxDMPloidy));

            mFileWriter.write(String.format(",%s,%d,%d,%d,%.2f",
                    amplifiedGenesStr, cluster.getFoldbacks().size(),
                    chromosomes.size() == 1 ? minPosition : 0, chromosomes.size() == 1 ? maxPosition : 0, maxCentroTeloCopyNumber));

            mFileWriter.write(String.format("%s",getCopyNumberSegmentData(cluster, samplePloidy)));

            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing DM data: {}", e.toString());
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

    private static List<Integer> CN_SEGMENT_BUCKETS = Lists.newArrayList(3, 6, 10, 20, 50, 100);

    private static String getCopyNumberSegmentData(final SvCluster cluster, double samplePloidy)
    {
        final long[] cnSegmentLengths = new long[CN_SEGMENT_BUCKETS.size()];

        for(List<SvBreakend> breakendList : cluster.getChrBreakendMap().values())
        {
            long prevPosition = 0;
            double prevCN = 0;
            boolean inSegment = false;

            for(SvBreakend breakend : breakendList)
            {
                if(!inSegment)
                {
                    if(breakend.orientation() == 1)
                    {
                        if(prevPosition > 0)
                        {
                            double avgCopyCN = (prevCN + breakend.copyNumber()) * 0.5;
                            addCnSegmentData(cnSegmentLengths, breakend.position() - prevPosition, samplePloidy, avgCopyCN);
                        }

                        continue;
                    }

                    inSegment = true;
                    prevPosition = breakend.position();
                    prevCN = breakend.copyNumber();
                }
                else
                {
                    if(breakend.orientation() == -1)
                    {
                        // another breakend increasing CN - record segment CN up to this point
                        addCnSegmentData(cnSegmentLengths, breakend.position() - prevPosition, samplePloidy, prevCN);

                        // move position on
                        prevPosition = breakend.position();
                        prevCN = breakend.copyNumber();
                    }
                    else
                    {
                        double avgCopyCN = (prevCN + breakend.copyNumber()) * 0.5;
                        addCnSegmentData(cnSegmentLengths, breakend.position() - prevPosition, samplePloidy, avgCopyCN);
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
