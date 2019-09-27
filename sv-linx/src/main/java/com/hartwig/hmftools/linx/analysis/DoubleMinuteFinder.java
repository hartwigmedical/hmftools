package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.typeAsInt;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.addSvToChrBreakendMap;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
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
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.linx.chaining.ChainFinder;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.cn.TelomereCentromereCnData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class DoubleMinuteFinder
{
    private CnDataLoader mCnAnalyser;
    private SvGeneTranscriptCollection mGeneTransCache;
    private ChainFinder mChainFinder;

    private List<Integer> mProcessedClusters;
    private Map<Integer, SvChain> mClusterChains;
    private Map<Integer, List<SvVarData>> mClusterSVs;

    private String mOutputDir;
    private BufferedWriter mFileWriter;

    private static final double PLOIDY_THRESHOLD = 8;
    private static final double HIGH_PLOIDY_FACTOR = 1.1;
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

        double maxSvPloidy = 0;
        double maxInfPloidy = 0;

        for(SvVarData var : cluster.getSVs())
        {
            maxSvPloidy = max(maxSvPloidy, var.ploidy());

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

        if(!hasDominantFoldback && sumFoldbackPloidy > 0 && 2 * sumFoldbackPloidy + maxInfPloidy > maxSvPloidy)
        {
            LOGGER.debug(String.format("cluster(%s) maxPloidy(%.1f) foldbacks(%d ploidy(max=%.1f sum=%.1f) explains AMP",
                    cluster.id(), maxSvPloidy, cluster.getFoldbacks().size(), maxFoldbackPloidy, sumFoldbackPloidy));

            cluster.addAnnotation(CLUSTER_ANNOT_BFB_AMP);
            return;
        }

        if(maxSvPloidy < PLOIDY_THRESHOLD)
            return;

        // other the criteria to be a DM are:
        // - at least 1 NON SIMPLE DEL variant has ploidy > 2.3x neigbouring major allele ploidy on at least one side
        // - minPloidy > 8
        // - NOT (INF=2) < 50k bases
        // - NOT amplifying a centromere or telomere

        // cluster satisfies the ploidy requirements - now attempt to find its boundaries, ie the SVs which
        // formed the DM by looking at ploidy and its ratio to adjacent major AP
        List<SvVarData> highPloidySVs = Lists.newArrayList();

        if(cluster.getSvCount() == 1)
        {
            highPloidySVs.add(cluster.getSV(0));
        }
        else
        {
            for (final SvVarData var : cluster.getSVs())
            {
                // TEMP: scale max ploidy as it grows due to uncertainty
                double svMaxPloidy = pow(var.ploidyMax(), HIGH_PLOIDY_FACTOR);

                if (svMaxPloidy >= maxSvPloidy)
                    highPloidySVs.add(var);
            }
        }

        // at least one high-ploidy breakend must have a high ploidy relative to the
        // adjacent CN segment's major allele ploidy and a min ploidy above the threshold, not including DELs
        boolean hasHighMinPloidy = false;

        int index = 0;
        while(index < highPloidySVs.size())
        {
            SvVarData var = highPloidySVs.get(index);

            if(var.type() != DEL && var.ploidyMin() >= PLOIDY_THRESHOLD)
                hasHighMinPloidy = true;

            boolean hasHighPloidyVsMAP = false;

            for (int se = SE_START; se <= SE_END; ++se)
            {
                if(se == SE_END && var.isSglBreakend())
                    continue;

                final SvBreakend breakend = var.getBreakend(se);

                double adjacentMap = getAdjacentMajorAllelePloidy(breakend);

                if(Double.isNaN(adjacentMap))
                    continue;

                // check the major allele ploidy outside this breakend
                if(var.ploidyMin() >= adjacentMap * ADJACENT_PLOIDY_RATIO)
                {
                    hasHighPloidyVsMAP = true;
                    break;
                }
            }

            if(hasHighPloidyVsMAP)
                ++index;
            else
                highPloidySVs.remove(index);
        }

        if(highPloidySVs.isEmpty() || !hasHighMinPloidy)
            return;

        LOGGER.debug(String.format("cluster(%s) possible DM: maxPloidy(%.1f) dmSvCount(%d)",
                cluster.id(), maxSvPloidy, highPloidySVs.size()));

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

        List<String> processedArms = Lists.newArrayList();

        for(SvVarData var : highPloidySVs)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(se == SE_END && (var.isLocal() || var.isSglBreakend()))
                    continue;

                final String chrArm = makeChrArmStr(var.chromosome(isStart(se)), var.arm(isStart(se)));

                if(processedArms.contains(chrArm))
                    continue;

                processedArms.add(chrArm);

                if(!amplifiedVsSampleAndArm(cluster, var.chromosome(isStart(se)), var.arm(isStart(se)), maxSvPloidy))
                {
                    return;
                }
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
                cluster.id(), maxSvPloidy, highPloidySVs.size(), fullyChained));

        // reportPotentialGroup(sampleId, cluster, highPloidySVs, fullyChained, dmChain);
    }

    private static double getAdjacentMajorAllelePloidy(final SvBreakend breakend)
    {
        // gets the CN segment data on the lower side of the breakend (ie opposite to orientation)
        final SvCNData cnData = breakend.getSV().getCopyNumberData(breakend.usesStart(), breakend.orientation() == -1);

        if(cnData == null)
            return Double.NaN;

        return cnData.majorAllelePloidy();
    }

    private boolean amplifiedVsSampleAndArm(final SvCluster cluster, final String chromosome, final String arm, double dmPloidy)
    {
        final TelomereCentromereCnData tcData = mCnAnalyser.getChrTeleCentroData().get(chromosome);

        if(tcData == null)
            return false;

        final PurityContext samplePurity = mCnAnalyser.getPurityContext();
        if(samplePurity == null)
            return false;

        double samplePloidy = samplePurity.score().maxPloidy();
        double centromereCopyNumber = arm == CHROMOSOME_ARM_P ? tcData.CentromerePArm : tcData.CentromereQArm;
        double telomereCopyNumber = arm == CHROMOSOME_ARM_P ? tcData.TelomerePArm : tcData.TelomereQArm;

        if(dmPloidy > centromereCopyNumber * ADJACENT_PLOIDY_RATIO && dmPloidy > telomereCopyNumber * ADJACENT_PLOIDY_RATIO)
        {
            return true;
        }
        else if(dmPloidy > 5 * samplePloidy)
        {
            return true;
        }

        LOGGER.debug(String.format("cluster(%s) possible DM: not amplified vs telo(%s) centro(%s) sample(%s)",
                cluster.id(), formatPloidy(telomereCopyNumber),formatPloidy(telomereCopyNumber), formatPloidy(samplePloidy)));

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
        double maxDMCopyNumber = 0;
        double maxCentroTeloCopyNumber = 0;

        for(final SvVarData var : highPloidySVs)
        {
            ++typeCounts[typeAsInt(var.type())];

            if(minDMPloidy == 0 || var.ploidyMin() < minDMPloidy)
                minDMPloidy = var.ploidyMin();

            maxDMCopyNumber = max(maxDMCopyNumber, max(var.copyNumber(true), var.copyNumber(false)));

            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(var.isSglBreakend() && se== SE_END)
                    continue;

                final String chromosome = var.chromosome(isStart(se));

                if(!chromosomes.contains(chromosome))
                    chromosomes.add(chromosome);

                final TelomereCentromereCnData tcData = mCnAnalyser.getChrTeleCentroData().get(chromosome);

                if(tcData != null)
                {
                    if(var.arm(isStart(se)) == CHROMOSOME_ARM_P)
                    {
                        maxCentroTeloCopyNumber = max(maxCentroTeloCopyNumber, tcData.TelomerePArm);
                        maxCentroTeloCopyNumber = max(maxCentroTeloCopyNumber, tcData.CentromerePArm);
                    }
                    else
                    {
                        maxCentroTeloCopyNumber = max(maxCentroTeloCopyNumber, tcData.TelomereQArm);
                        maxCentroTeloCopyNumber = max(maxCentroTeloCopyNumber, tcData.CentromereQArm);
                    }
                }
            }

            svIds = appendStr(svIds, var.idStr(), ';');
        }

        final String dmTypesStr = getSvTypesStr(typeCounts);

        final PurityContext purityContext = mCnAnalyser.getPurityContext();
        double samplePurity = purityContext != null ? purityContext.bestFit().purity() : 0;
        double samplePloidy = purityContext != null ? purityContext.bestFit().ploidy() : 0;

        long dmChainLength = chain != null ? chain.getLength(false) : 0;
        int chainSvCount = chain != null ? chain.getSvCount() : 0;

        // get amplified genes list by looking at all section traversed by this chain or breakends with genes in them?
        String amplifiedGenesStr = chain != null ? getAmplifiedGenesList(chain) : "";

        String chromosomeStr = "";
        for(String chr : chromosomes)
        {
            chromosomeStr = appendStr(chromosomeStr, chr, ';');
        }

        try
        {
            if (mFileWriter == null)
            {
                String outputFileName = mOutputDir;

                outputFileName += "LNX_DOUBLE_MINUTES.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,ClusterId,ClusterDesc,ResolvedType,ClusterCount,SamplePurity,SamplePloidy,DMSvCount,DMSvTypes");
                mFileWriter.write(",FullyChained,ChainLength,ChainCount,SvIds,Chromosomes,MaxTeloCentroCn");
                mFileWriter.write(",MaxCopyNumber,MinPloidy,AmpGenes");
                mFileWriter.newLine();
            }

            mFileWriter.write(String.format("%s,%d,%s,%s,%d",
                    sampleId, cluster.id(), cluster.getDesc(), cluster.getResolvedType(), cluster.getSvCount()));

            mFileWriter.write(String.format(",%.2f,%.2f,%d,%s,%s,%d,%d",
                    samplePurity, samplePloidy, highPloidySVs.size(), dmTypesStr, fullyChained, dmChainLength, chainSvCount));

            mFileWriter.write(String.format(",%s,%s,%.2f,%.2f,%.2f,%s",
                    svIds, chromosomeStr, maxCentroTeloCopyNumber, maxDMCopyNumber, minDMPloidy, amplifiedGenesStr));

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

    public void close()
    {
        closeBufferedWriter(mFileWriter);
    }


}
