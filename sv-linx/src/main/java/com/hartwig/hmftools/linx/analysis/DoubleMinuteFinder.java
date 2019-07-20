package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.typeAsInt;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.areLinkedSection;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_CENTROMERE;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getSvTypesStr;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.CENTROMERE_CN;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.P_ARM_TELOMERE_CN;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.Q_ARM_TELOMERE_CN;
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
    private static final double DOMINANT_FOLDBACK_THRESHOLD = 0.75;
    private static final double DOMINANT_FOLDBACK_INF_THRESHOLD = 0.5;
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

            if(var.isNoneSegment())
                maxInfPloidy = max(maxInfPloidy, var.ploidy());
        }

        if(maxSvPloidy < PLOIDY_THRESHOLD)
            return;

        double sumFoldbackPloidy = 0;
        double maxFoldbackPloidy = 0;

        if(!cluster.getFoldbacks().isEmpty())
        {
            List<SvVarData> chainedFoldbacks = Lists.newArrayList();

            for (SvVarData var : cluster.getFoldbacks())
            {
                if(var.isChainedFoldback() && !chainedFoldbacks.contains(var))
                {
                    chainedFoldbacks.add(var);
                    chainedFoldbacks.add(var.getChainedFoldbackSv());
                }

                maxFoldbackPloidy = max(maxFoldbackPloidy, var.ploidy());
                sumFoldbackPloidy += var.ploidy();
            }
        }

        // check plausibility of BFB explaining max observed ploidy allowing for INFs:
        // - BFB plausible if 2*sum(foldback ploidy) + maxINFPloidy > maxPloidy and maxFBPloidy > 0.15 * maxPloidy
        // - If foldbacks exist but one foldback is a ‘dominant foldback’  (>75% of the total foldback and 0.5 * maxINFPloidy) then BFB is still implausible

        boolean hasDominantFoldback = (maxFoldbackPloidy > DOMINANT_FOLDBACK_THRESHOLD * sumFoldbackPloidy)
                && (maxFoldbackPloidy > DOMINANT_FOLDBACK_INF_THRESHOLD * maxInfPloidy);

        if(!hasDominantFoldback)
        {
            if(2 * sumFoldbackPloidy + maxInfPloidy > maxSvPloidy)
            {
                LOGGER.debug(String.format("cluster(%s) maxPloidy(%.1f) foldbacks(%d ploidy(max=%.1f sum=%.1f) explains AMP",
                        cluster.id(), maxSvPloidy, cluster.getFoldbacks().size(), maxFoldbackPloidy, sumFoldbackPloidy));

                cluster.addAnnotation(CLUSTER_ANNOT_BFB_AMP);
                return;
            }
        }

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
                if(se == SE_END && var.isNullBreakend())
                    continue;

                final SvBreakend breakend = var.getBreakend(isStart(se));

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

        if(highPloidySVs.size() == 1 && highPloidySVs.get(0).type() != DUP)
            return;

        if(highPloidySVs.size() == 2 && highPloidySVs.get(0).isNoneSegment() && highPloidySVs.get(1).isNoneSegment())
        {
            long distance = abs(highPloidySVs.get(0).position(true) - highPloidySVs.get(1).position(true));

            if(distance < INF_PAIR_MIN_DISTANCE)
                return;
        }

        boolean hasTorCAmplication = false;

        for(SvVarData var : highPloidySVs)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(se == SE_END && var.isLocal() || var.isNullBreakend())
                    continue;

                if(hasTelomereOrCentromereAmplication(var.chromosome(isStart(se)), var.arm(isStart(se))))
                {
                    hasTorCAmplication = true;
                    break;
                }
            }
        }

        if(hasTorCAmplication)
            return;

        final SvChain dmChain = createDMChain(cluster, highPloidySVs);

        // a single DUP or a chain involving all high-ploidy SVs which can be made into a loop
        boolean fullyChained = false;

        if(dmChain != null)
        {
            if(highPloidySVs.size() == 1)
                fullyChained = true; // the single DUP case
            else
                fullyChained = dmChain.getSvCount() == highPloidySVs.size() && dmChain.isClosedLoop();
        }

        mClusterChains.put(cluster.id(), dmChain);
        mClusterSVs.put(cluster.id(), highPloidySVs);

        if(fullyChained)
            cluster.setDoubleMinuteSVs(highPloidySVs);

        cluster.addAnnotation(CLUSTER_ANNOT_DM);

        if(highPloidySVs.size() == 1 && fullyChained)
        {
            // single DUPs won't go through the chaining routine so cache this chain here
            cluster.addChain(dmChain, false);
        }

        LOGGER.debug(String.format("cluster(%s) maxPloidy(%.1f) dmSvCount(%d) fullyChained(%s)",
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

    private static final double PLOIDY_VS_T_AND_C = 3;

    private boolean hasTelomereOrCentromereAmplication(final String chromosome, final String arm)
    {
        final double[] cnData = mCnAnalyser.getChrCopyNumberMap().get(chromosome);

        if(cnData == null)
            return false;

        //  we want to check telomere / centromere CN < max (8, 3x sample ploidy).

        final PurityContext samplePurity = mCnAnalyser.getPurityContext();
        if(samplePurity == null)
            return false;

        double telomereCN = (arm == CHROMOSOME_ARM_CENTROMERE) ? cnData[P_ARM_TELOMERE_CN] : cnData[Q_ARM_TELOMERE_CN];

        if(cnData[CENTROMERE_CN] > PLOIDY_THRESHOLD || telomereCN > PLOIDY_THRESHOLD)
            return true;

        if(telomereCN > samplePurity.score().maxPloidy() * PLOIDY_VS_T_AND_C
        || cnData[CENTROMERE_CN] > samplePurity.score().maxPloidy() * PLOIDY_VS_T_AND_C)
            return true;

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
            return chain;
        }

        mChainFinder.initialise(cluster, dmSVList);
        mChainFinder.formChains(false);

        if(mChainFinder.getUniqueChains().size() != 1)
            return null;

        SvChain chain = mChainFinder.getUniqueChains().get(0);

        // check whether the chain could form a loop
        SvBreakend chainStart = chain.getOpenBreakend(true);
        SvBreakend chainEnd = chain.getOpenBreakend(false);

        if(chainStart != null && !chainStart.getSV().isNullBreakend() && chainEnd != null && !chainEnd.getSV().isNullBreakend())
        {
            if (areLinkedSection(chainStart.getSV(), chainEnd.getSV(), chainStart.usesStart(), chainEnd.usesStart()))
            {
                SvLinkedPair pair = SvLinkedPair.from(chainStart, chainEnd);

                if (chain.linkWouldCloseChain(pair))
                {
                    chain.addLink(pair, true);
                }
            }
        }

        mChainFinder.clear();
        return chain;
    }

    public void reportCluster(final String sampleId, final SvCluster cluster)
    {
        // check whether this cluster has been processed
        analyseCluster(cluster, false);

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

        for(final SvVarData var : highPloidySVs)
        {
            ++typeCounts[typeAsInt(var.type())];

            if(minDMPloidy == 0 || var.ploidyMin() < minDMPloidy)
                minDMPloidy = var.ploidyMin();

            maxDMCopyNumber = max(maxDMCopyNumber, max(var.copyNumber(true), var.copyNumber(false)));

            if(!chromosomes.contains(var.chromosome(true)))
                chromosomes.add(var.chromosome(true));

            if(!var.isNullBreakend() && !chromosomes.contains(var.chromosome(false)))
                chromosomes.add(var.chromosome(false));

            svIds = appendStr(svIds, var.id(), ';');
        }

        final String dmTypesStr = getSvTypesStr(typeCounts);

        final PurityContext purityContext = mCnAnalyser.getPurityContext();
        double samplePurity = purityContext != null ? purityContext.bestFit().purity() : 0;
        double samplePloidy = purityContext != null ? purityContext.bestFit().ploidy() : 0;

        long dmChainLength = chain != null ? chain.getLength(false) : 0;
        int chainSvCount = chain != null ? chain.getSvCount() : 0;

        // get amplified genes list by looking at all section traversed by this chain or breakends with genes in them?
        String amplifiedGenesStr = chain != null ? getAmplifiedGenesList(chain) : "";

        long posStart = 0;
        long posEnd = 0;

        if(fullyChained && highPloidySVs.size() == 1)
        {
            SvVarData dup = highPloidySVs.get(0);
            posStart = dup.position(true);
            posEnd = dup.position(false);
        }

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

                outputFileName += "SVA_DM.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,ClusterId,ClusterDesc,ResolvedType,ClusterCount,SamplePurity,SamplePloidy,DMSvCount,DMSvTypes");
                mFileWriter.write(",FullyChained,ChainLength,ChainCount,SvIds,Chromosomes,DupPosStart,DupPosEnd");
                mFileWriter.write(",MaxCopyNumber,MinPloidy,AmpGenes");
                mFileWriter.newLine();
            }

            mFileWriter.write(String.format("%s,%d,%s,%s,%d",
                    sampleId, cluster.id(), cluster.getDesc(), cluster.getResolvedType(), cluster.getSvCount()));

            mFileWriter.write(String.format(",%.2f,%.2f,%d,%s,%s,%d,%d",
                    samplePurity, samplePloidy, highPloidySVs.size(), dmTypesStr, fullyChained, dmChainLength, chainSvCount));

            mFileWriter.write(String.format(",%s,%s,%d,%d,%.2f,%.2f,%s",
                    svIds, chromosomeStr, posStart, posEnd, maxDMCopyNumber, minDMPloidy, amplifiedGenesStr));

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

        /*
                List<Integer> ploidyBuckets = Lists.newArrayList();
            int maxPermittedBuckets = calcMaxBuckets(svCount, clusterMaxPloidy);

            for (final SvVarData var : cluster.getSVs())
            {
                int ploidyBucketMin = (int)round(var.ploidyMin());
                int ploidyBucketMax = (int)round(var.ploidyMax());

                boolean matched = false;
                for(Integer ploidyBucket : ploidyBuckets)
                {
                    if(ploidyBucketMin <= ploidyBucket && ploidyBucketMax >= ploidyBucket)
                    {
                        matched = true;
                        break;
                    }
                }

                if (!matched)
                {
                    int avgPloidyBucket = (int)round((ploidyBucketMin + ploidyBucketMax) * 0.5);
                    ploidyBuckets.add(avgPloidyBucket);
                }
            }

     */

    public void close()
    {
        closeBufferedWriter(mFileWriter);
    }

    private static final double PLOIDY_STEPWISE_FACTOR = 0.8;

    private static int calcMaxBuckets(int svCount, double maxPloidy)
    {
        double expectedBuckets = log(svCount) * log(maxPloidy) * PLOIDY_STEPWISE_FACTOR;
        return min(svCount,(int)round(expectedBuckets));
    }

}
