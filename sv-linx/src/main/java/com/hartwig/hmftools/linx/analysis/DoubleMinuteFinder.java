package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.typeAsInt;
import static com.hartwig.hmftools.linx.analysis.LinkFinder.areLinkedSection;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getSvTypesStr;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNONTATION_DM;
import static com.hartwig.hmftools.linx.types.SvCluster.isSpecificCluster;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LINK_TYPE_TI;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isSpecificSV;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;

import java.io.BufferedWriter;
import java.io.IOException;
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
    private ChainFinder mChainFinder;

    private String mOutputDir;
    private BufferedWriter mFileWriter;

    private static double PLOIDY_THRESHOLD = 3.5;
    private static double FOLDBACK_PLOIDY_RATIO = 0.5;
    private static double PLOIDY_STEPWISE_FACTOR = 0.8;
    private static double HIGH_PLOIDY_FACTOR = 1.1;
    private static double ADJACENT_PLOIDY_RATIO = 2.3;

    // old constants
    private static double DM_PLOIDY_MIN_RATIO = 2.3;
    private static double DM_MIN_PLOIDY = 3;
    private static double DM_PLOIDY_INCOMPLETE_MIN_RATIO = 4;
    private static double DM_INCOMPLETE_MIN_PLOIDY = 10;
    private static int DM_MAX_SV_COUNT = 16;

    private static final Logger LOGGER = LogManager.getLogger(DoubleMinuteFinder.class);

    public DoubleMinuteFinder()
    {
        mChainFinder = new ChainFinder();
        mCnAnalyser = null;
        mGeneTransCache = null;
        mOutputDir = null;
        mFileWriter = null;
    }

    public void setGeneTransCache(final SvGeneTranscriptCollection geneTransCache) { mGeneTransCache = geneTransCache; }
    public void setCopyNumberAnalyser(CnDataLoader cnAnalyser) { mCnAnalyser = cnAnalyser; }

    public void setOutputDir(final String outputDir)
    {
        mOutputDir = outputDir;
    }

    public void analyseCluster(final String sampleId, SvCluster cluster)
    {
        if(cluster.hasAnnotation(CLUSTER_ANNONTATION_DM))
            return;

        if(cluster.getSvCount() == 1 && cluster.getSV(0).type() != DUP)
            return;

        // isSpecificCluster(cluster);

        double clusterMaxPloidy = cluster.getSVs().stream().mapToDouble(x -> x.ploidy()).max().getAsDouble();

        // clusterMaxPloidy = cluster.getMaxPloidy(); - no longer used since driven by replication logic

        if(clusterMaxPloidy < PLOIDY_THRESHOLD)
            return;

        // dismiss the cluster if explained by foldbacks
        double foldbackPloidyTotal = cluster.getFoldbacks().stream().mapToDouble(x -> x.ploidy()).sum();

        if(foldbackPloidyTotal >= FOLDBACK_PLOIDY_RATIO * clusterMaxPloidy)
        {
            LOGGER.debug(String.format("cluster(%s) maxPloidy(%.1f) foldbacks(count=%d ploidyTotal=%.1f) invalidates DM",
                    cluster.id(), clusterMaxPloidy, cluster.getFoldbacks().size(), foldbackPloidyTotal));
            return;
        }

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

                if (svMaxPloidy >= clusterMaxPloidy)
                    highPloidySVs.add(var);
            }
        }

        // at least one high-ploidy breakend must have a high ploidy relative to the
        // adjacent CN segment's major allele ploidy
        int index = 0;
        while(index < highPloidySVs.size())
        {
            SvVarData var = highPloidySVs.get(index);
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

        if(highPloidySVs.isEmpty())
            return;

        if(highPloidySVs.size() == 1 && highPloidySVs.get(0).type() != DUP)
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

        if(fullyChained)
            cluster.setDoubleMinuteSVs(highPloidySVs);

        cluster.addAnnotation(CLUSTER_ANNONTATION_DM);

        if(highPloidySVs.size() == 1 && fullyChained)
        {
            // single DUPs won't go through the chaining routine so cache this chain here
            cluster.addChain(dmChain, false);
        }

        LOGGER.debug(String.format("cluster(%s) maxPloidy(%.1f) dmSvCount(%d) fullyChained(%s)",
                cluster.id(), clusterMaxPloidy, highPloidySVs.size(), fullyChained));

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
        mChainFinder.formClusterChains(false);

        if(mChainFinder.getChains().size() != 1)
            return null;

        SvChain chain = mChainFinder.getChains().get(0);

        // check whether the chain could form a loop
        SvBreakend chainStart = chain.getOpenBreakend(true);
        SvBreakend chainEnd = chain.getOpenBreakend(false);

        if(chainStart != null && !chainStart.getSV().isNullBreakend() && chainEnd != null && !chainEnd.getSV().isNullBreakend())
        {
            if (areLinkedSection(chainStart.getSV(), chainEnd.getSV(), chainStart.usesStart(), chainEnd.usesStart(), false))
            {
                SvLinkedPair pair = SvLinkedPair.from(chainStart, chainEnd, LINK_TYPE_TI);

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
        if(mOutputDir.isEmpty())
            return;

        List<SvVarData> highPloidySVs = cluster.getDoubleMinuteSVs();

        if(highPloidySVs.isEmpty())
            return;

        final SvChain chain = createDMChain(cluster, highPloidySVs);

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

    private static int calcMaxBuckets(int svCount, double maxPloidy)
    {
        double expectedBuckets = log(svCount) * log(maxPloidy) * PLOIDY_STEPWISE_FACTOR;
        return min(svCount,(int)round(expectedBuckets));
    }

}
