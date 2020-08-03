package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionsOverlap;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.DoubleMinuteFinder.getAdjacentMajorAPRatio;
import static com.hartwig.hmftools.linx.types.LinxConstants.ADJACENT_JCN_RATIO;
import static com.hartwig.hmftools.linx.types.SvVarData.RELATION_TYPE_NEIGHBOUR;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class DoubleMinuteData
{
    public final SvCluster Cluster;
    public final List<SvVarData> SVs; // identified double minute SVs
    public final List<SvVarData> UnchainedSVs; // subset which could not be chained

    public double MaxJcn;
    public double MinAdjacentMARatio;
    public double MinAdjMAJcnRatio;

    public final List<SvVarData> CandidateSVs;
    public final List<SvChain> Chains;
    public final List<SvChain> ValidChains;
    public boolean FullyChained;

    // annotations relating to chained segments
    public long ClosedSegmentLength; // Sum of closed segment length
    public int ClosedBreakends; // # of closed breakends
    public double ClosedJcnTotal; // sum of JCN of closed breakends
    public int OpenBreakends; // # of open breakends
    public double OpenJcnTotal; // Sum of JCN of open breakends
    public double OpenJcnMax; // max JCN of open breakends
    public int SimpleDels; // # of non-overlapping DELs

    public int NonSegmentFoldbacks;
    public double NonSegmentFoldbackJcnTotal;

    public double TotalSegmentCnChange;
    public boolean ChainsCentromere;

    public double[] FbInternalData;
    public Map<Integer,double[]> ChainIntExtData; // data for SVs going from a chained segment to outside any chained segment
    public Map<Integer,double[]> ChainSglInternalData;
    public Map<Integer,double[]> ChainInfInternalData;

    private static final int OPEN_CHAIN_ID = -1;

    private boolean mIsDoubleMinute;

    public DoubleMinuteData(final SvCluster cluster, final List<SvVarData> svList)
    {
        Cluster = cluster;
        SVs = svList;
        UnchainedSVs = Lists.newArrayList();

        MinAdjacentMARatio = 0;
        MaxJcn = 0;

        Chains = Lists.newArrayList();
        ValidChains = Lists.newArrayList();
        CandidateSVs = Lists.newArrayList();
        FullyChained = false;

        MinAdjMAJcnRatio = 0;

        ClosedSegmentLength = 0;
        ClosedBreakends = 0;
        ClosedJcnTotal = 0;
        OpenBreakends = 0;
        OpenJcnTotal = 0;
        OpenJcnMax = 0;
        SimpleDels = 0;
        TotalSegmentCnChange = 0;
        ChainsCentromere = false;

        FbInternalData = new double[SEG_DATA_SUM +1];
        ChainIntExtData = Maps.newHashMap();
        ChainSglInternalData = Maps.newHashMap();
        ChainInfInternalData = Maps.newHashMap();

        mIsDoubleMinute = false;
    }

    public boolean isDoubleMinute() { return mIsDoubleMinute; }

    public void annotate(final Map<String,List<SvBreakend>> chrBreakendMap)
    {
        for(SvVarData var : SVs)
        {
            if(!Chains.stream().anyMatch(x -> x.getSvList().contains(var)))
            {
                UnchainedSVs.add(var);
            }

            if(MinAdjMAJcnRatio == 0)
                MinAdjMAJcnRatio = getAdjacentMajorAPRatio(var);
            else
                MinAdjMAJcnRatio = min(getAdjacentMajorAPRatio(var), MinAdjMAJcnRatio);

            MaxJcn = max(MaxJcn, var.jcn());
        }

        final List<LinkedPair> allLinkedPairs = Lists.newArrayList();
        Chains.forEach(x -> allLinkedPairs.addAll(x.getLinkedPairs()));

        setChainCharacteristics(chrBreakendMap, allLinkedPairs);

        for(SvVarData var : Cluster.getFoldbacks())
        {
            if(var.isChainedFoldback())
                continue;

            if(allLinkedPairs.stream()
                    .anyMatch(x -> x.chromosome().equals(var.chromosome(true))
                            && positionsOverlap(
                                    var.position(true), var.position(false),
                                    x.getBreakend(SE_START).position(), x.getBreakend(SE_END).position())))
            {
                continue;
            }

            ++NonSegmentFoldbacks;
            NonSegmentFoldbackJcnTotal += var.jcn();
        }

        mIsDoubleMinute = checkCriteria();
    }

    private void setChainCharacteristics(final Map<String,List<SvBreakend>> chrBreakendMap, final List<LinkedPair> allLinkedPairs)
    {
        final Set<SvBreakend> observedBreakends = Sets.newHashSet();
        final Set<SvVarData> observedSVs = Sets.newHashSet();

        for(SvChain chain : Chains)
        {
            // simple DELs don't contribute towards the closed breakend count
            final List<SvVarData> simpleDels = chain.getSvList().stream()
                    .filter(x -> x.type() == DEL)
                    .filter(x -> x.getNearestSvRelation().equals(RELATION_TYPE_NEIGHBOUR))
                    .collect(Collectors.toList());

            SimpleDels += simpleDels.size();

            for(LinkedPair pair : chain.getLinkedPairs())
            {
                ClosedSegmentLength += pair.length();

                if(!simpleDels.contains(pair.first()))
                {
                    ++ClosedBreakends;
                    ClosedJcnTotal += pair.first().jcn();
                }

                if(!simpleDels.contains(pair.second()))
                {
                    ++ClosedBreakends;
                    ClosedJcnTotal += pair.second().jcn();
                }

                if(pair.firstBreakend().arm() != pair.secondBreakend().arm())
                    ChainsCentromere = true;

                setPairCharacteristics(chain, pair, chrBreakendMap, allLinkedPairs, observedBreakends, observedSVs);
            }

            if(!chain.isClosedLoop())
            {
                if(!chain.getFirstSV().isSglBreakend())
                {
                    OpenJcnTotal += chain.getFirstSV().jcn();
                    ++OpenBreakends;
                    OpenJcnMax = max(OpenJcnMax, chain.getFirstSV().jcn());
                }

                if(!chain.getLastSV().isSglBreakend())
                {
                    OpenJcnTotal += chain.getLastSV().jcn();
                    ++OpenBreakends;
                    OpenJcnMax = max(OpenJcnMax, chain.getLastSV().jcn());
                }
            }
        }

        // also account for SVs not in chains
        for(SvVarData var : UnchainedSVs)
        {
            OpenJcnTotal += var.jcn() * 2;
            OpenBreakends += 2;
            OpenJcnMax = max(OpenJcnMax, var.jcn());
        }
    }

    public static final int SEG_DATA_COUNT = 0;
    public static final int SEG_DATA_MAX = 1;
    public static final int SEG_DATA_SUM = 2;

    private void setPairCharacteristics(
            final SvChain chain, final LinkedPair pair, final Map<String,List<SvBreakend>> chrBreakendMap, final List<LinkedPair> allLinkedPairs,
            final Set<SvBreakend> observedBreakends, final Set<SvVarData> observedSVs)
    {
        final List<SvBreakend> breakendList = chrBreakendMap.get(pair.chromosome());

        // If 2 variants are assembled from Internal-External-Internal or External-Internal-External (ie short templated insertions),
        // exclude from the Internal-External count

        int startIndex = pair.getBreakend(SE_START).getChrPosIndex();
        int endIndex = pair.getBreakend(SE_END).getChrPosIndex();

        final Set<SvBreakend> assembledBreakends = Sets.newHashSet();

        for(int index = startIndex + 1; index < endIndex; ++index)
        {
            final SvBreakend breakend = breakendList.get(index);

            if(SVs.contains(breakend.getSV()))
                continue;

            if(observedBreakends.contains(breakend))
                continue;

            observedBreakends.add(breakend);

            if(observedSVs.contains(breakend.getSV()))
                continue;

            observedSVs.add(breakend.getSV());

            if(!breakend.getSV().isSglBreakend())
            {
                // look for breakends going from within a segment to outside all segments
                final SvBreakend otherBreakend = breakend.getOtherBreakend();

                boolean otherBreakendInSegment = allLinkedPairs.stream()
                        .anyMatch(x -> x.chromosome().equals(otherBreakend.chromosome())
                                && positionWithin(
                                otherBreakend.position(), x.getBreakend(SE_START).position(), x.getBreakend(SE_END).position()));

                if(!otherBreakendInSegment)
                {
                    //  check for the other breakend forming an assembled TI with the next breakend
                    boolean inShortLink = assembledBreakends.contains(breakend) ||
                            (index < endIndex - 1 ? inShortAssembledLink(breakend, breakendList.get(index + 1)) : false);

                    if(inShortLink)
                    {
                        assembledBreakends.add(breakend);
                        assembledBreakends.add(breakendList.get(index + 1));
                    }
                    else
                    {
                        LNX_LOGGER.debug("cluster({}) pair({}) has in-out SV({})", Cluster.id(), pair.toString(), breakend.getSV());

                        double[] segmentData = getOrAddSegmentData(chain, ChainIntExtData);

                        segmentData[SEG_DATA_COUNT] += 1;
                        segmentData[SEG_DATA_SUM] += breakend.jcn();
                        segmentData[SEG_DATA_MAX] = max(segmentData[SEG_DATA_MAX], breakend.jcn());
                    }
                }
            }

            StructuralVariantType svType = breakend.type();
            if(svType != SGL && svType != INF)
            {
                if(breakend.isFoldback() && !breakend.getSV().isChainedFoldback())
                    svType = INV;
                else
                    continue;
            }

            double[] segmentData;

            if(svType == INV)
            {
                segmentData = FbInternalData;
            }
            else
            {
                Map<Integer,double[]> chainMap = svType == SGL ? ChainSglInternalData : ChainInfInternalData;
                segmentData = getOrAddSegmentData(chain, chainMap);
            }

            segmentData[SEG_DATA_COUNT] += 1;
            segmentData[SEG_DATA_SUM] += breakend.jcn();
            segmentData[SEG_DATA_MAX] = max(segmentData[SEG_DATA_MAX], breakend.jcn());
        }
    }

    private double[] getOrAddSegmentData(final SvChain chain, Map<Integer,double[]> segmentMap)
    {
        int chainId = chain.isClosedLoop() ? chain.id() : OPEN_CHAIN_ID;

        double[] segmentData = segmentMap.get(chainId);

        if(segmentData == null)
        {
            segmentData = new double[SEG_DATA_SUM + 1];
            segmentMap.put(chainId, segmentData);
        }

        return segmentData;
    }

    private boolean inShortAssembledLink(final SvBreakend breakend, final SvBreakend nextBreakend)
    {
        if(breakend.getLinkedPairs().stream().filter(x -> x.isAssembled()).anyMatch(x -> x.hasBreakend(nextBreakend)))
        {
            return true;
        }
        else
        {
            final SvBreakend otherBreakend = breakend.getOtherBreakend();
            final SvBreakend nextOtherBreakend = nextBreakend.getOtherBreakend();

            if(otherBreakend == null || nextOtherBreakend == null)
                return false;

            if(otherBreakend.getLinkedPairs().stream().filter(x -> x.isAssembled()).anyMatch(x -> x.hasBreakend(nextOtherBreakend)))
            {
                return true;
            }
        }

        return false;
    }

    private static final int MIN_CLOSED_SEG_LENGTH = 1500;
    private static final double MIN_CLOSED_SEG_RATIO = 0.66;
    private static final double LOW_JCN_BUFFER = 4;

    private boolean checkCriteria()
    {
        /* The following criteria must be met:
        - The total length of closed segments > 1500 bases
        -  Either a complete closed chain is formed OR
            - Closed breakends / (OpenBreakends + ClosedBreakends >= 2/3
        - If only a single DUP, then both sides must meet DM criteria
        - (Max JCN - 4 )  / [sum(intExtJCN excluding assembled TI) + max(max_INT_INF_JCN,max_INT_SGL_JCN) + sum(FB JN in cluster)] > 1
            For closed chains all INT JCN are for breakends on that closed chain only
            FB are always at the full cluster level
        */

        if(ClosedSegmentLength < MIN_CLOSED_SEG_LENGTH)
            return false;

        double closedSegmentRatio = ClosedBreakends / (double)(ClosedBreakends + OpenBreakends);

        if(closedSegmentRatio < MIN_CLOSED_SEG_RATIO)
            return false;

        if(SVs.size() == 1)
        {
            final SvVarData var = SVs.get(0);
            if(var.type() != DUP)
                return false;

            double maxAdjacentMaJcn = max(
                    var.getBreakend(true).majorAlleleJcn(true),
                    var.getBreakend(false).majorAlleleJcn(false));

            if(var.jcn() / max(maxAdjacentMaJcn, 0.01) < ADJACENT_JCN_RATIO)
                return false;
        }

        double foldbackJcn = Cluster.getFoldbacks().stream()
                .filter(x -> !SVs.contains(x))
                .mapToDouble(x -> x.isChainedFoldback() ? x.jcn() * 0.5 : x.jcn())
                .sum();

        boolean hasValidCriteria = false;

        //  check closed chains and open chains / SVs as well
        final List<Integer> chainIds = Lists.newArrayList(OPEN_CHAIN_ID);
        Chains.stream().filter(x -> x.isClosedLoop()).forEach(x -> chainIds.add(x.id()));

        final List<SvVarData> nonClosedChainSVs = Lists.newArrayList(UnchainedSVs);
        Chains.stream().filter(x -> !x.isClosedLoop()).forEach(x -> nonClosedChainSVs.addAll(x.getSvList()));

        for(Integer chainId : chainIds)
        {
            final SvChain chain = Chains.stream().filter(x -> x.id() == chainId).findFirst().orElse(null);

            if(nonClosedChainSVs.isEmpty() && chainId == OPEN_CHAIN_ID)
                continue;

            // find the highest foldback JCN or if none the highest JCN
            double maxChainJcn = 0;
            double maxChainFoldbackJcn = 0;

            if(chain != null)
            {
                for(SvVarData var : chain.getSvList())
                {
                    maxChainJcn = max(maxChainJcn, var.jcn());

                    if(var.isFoldback())
                        maxChainFoldbackJcn = maxChainFoldbackJcn == 0 ? var.jcn() : min(maxChainFoldbackJcn, var.jcn());
                }
            }
            else
            {
                for(SvVarData var : nonClosedChainSVs)
                {
                    maxChainJcn = max(maxChainJcn, var.jcn());

                    if(var.isFoldback())
                        maxChainFoldbackJcn = maxChainFoldbackJcn == 0 ? var.jcn() : min(maxChainFoldbackJcn, var.jcn());
                }
            }

            double maxJcn = maxChainFoldbackJcn > 0 ? maxChainFoldbackJcn : maxChainJcn;

            final double[] sglInternalValues = ChainSglInternalData.get(chainId);
            final double[] infInternalValues = ChainSglInternalData.get(chainId);
            final double[] intExtValues = ChainIntExtData.get(chainId);

            double sglInfMaxJcn = max(
                    sglInternalValues != null ? sglInternalValues[SEG_DATA_MAX] : 0,
                    infInternalValues != null ? infInternalValues[SEG_DATA_MAX] : 0);

            double intExtJcnTotal = intExtValues != null ? intExtValues[SEG_DATA_SUM] : 0;

            double opposingJcn = intExtJcnTotal + sglInfMaxJcn + foldbackJcn;

            if(opposingJcn + LOW_JCN_BUFFER <= maxJcn)
            {
                hasValidCriteria = true;

                if(chain != null)
                    ValidChains.add(chain);
            }
        }

        return hasValidCriteria;
    }

    public String internalTypeCountsAsStr()
    {
        return String.format("%.0f,%.1f,%.1f,%.0f,%.1f,%.1f,%.0f,%.1f,%.1f,%.0f,%.1f,%.1f",
                ChainIntExtData.values().stream().mapToDouble(x -> x[SEG_DATA_COUNT]).sum(),
                ChainIntExtData.values().stream().mapToDouble(x -> x[SEG_DATA_SUM]).sum(),
                ChainIntExtData.values().stream().mapToDouble(x -> x[SEG_DATA_MAX]).max().orElse(0),
                FbInternalData[SEG_DATA_COUNT], FbInternalData[SEG_DATA_SUM], FbInternalData[SEG_DATA_MAX],
                ChainSglInternalData.values().stream().mapToDouble(x -> x[SEG_DATA_COUNT]).sum(),
                ChainSglInternalData.values().stream().mapToDouble(x -> x[SEG_DATA_SUM]).sum(),
                ChainSglInternalData.values().stream().mapToDouble(x -> x[SEG_DATA_MAX]).max().orElse(0),
                ChainInfInternalData.values().stream().mapToDouble(x -> x[SEG_DATA_COUNT]).sum(),
                ChainInfInternalData.values().stream().mapToDouble(x -> x[SEG_DATA_SUM]).sum(),
                ChainInfInternalData.values().stream().mapToDouble(x -> x[SEG_DATA_MAX]).max().orElse(0));
    }

}
