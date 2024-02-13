package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.classifySinglePairResolvedType;
import static com.hartwig.hmftools.linx.analysis.DoubleMinuteFinder.JCN_UPPER_THRESHOLD;
import static com.hartwig.hmftools.linx.analysis.DoubleMinuteFinder.MIN_SEGMENT_DEPTH_WINDOW_COUNT;
import static com.hartwig.hmftools.linx.analysis.DoubleMinuteFinder.getAdjacentMajorAPRatio;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.types.LinxConstants.ADJACENT_JCN_RATIO;
import static com.hartwig.hmftools.linx.types.SvVarData.RELATION_TYPE_NEIGHBOUR;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
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
    public double MinAdjMAJcnRatio;

    public final List<SvVarData> CandidateSVs;
    public final List<SvChain> Chains;
    public boolean FullyChained;

    public final List<SvChain> ValidChains;
    public final List<SvVarData> ValidSVs; // those in valid chains

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
    public boolean LowDepthCount;

    public double[] FbInternalData;
    public Map<Integer,double[]> ChainIntExtData; // data for SVs going from a chained segment to outside any chained segment
    public Map<Integer,double[]> ChainSglInternalData;
    public Map<Integer,double[]> ChainInfInternalData;

    private final List<ChrBaseRegion> mDmRegions;

    private static final int OPEN_CHAIN_ID = -1;

    private boolean mIsDoubleMinute;

    public DoubleMinuteData(final SvCluster cluster, final List<SvVarData> svList)
    {
        Cluster = cluster;
        SVs = svList;
        UnchainedSVs = Lists.newArrayList();

        MaxJcn = 0;

        Chains = Lists.newArrayList();
        ValidChains = Lists.newArrayList();
        ValidSVs = Lists.newArrayList();
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
        LowDepthCount = false;

        FbInternalData = new double[SEG_DATA_SUM +1];
        ChainIntExtData = Maps.newHashMap();
        ChainSglInternalData = Maps.newHashMap();
        ChainInfInternalData = Maps.newHashMap();

        mDmRegions = Lists.newArrayList();
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

        buildDmRegions();

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

        mIsDoubleMinute = checkCriteria(allLinkedPairs);
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
                ClosedSegmentLength += pair.baseLength();

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
                        && positionWithin(otherBreakend.position(), x.getBreakend(SE_START).position(), x.getBreakend(SE_END).position()));

                if(!otherBreakendInSegment)
                {
                    otherBreakendInSegment =
                            mDmRegions.stream().anyMatch(x -> x.containsPosition(otherBreakend.chromosome(), otherBreakend.position()));

                    // breakend must not be within 1K of any DM SV
                    // otherBreakendInSegment = isCloseToDmSV(breakend);
                }

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
                        LNX_LOGGER.debug("cluster({}) pair({}) has in-out SV({}) with JCN({})",
                                Cluster.id(), pair.toString(), breakend.getSV(), String.format("%.1f", breakend.jcn()));

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
                // ignore SGLs in a pair
                if(svType == SGL && index < endIndex - 1)
                {
                    final SvBreakend nextBreakend = breakendList.get(index + 1);
                    if(inSglPairCluster(breakend, nextBreakend))
                    {
                        observedBreakends.add(nextBreakend);
                        observedSVs.add(nextBreakend.getSV());
                        continue;
                    }
                }

                Map<Integer,double[]> chainMap = svType == SGL ? ChainSglInternalData : ChainInfInternalData;
                segmentData = getOrAddSegmentData(chain, chainMap);
            }

            segmentData[SEG_DATA_COUNT] += 1;
            segmentData[SEG_DATA_SUM] += breakend.jcn();
            segmentData[SEG_DATA_MAX] = max(segmentData[SEG_DATA_MAX], breakend.jcn());
        }
    }

    private static boolean inSglPairCluster(final SvBreakend breakend, final SvBreakend nextBreakend)
    {
        if(nextBreakend.type() != SGL || nextBreakend.getCluster() != breakend.getCluster())
            return false;

        return breakend.getCluster().getResolvedType().isSimple()
                || classifySinglePairResolvedType(breakend.getSV(), nextBreakend.getSV()).isSimple();
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

    private void buildDmRegions()
    {
        if(SVs.size() == 1)
        {
            final SvVarData var = SVs.get(0);

            if(var.type() == DUP)
            {
                mDmRegions.add(new ChrBaseRegion(var.chromosome(true), var.position(true), var.position(false)));
            }

            return;
        }

        for(Map.Entry<String,List<SvBreakend>> entry : Cluster.getChrBreakendMap().entrySet())
        {
            final String chromosome = entry.getKey();

            if(SVs.stream().noneMatch(x -> x.chromosome(true).equals(chromosome) || x.chromosome(false).equals(chromosome)))
                continue;

            ChrBaseRegion currentRegion = null;
            double regionStartCopyNumber = 0;

            for(SvBreakend breakend : entry.getValue())
            {
                if(!SVs.contains(breakend.getSV()))
                    continue;

                if(breakend.orientation() == NEG_ORIENT)
                {
                    if(currentRegion != null)
                        continue;

                    if(getMajorAlleleJcnRatio(breakend) >= ADJACENT_JCN_RATIO)
                    {
                        regionStartCopyNumber = breakend.copyNumberLowSide();
                        currentRegion = new ChrBaseRegion(chromosome, breakend.position(), 0);
                    }
                }
                else
                {
                    if(currentRegion == null)
                    {
                        // last region was ended but now extend it onto this additional closing breakend
                        if(!mDmRegions.isEmpty())
                        {
                            ChrBaseRegion lastRegion = mDmRegions.get(mDmRegions.size() - 1);
                            lastRegion.setEnd(breakend.position());
                        }

                        continue;
                    }

                    double lowSideCn = breakend.copyNumberLowSide();

                    if(lowSideCn < regionStartCopyNumber || copyNumbersEqual(regionStartCopyNumber, lowSideCn)
                    || getMajorAlleleJcnRatio(breakend) >= ADJACENT_JCN_RATIO)
                    {
                        regionStartCopyNumber = 0;
                        currentRegion.setEnd(breakend.position());
                        mDmRegions.add(currentRegion);
                        currentRegion = null;
                    }
                }
            }
        }

        if(LNX_LOGGER.isDebugEnabled())
        {
            mDmRegions.forEach(x -> LNX_LOGGER.debug("cluster({}) dmRegion({})", Cluster.id(), x));
        }
    }

    private static final int MIN_CLOSED_SEG_LENGTH = 1500;
    private static final double MIN_CLOSED_SEG_RATIO = 0.66;
    private static final double LOW_JCN_BUFFER = 4;

    private boolean checkCriteria(final List<LinkedPair> allLinkedPairs)
    {
        /* The following criteria must be met:
        - The total length of closed segments > 1500 bases
        - Either a complete closed chain is formed OR
            - Closed breakends / (OpenBreakends + ClosedBreakends >= 2/3
        - If only a single DUP, then both sides must meet DM criteria
        - (Max JCN - 4 )  / [sum(intExtJCN excluding assembled TI) + max(max_INT_INF_JCN,max_INT_SGL_JCN) + sum(FB JN in cluster)] > 1
            For closed chains all INT JCN are for breakends on that closed chain only
            FB are always at the full cluster level
         - If the cluster includes only closed segments enclosed by single or inferred breakends on both sides,
            then at least one closed segment must have PURPLE depthWindowCount > 5
        */

        double closedSegmentRatio = ClosedBreakends / (double)(ClosedBreakends + OpenBreakends);

        if(closedSegmentRatio < MIN_CLOSED_SEG_RATIO)
            return false;

        if(SVs.size() == 1)
        {
            final SvVarData var = SVs.get(0);
            if(var.type() != DUP)
                return false;

            if(!variantExceedsBothAdjacentJcn(var))
                return false;
        }

        if(MaxJcn < JCN_UPPER_THRESHOLD)
        {
            if(SVs.stream().anyMatch(x -> x.isSglBreakend()))
                return false;

            if(Chains.isEmpty() || Chains.stream().anyMatch(x -> !x.isClosedLoop()))
                return false;

            // require 1 SV to satisfy the adjacent major allele JCN rule on both sides
            if(SVs.size() > 1)
            {
                if(SVs.stream().noneMatch(x -> variantExceedsBothAdjacentJcn(x)))
                    return false;
            }
        }

        if(allLinkedPairs.stream().allMatch(x -> haveLowAdjacentDepthWindowCount(x)))
        {
            LowDepthCount = true;
            return false;
        }

        boolean hasValidCriteria = setValidChains();

        return hasValidCriteria;
    }

    private static boolean haveLowAdjacentDepthWindowCount(final LinkedPair pair)
    {
        if(!pair.first().isSglBreakend() || !pair.second().isSglBreakend())
            return false;

        if(abs(pair.firstBreakend().getChrPosIndex() - pair.secondBreakend().getChrPosIndex()) != 1)
            return false;

        final SvBreakend first = pair.firstBreakend();
        final SvBreakend second = pair.secondBreakend();

        int depthWindowCount1 = first.getSV().getCopyNumberData(true, first.orientation() == 1).DepthWindowCount;
        int depthWindowCount2 = second.getSV().getCopyNumberData(true, second.orientation() == 1).DepthWindowCount;

        return min(depthWindowCount1, depthWindowCount2) < MIN_SEGMENT_DEPTH_WINDOW_COUNT;
    }

    protected static boolean variantExceedsBothAdjacentJcn(final SvVarData var)
    {
        if(var.isSglBreakend())
            return false;

        double minRatio =
                min(getMajorAlleleJcnRatio(var.getBreakend(true)), getMajorAlleleJcnRatio(var.getBreakend(false)));

        return minRatio >= ADJACENT_JCN_RATIO;
    }

    protected static double getMajorAlleleJcnRatio(final SvBreakend breakend)
    {
        return breakend.jcn() / max(breakend.majorAlleleJcn(breakend.orientation() == NEG_ORIENT), 0.01);
    }

    private boolean setValidChains()
    {
        boolean hasValidCriteria = false;

        double foldbackJcn = Cluster.getFoldbacks().stream()
                .filter(x -> !SVs.contains(x))
                .mapToDouble(x -> x.isChainedFoldback() ? x.jcn() * 0.5 : x.jcn())
                .sum();

        //  check closed chains and open chains / SVs as well
        final List<Integer> chainIds = Lists.newArrayList(OPEN_CHAIN_ID);
        Chains.stream().filter(x -> x.isClosedLoop()).forEach(x -> chainIds.add(x.id()));

        final List<SvVarData> nonClosedChainSVs = Lists.newArrayList(UnchainedSVs);
        final List<SvChain> openChains = Chains.stream().filter(x -> !x.isClosedLoop()).collect(Collectors.toList());
        openChains.forEach(x -> nonClosedChainSVs.addAll(x.getSvList()));

        long nonClosedSegmentLength = Chains.stream().filter(x -> !x.isClosedLoop()).mapToLong(x -> x.getLength(false)).sum();

        for(Integer chainId : chainIds)
        {
            final SvChain chain = Chains.stream().filter(x -> x.id() == chainId).findFirst().orElse(null);

            if(nonClosedChainSVs.isEmpty() && chainId == OPEN_CHAIN_ID)
                continue;

            if(chain != null && chain.getLength(false) < MIN_CLOSED_SEG_LENGTH)
                continue;
            else if(chain == null && nonClosedSegmentLength < MIN_CLOSED_SEG_LENGTH)
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
            final double[] infInternalValues = ChainInfInternalData.get(chainId);
            final double[] intExtValues = ChainIntExtData.get(chainId);

            double sglInfMaxJcn = max(
                    sglInternalValues != null ? sglInternalValues[SEG_DATA_MAX] : 0,
                    infInternalValues != null ? infInternalValues[SEG_DATA_MAX] : 0);

            double intExtJcnTotal = intExtValues != null ? intExtValues[SEG_DATA_SUM] : 0;

            double opposingJcn = intExtJcnTotal + sglInfMaxJcn + foldbackJcn;

            LNX_LOGGER.debug(String.format("cluster(%s) DM chain(%s) maxJcn(%.1f) opposingJcn(%.1f intExtJcnTotal=%.1f sglInfMaxJcn=%.1f foldbackJcn=%.1f)",
                    Cluster.id(), chain != null ? String.valueOf(chain.id()) : "open",
                    maxJcn, opposingJcn, intExtJcnTotal, sglInfMaxJcn, foldbackJcn));

            if(opposingJcn + LOW_JCN_BUFFER <= maxJcn)
            {
                hasValidCriteria = true;

                if(chain != null)
                {
                    ValidChains.add(chain);
                    ValidSVs.addAll(chain.getSvList());
                }
                else
                {
                    ValidSVs.addAll(nonClosedChainSVs);
                    ValidChains.addAll(openChains);
                }
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
