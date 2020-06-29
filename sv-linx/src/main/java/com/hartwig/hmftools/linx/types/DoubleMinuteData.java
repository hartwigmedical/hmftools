package com.hartwig.hmftools.linx.types;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.positionWithin;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.linx.chaining.SvChain;

public class DoubleMinuteData
{
    public final SvCluster Cluster;
    public final List<SvVarData> SVs;

    public double MaxBFBJcn;
    public double MinAdjacentMARatio;

    public final List<SvVarData> CandidateSVs;
    public final List<SvChain> Chains;
    public boolean FullyChained;

    // annotations relating to chained segments
    public long ClosedSegmentLength; // Sum of closed segment length
    public int ClosedBreakends; // # of closed breakends
    public double ClosedJcnTotal; // sum of JCN of closed breakends
    public int OpenBreakends; // # of open breakends
    public double OpenJcnTotal; // Sum of JCN of open breakends
    public double OpenJcnMax; // max JCN of open breakends

    public double TotalSegmentCnChange;
    public boolean ChainsCentromere;

    public Map<StructuralVariantType,double[]> InternalTypeData;

    public DoubleMinuteData(final SvCluster cluster, final List<SvVarData> svList)
    {
        Cluster = cluster;
        SVs = svList;

        MinAdjacentMARatio = 0;
        MaxBFBJcn = 0;

        Chains = Lists.newArrayList();
        CandidateSVs = Lists.newArrayList();
        FullyChained = false;

        ClosedSegmentLength = 0;
        ClosedBreakends = 0;
        ClosedJcnTotal = 0;
        OpenBreakends = 0;
        OpenJcnTotal = 0;
        OpenJcnMax = 0;
        TotalSegmentCnChange = 0;
        ChainsCentromere = false;

        InternalTypeData = Maps.newHashMap();
        InternalTypeData.put(INV, new double[INT_SEG_SUM+1]);
        InternalTypeData.put(SGL, new double[INT_SEG_SUM+1]);
        InternalTypeData.put(INF, new double[INT_SEG_SUM+1]);
    }

    public void setChainCharacteristics(final Map<String,List<SvBreakend>> chrBreakendMap)
    {
        final Set<SvBreakend> observedBreakends = Sets.newHashSet();
        for(SvChain chain : Chains)
        {
            for(SvLinkedPair pair : chain.getLinkedPairs())
            {
                ClosedSegmentLength += pair.length();
                ClosedBreakends += 2;
                ClosedJcnTotal += pair.first().jcn();
                ClosedJcnTotal += pair.second().jcn();

                if(pair.firstBreakend().arm() != pair.secondBreakend().arm())
                    ChainsCentromere = true;

                setPairCharacteristics(pair, chrBreakendMap, observedBreakends);

                // only count towards total if link ends aren't contained within another link
                TotalSegmentCnChange += getPairCnChange(pair);
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
    }

    private double getPairCnChange(final SvLinkedPair pair)
    {
        // checks if the pair has neither end within another chain link
        double cnChange = 0;
        for(int se = SE_START; se <= SE_END; ++se)
        {
            SvBreakend pairBE = pair.getBreakend(se);
            boolean isOuter = true;


            for(SvChain chain : Chains)
            {
                for(SvLinkedPair chainLink : chain.getLinkedPairs())
                {
                    if(chainLink.chromosome().equals(pair))
                    {
                        if(positionWithin(pairBE.position(), chainLink.getBreakend(true).position(), chainLink.getBreakend(false).position()))
                        {
                            isOuter = false;
                            break;
                        }
                    }
                }

                if(!isOuter)
                    break;
            }

            if(isOuter)
            {
                cnChange += pairBE.copyNumberLowSide() * (pairBE.orientation() == -1 ? 1 : -1);
            }
        }

        return cnChange;
    }

    public static final int INT_SEG_COUNT = 0;
    public static final int INT_SEG_MAX = 1;
    public static final int INT_SEG_SUM = 2;

    private void setPairCharacteristics(final SvLinkedPair pair, final Map<String,List<SvBreakend>> chrBreakendMap, final Set<SvBreakend> observedBreakends)
    {
        final List<SvBreakend> breakendList = chrBreakendMap.get(pair.chromosome());

        for(int index = pair.getBreakend(SE_START).getChrPosIndex() + 1; index < pair.getBreakend(SE_END).getChrPosIndex(); ++index)
        {
            final SvBreakend breakend = breakendList.get(index);

            if(observedBreakends.contains(breakend))
                continue;

            observedBreakends.add(breakend);

            StructuralVariantType svType = breakend.getSV().type();
            if(svType != SGL && svType != INF)
            {
                if(breakend.isFoldback())
                    svType = INV;
                else
                    continue;
            }

            double[] counts = InternalTypeData.get(svType);
            counts[INT_SEG_COUNT] += 1;
            counts[INT_SEG_SUM] += breakend.jcn();
            counts[INT_SEG_MAX] = max(counts[INT_SEG_MAX], breakend.jcn());
        }
    }

    public String internalTypeCountsAsStr()
    {
        return String.format("%.0f,%.1f,%.1f,%.0f,%.1f,%.1f,%.0f,%.1f,%.1f",
                InternalTypeData.get(INV)[INT_SEG_COUNT], InternalTypeData.get(INV)[INT_SEG_SUM], InternalTypeData.get(INV)[INT_SEG_MAX],
                InternalTypeData.get(SGL)[INT_SEG_COUNT], InternalTypeData.get(SGL)[INT_SEG_SUM], InternalTypeData.get(SGL)[INT_SEG_MAX],
                InternalTypeData.get(INF)[INT_SEG_COUNT], InternalTypeData.get(INF)[INT_SEG_SUM], InternalTypeData.get(INF)[INT_SEG_MAX]);
    }

}
