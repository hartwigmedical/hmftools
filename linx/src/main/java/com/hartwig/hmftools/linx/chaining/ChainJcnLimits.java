package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.cn.JcnCalcData;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

public class ChainJcnLimits
{
    private Map<String,List<SegmentJcn>> mChrAlleleJCNs;
    private double mValidAlleleJcnSegmentPerc;

    // referenced state
    private Map<String, List<SvBreakend>> mChrBreakendMap;
    private int mClusterId;

    // JCN level below which a chain segment cannot cross
    public static final double CLUSTER_ALLELE_JCN_MIN = 0.15;

    public ChainJcnLimits()
    {
        mChrAlleleJCNs = Maps.newHashMap();
        mClusterId = -1;
        mChrBreakendMap = null;
    }

    public final Map<String, List<SegmentJcn>> getChrAlleleJCNs() { return mChrAlleleJCNs; }
    public double getValidAlleleJcnSegmentPerc() { return mValidAlleleJcnSegmentPerc; }

    public void initialise(int clusterId, final Map<String,List<SvBreakend>> chrBreakendMap)
    {
        mChrBreakendMap = chrBreakendMap;
        mClusterId = clusterId;

        mChrAlleleJCNs.clear();
        mValidAlleleJcnSegmentPerc = 0;
    }

    public static final int RANGE_TOTAL = 0;
    public static final int DELETED_TOTAL = 1;

    public final long[] calcRangeData()
    {
        if(mChrAlleleJCNs.isEmpty())
            return null;

        long[] rangeData = {0, 0};

        for (final List<SegmentJcn> segmentList : mChrAlleleJCNs.values())
        {
            for(final SegmentJcn segment : segmentList)
            {
                if(segment.clusterJcn() >= CLUSTER_ALLELE_JCN_MIN)
                    rangeData[RANGE_TOTAL] += segment.length();
                else
                    rangeData[DELETED_TOTAL] += segment.length();
            }
        }

        return rangeData;
    }

    public void determineBreakendJCNs()
    {
        int totalSegCount = 0;
        int totalValidSegCount = 0;

        for (final Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            final String chromosome = entry.getKey();
            final List<SvBreakend> breakendList = entry.getValue();
            int breakendCount = breakendList.size();

            // a multi-dim array of breakend index for this arm to A allele JCN, B non-disrupted JCN, and cluster JCN
            List<SegmentJcn> alleleJCNs = Lists.newArrayListWithCapacity(breakendCount);

            mChrAlleleJCNs.put(chromosome, alleleJCNs);

            boolean inSegment = false;
            SvBreakend segStartBreakend = null;
            int segStartIndex = 0;

            for (int i = 0; i < breakendList.size(); ++i)
            {
                SvBreakend breakend = breakendList.get(i);
                SvBreakend nextBreakend = i < breakendList.size() - 1 ? breakendList.get(i + 1) : null;

                alleleJCNs.add(
                        new SegmentJcn(breakend.position(), nextBreakend != null ? nextBreakend.position() : breakend.position(),
                                breakend.majorAlleleJcn(false), breakend.minorAlleleJcn(false)));

                if(!inSegment)
                {
                    inSegment = true;
                    segStartBreakend = breakend;
                    segStartIndex = i;
                }
                else
                {
                    if(nextBreakend == null || nextBreakend.getChrPosIndex() > breakend.getChrPosIndex() + 1)
                    {
                        // a gap in the cluster so need to evaluate this contiguous section
                        inSegment = false;

                        boolean validCalcs = calculateNoneClusterSegmentPloidies(
                                alleleJCNs, segStartIndex, i, segStartBreakend, breakend);
                        ++totalSegCount;

                        if(validCalcs)
                            ++totalValidSegCount;
                    }
                }
            }

            calculateClusterSegmentJCNs(alleleJCNs);
        }

        LNX_LOGGER.debug("cluster({}) chromosomes({}) AP totalSegments({} valid={})",
                mClusterId, mChrBreakendMap.size(), totalSegCount, totalValidSegCount);
    }

    private boolean calculateNoneClusterSegmentPloidies(
            final List<SegmentJcn> alleleJCNs, int startIndex, int endIndex, SvBreakend startBreakend, SvBreakend endBreakend)
    {
        double startMajorAP = startBreakend.majorAlleleJcn(true);
        double startMinorAP = startBreakend.minorAlleleJcn(true);
        int startMajorAPR = (int)round(startMajorAP);
        int startMinorAPR = (int)round(startMinorAP);
        boolean startMajorAPIsInt = isJcnCloseToInteger(startMajorAP);
        boolean startMinorAPIsInt = isJcnCloseToInteger(startMinorAP);

        // map each major and minor AP into a frequency
        Map<Integer,Integer> jcnFrequency = Maps.newHashMap();

        int segCount = endIndex - startIndex + 1;
        for (int i = startIndex; i <= endIndex; ++i)
        {
            int majorAP = (int)round(alleleJCNs.get(i).MajorAP);
            int minorAP = (int)round(alleleJCNs.get(i).MinorAP);
            Integer repeatCount = jcnFrequency.get(majorAP);
            if(repeatCount == null)
                jcnFrequency.put(majorAP, 1);
            else
                jcnFrequency.put(majorAP, repeatCount+1);

            if(minorAP != majorAP)
            {
                repeatCount = jcnFrequency.get(minorAP);
                if(repeatCount == null)
                    jcnFrequency.put(minorAP, 1);
                else
                    jcnFrequency.put(minorAP, repeatCount+1);
            }
        }

        int maxJcnCount = 0;
        double aJcn = 0;

        for(Map.Entry<Integer,Integer> entry : jcnFrequency.entrySet())
        {
            if(entry.getValue() > maxJcnCount)
            {
                maxJcnCount = entry.getValue();
                aJcn = entry.getKey();
            }
        }

        boolean aJcnValid = maxJcnCount >= segCount * 0.9;

        if(!aJcnValid)
            return false;

        // use knowledge of A to find a minimum for non-disrupted B
        double bJcnMin = -1;

        double startClusterJcn = startBreakend.orientation() == 1 ? startBreakend.jcn() : 0;

        if(startMajorAPIsInt && startMajorAPR == aJcn)
        {
            bJcnMin = startMinorAP - startClusterJcn;
        }
        else if(startMinorAPIsInt && startMinorAPR == aJcn)
        {
            bJcnMin = startMajorAP - startClusterJcn;
        }

        double endClusterJcn = endBreakend.orientation() == -1 ? endBreakend.jcn() : 0;

        for (int i = startIndex; i <= endIndex; ++i)
        {
            int majorAP = (int) round(alleleJCNs.get(i).MajorAP);
            int minorAP = (int) round(alleleJCNs.get(i).MinorAP);

            if(majorAP == aJcn)
            {
                if(i == endIndex)
                    minorAP -= endClusterJcn;

                bJcnMin = bJcnMin == -1 ? minorAP : min(minorAP, bJcnMin);
            }
            else if(minorAP == aJcn)
            {
                if(i == endIndex)
                    majorAP -= endClusterJcn;

                bJcnMin = bJcnMin == -1 ? majorAP : min(majorAP, bJcnMin);
            }
        }

        // set these values int each segment
        for (int i = startIndex; i <= endIndex; ++i)
        {
            alleleJCNs.get(i).AFixedAP = aJcn;
            alleleJCNs.get(i).BUndisruptedAP = bJcnMin;
            alleleJCNs.get(i).setValid(true);
        }

        return true;
    }

    private boolean calculateClusterSegmentJCNs(List<SegmentJcn> alleleJCNs)
    {
        // first establish the non-disrupted B JCN across all segments
        double bNonClusterJcnMin = alleleJCNs.stream()
                .filter(x -> x.isValid())
                .mapToDouble(x -> x.BUndisruptedAP).min().orElse(0);

        // finally set a cluster JCN using knowledge of the other 2 ploidies
        int validSegments = 0;
        int segmentCount = alleleJCNs.size();

        for(int i = 0; i < segmentCount; ++i)
        {
            if(!alleleJCNs.get(i).isValid())
                continue;

            ++validSegments;

            double majorAP = alleleJCNs.get(i).MajorAP;
            double minorAP = alleleJCNs.get(i).MinorAP;
            double aJcn = alleleJCNs.get(i).AFixedAP;

            alleleJCNs.get(i).BUndisruptedAP = bNonClusterJcnMin;

            double clusterJcn = 0;
            if(majorAP > aJcn + 0.5)
            {
                clusterJcn = majorAP - bNonClusterJcnMin;
            }
            else
            {
                clusterJcn = minorAP - bNonClusterJcnMin;
            }

            alleleJCNs.get(i).setClusterJcn(max(clusterJcn, 0.0));
            alleleJCNs.get(i).setValid(true);
        }

        mValidAlleleJcnSegmentPerc = validSegments / (double)segmentCount;

        return true;
    }

    private static final double JCN_INTEGER_PROXIMITY = 0.25;

    private static boolean isJcnCloseToInteger(double jcn)
    {
        double remainder = abs(jcn % 1.0);
        return remainder <= JCN_INTEGER_PROXIMITY || remainder >= (1 - JCN_INTEGER_PROXIMITY);
    }

    public static boolean hasValidAlleleJcnData(int breakendIndex, final List<SegmentJcn> allelePloidies)
    {
        if(allelePloidies == null)
            return false;

        if(allelePloidies.size() < breakendIndex)
            return false;

        return allelePloidies.get(breakendIndex).isValid();
    }

    public boolean linkHasJcnSupport(final LinkedPair pair, double jcn)
    {
        final SvBreakend lowerBreakend = pair.getBreakend(true);
        final SvBreakend upperBreakend = pair.getBreakend(false);

        final List<SegmentJcn> segments = mChrAlleleJCNs.get(pair.chromosome());

        if(segments == null || segments.isEmpty())
            return false;

        int startIndex = lowerBreakend.getClusterChrPosIndex();

        for(int i = startIndex + 1; i < upperBreakend.getClusterChrPosIndex(); ++i)
        {
            final SegmentJcn segment = segments.get(i);

            if(!segment.isValid())
                continue;

            double segmentJcn = segment.unlinkedJcn();

            if(segmentJcn < jcn && !copyNumbersEqual(segmentJcn, jcn))
                return false;
        }

        return true;
    }

    public void assignLinkJcn(final LinkedPair pair, double jcn)
    {
        final SvBreakend lowerBreakend = pair.getBreakend(true);
        final SvBreakend upperBreakend = pair.getBreakend(false);

        final List<SegmentJcn> segments = mChrAlleleJCNs.get(pair.chromosome());

        if(segments == null || segments.isEmpty())
            return;

        for(int i = lowerBreakend.getClusterChrPosIndex(); i < upperBreakend.getClusterChrPosIndex(); ++i)
        {
            SegmentJcn segment = segments.get(i);

            if(!segment.isValid())
                continue;

            segment.addLinkJcn(jcn);
        }
    }

    public static JcnCalcData calcJcnUncertainty(final JcnCalcData data1, final JcnCalcData data2)
    {
        if(!data1.Valid || !data2.Valid)
            return new JcnCalcData(0, 0, false);

        if(data1.JcnUncertainty == 0)
            return new JcnCalcData(data1.JcnEstimate, 0, true);
        else if(data2.JcnUncertainty == 0)
            return new JcnCalcData(data2.JcnEstimate, 0, true);

        double uncertInvSqrd1 = 1 / pow(data1.JcnUncertainty, 2);
        double uncertInvSqrd2 = 1 / pow(data2.JcnUncertainty, 2);

        double sumUncertainty = uncertInvSqrd1 + uncertInvSqrd2;
        double sumObservedUncertainty = data1.JcnEstimate * uncertInvSqrd1 + data2.JcnEstimate * uncertInvSqrd2;

        double estJcn = sumObservedUncertainty / sumUncertainty;

        double adjUncertainty = uncertInvSqrd1 * pow(max(data1.JcnEstimate - estJcn, data1.JcnUncertainty /2),2);
        adjUncertainty += uncertInvSqrd2 * pow(max(data2.JcnEstimate - estJcn, data2.JcnUncertainty /2),2);

        double estUncertainty = sqrt(2 * adjUncertainty / sumUncertainty);

        return new JcnCalcData(estJcn, estUncertainty, true);
    }

    public static final double UNCERTAINTY_SCALE_FACTOR = sqrt(2);

    public static boolean jcnMatchForSplits(double jcnLower, double uncertaintyLower, double jcn, double uncertainty)
    {
        return copyNumbersEqual(jcnLower * 2, jcn)
            || jcnOverlap(jcnLower * 2, uncertaintyLower * UNCERTAINTY_SCALE_FACTOR, jcn, uncertainty);
    }

    public static boolean jcnMatch(final SvVarData var1, final SvVarData var2)
    {
        return jcnMatch(var1.jcn(), var1.jcnUncertainty(), var2.jcn(), var2.jcnUncertainty());
    }

    public static boolean jcnMatch(final SvBreakend breakend1, final SvBreakend breakend2)
    {
        return jcnMatch(breakend1.jcn(), breakend1.jcnUncertainty(), breakend2.jcn(), breakend2.jcnUncertainty());
    }

    public static boolean jcnMatch(double jcn1, double uncertainty1, double jcn2, double uncertainty2)
    {
        return copyNumbersEqual(jcn1, jcn2) || jcnOverlap(jcn1, uncertainty1, jcn2, uncertainty2);
    }

    public static boolean jcnOverlap(double jcn1, double uncertainty1, double jcn2, double uncertainty2)
    {
        if(jcn1 + uncertainty1 < jcn2 - uncertainty2)
            return false;

        if(jcn2 + uncertainty2 < jcn1 - uncertainty1)
            return false;

        return true;
    }

    public static boolean jcnExceedsMajorAlleleJcn(final SvBreakend breakend)
    {
        double adjacentMaJcn = breakend.majorAlleleJcn(breakend.orientation() == -1);
        return breakend.jcn() > adjacentMaJcn && !copyNumbersEqual(breakend.jcn(), adjacentMaJcn);
    }

}
