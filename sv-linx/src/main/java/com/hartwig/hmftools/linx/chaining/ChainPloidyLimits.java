package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.cn.PloidyCalcData;
import com.hartwig.hmftools.linx.types.SvBreakend;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ChainPloidyLimits
{
    private Map<String, double[][]> mChrAllelePloidies;
    private double mValidAllelePloidySegmentPerc;

    // referenced state
    private Map<String, List<SvBreakend>> mChrBreakendMap;
    private int mClusterId;

    // ploidy level below which a chain segment cannot cross
    public static final double CLUSTER_ALLELE_PLOIDY_MIN = 0.15;

    public static final int MAJOR_AP = 0;
    public static final int MINOR_AP = 1;
    public static final int A_FIXED_AP = 2;
    public static final int B_NON_DIS_AP = 3;
    public static final int CLUSTER_AP = 4;
    public static final int AP_DATA_VALID = 5;
    public static final int AP_IS_VALID = 1;

    private static final Logger LOGGER = LogManager.getLogger(ChainPloidyLimits.class);

    public ChainPloidyLimits()
    {
        mChrAllelePloidies = Maps.newHashMap();
        mClusterId = -1;
        mChrBreakendMap = null;
    }

    public final Map<String, double[][]> getChrAllelePloidies() { return mChrAllelePloidies; }
    public double getValidAllelePloidySegmentPerc() { return mValidAllelePloidySegmentPerc; }

    public void initialise(int clusterId, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        mChrBreakendMap = chrBreakendMap;
        mClusterId = clusterId;

        mChrAllelePloidies.clear();
        mValidAllelePloidySegmentPerc = 0;
    }

    public void determineBreakendPloidies()
    {
        int totalSegCount = 0;
        int totalValidSegCount = 0;

        for (final Map.Entry<String, List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            final String chromosome = entry.getKey();
            final List<SvBreakend> breakendList = entry.getValue();
            int breakendCount = breakendList.size();

            // a multi-dim array of breakend index for this arm to A allele ploidy, B non-disrupted ploidy, and cluster ploidy
            double[][] allelePloidies = new double[breakendCount][AP_DATA_VALID +1];

            mChrAllelePloidies.put(chromosome, allelePloidies);

            boolean inSegment = false;
            SvBreakend segStartBreakend = null;
            int segStartIndex = 0;

            for (int i = 0; i < breakendList.size(); ++i)
            {
                SvBreakend breakend = breakendList.get(i);

                allelePloidies[i][MAJOR_AP] = breakend.majorAllelePloidy(false);
                allelePloidies[i][MINOR_AP] = breakend.minorAllelePloidy(false);

                if(!inSegment)
                {
                    inSegment = true;
                    segStartBreakend = breakend;
                    segStartIndex = i;
                }
                else
                {
                    SvBreakend nextBreakend = (i < breakendList.size() - 1) ? breakendList.get(i + 1) : null;
                    if(nextBreakend == null || nextBreakend.getChrPosIndex() > breakend.getChrPosIndex() + 1)
                    {
                        // a gap in the cluster so need to evaluate this contiguous section
                        inSegment = false;
                        boolean validCalcs = calculateNoneClusterSegmentPloidies(
                                allelePloidies, segStartIndex, i, segStartBreakend, breakend);
                        ++totalSegCount;

                        if(validCalcs)
                            ++totalValidSegCount;
                    }
                }
            }

            calculateClusterSegmentPloidies(allelePloidies);
        }

        LOGGER.debug("cluster({}) chromosomes({}) AP totalSegments({} valid={})",
                mClusterId, mChrBreakendMap.size(), totalSegCount, totalValidSegCount);
    }

    private boolean calculateNoneClusterSegmentPloidies(
            double[][] allelePloidies, int startIndex, int endIndex, SvBreakend startBreakend, SvBreakend endBreakend)
    {
        double startMajorAP = startBreakend.majorAllelePloidy(true);
        double startMinorAP = startBreakend.minorAllelePloidy(true);
        int startMajorAPR = (int)round(startMajorAP);
        int startMinorAPR = (int)round(startMinorAP);
        boolean startMajorAPIsInt = isPloidyCloseToInteger(startMajorAP);
        boolean startMinorAPIsInt = isPloidyCloseToInteger(startMinorAP);

        // map each major and minor AP into a frequency
        Map<Integer,Integer> ploidyFrequency = new HashMap();

        int segCount = endIndex - startIndex + 1;
        for (int i = startIndex; i <= endIndex; ++i)
        {
            int majorAP = (int)round(allelePloidies[i][MAJOR_AP]);
            int minorAP = (int)round(allelePloidies[i][MINOR_AP]);
            Integer repeatCount = ploidyFrequency.get(majorAP);
            if(repeatCount == null)
                ploidyFrequency.put(majorAP, 1);
            else
                ploidyFrequency.put(majorAP, repeatCount+1);

            if(minorAP != majorAP)
            {
                repeatCount = ploidyFrequency.get(minorAP);
                if(repeatCount == null)
                    ploidyFrequency.put(minorAP, 1);
                else
                    ploidyFrequency.put(minorAP, repeatCount+1);
            }
        }

        int maxPloidyCount = 0;
        double aPloidy = 0;

        for(Map.Entry<Integer,Integer> entry : ploidyFrequency.entrySet())
        {
            if(entry.getValue() > maxPloidyCount)
            {
                maxPloidyCount = entry.getValue();
                aPloidy = entry.getKey();
            }
        }

        boolean aPloidyValid = maxPloidyCount >= segCount * 0.9;

        if(!aPloidyValid)
            return false;

        // use knowledge of A to find a minimum for non-disrupted B
        double bPloidyMin = -1;

        double startClusterPloidy = startBreakend.orientation() == 1 ? startBreakend.ploidy() : 0;

        if(startMajorAPIsInt && startMajorAPR == aPloidy)
        {
            bPloidyMin = startMinorAP - startClusterPloidy;
        }
        else if(startMinorAPIsInt && startMinorAPR == aPloidy)
        {
            bPloidyMin = startMajorAP - startClusterPloidy;
        }

        double endClusterPloidy = endBreakend.orientation() == -1 ? endBreakend.ploidy() : 0;

        for (int i = startIndex; i <= endIndex; ++i)
        {
            int majorAP = (int) round(allelePloidies[i][MAJOR_AP]);
            int minorAP = (int) round(allelePloidies[i][MINOR_AP]);

            if(majorAP == aPloidy)
            {
                if(i == endIndex)
                    minorAP -= endClusterPloidy;

                bPloidyMin = bPloidyMin == -1 ? minorAP : min(minorAP, bPloidyMin);
            }
            else if(minorAP == aPloidy)
            {
                if(i == endIndex)
                    majorAP -= endClusterPloidy;

                bPloidyMin = bPloidyMin == -1 ? majorAP : min(majorAP, bPloidyMin);
            }
        }

        // set these values int each segment
        for (int i = startIndex; i <= endIndex; ++i)
        {
            allelePloidies[i][A_FIXED_AP] = aPloidy;
            allelePloidies[i][B_NON_DIS_AP] = bPloidyMin;
            allelePloidies[i][AP_DATA_VALID] = AP_IS_VALID;
        }

        return true;
    }

    private boolean calculateClusterSegmentPloidies(double[][] allelePloidies)
    {
        // first establish the non-disrupted B ploidy across all segments
        double bNonClusterPloidyMin = allelePloidies[0][B_NON_DIS_AP];
        for(int i = 1; i < allelePloidies.length; ++i)
        {
            bNonClusterPloidyMin = min(bNonClusterPloidyMin, allelePloidies[i][B_NON_DIS_AP]);
        }

        // finally set a cluster ploidy using knowledge of the other 2 ploidies
        int validSegments = 0;
        int segmentCount = allelePloidies.length;
        for(int i = 0; i < segmentCount; ++i)
        {
            if(allelePloidies[i][AP_DATA_VALID] != AP_IS_VALID)
                continue;

            ++validSegments;

            double majorAP = allelePloidies[i][MAJOR_AP];
            double minorAP = allelePloidies[i][MINOR_AP];
            double aPloidy = allelePloidies[i][A_FIXED_AP];

            allelePloidies[i][B_NON_DIS_AP] = bNonClusterPloidyMin;

            double clusterPloidy = 0;
            if(majorAP > aPloidy + 0.5)
            {
                clusterPloidy = majorAP - bNonClusterPloidyMin;
            }
            else
            {
                clusterPloidy = minorAP - bNonClusterPloidyMin;
            }

            allelePloidies[i][CLUSTER_AP] = max(clusterPloidy, 0.0);
            allelePloidies[i][AP_DATA_VALID] = AP_IS_VALID;
        }

        mValidAllelePloidySegmentPerc = validSegments / (double)segmentCount;

        return true;
    }

    private static double PLOIDY_INTEGER_PROXIMITY = 0.25;

    private static boolean isPloidyCloseToInteger(double ploidy)
    {
        double remainder = abs(ploidy % 1.0);
        return remainder <= PLOIDY_INTEGER_PROXIMITY || remainder >= (1 - PLOIDY_INTEGER_PROXIMITY);
    }

    public static boolean hasValidAllelePloidyData(int breakendIndex, final double[][] allelePloidies)
    {
        if(allelePloidies == null)
            return false;

        if(allelePloidies.length < breakendIndex)
            return false;

        final double[] beAllelePloidies = allelePloidies[breakendIndex];

        return beAllelePloidies[AP_DATA_VALID] == AP_IS_VALID;
    }

    public static PloidyCalcData calcPloidyUncertainty(final PloidyCalcData data1, final PloidyCalcData data2)
    {
        if(!data1.Valid || !data2.Valid)
            return new PloidyCalcData(0, 0, false);

        if(data1.PloidyUncertainty == 0)
            return new PloidyCalcData(data1.PloidyEstimate, 0, true);
        else if(data2.PloidyUncertainty == 0)
            return new PloidyCalcData(data2.PloidyEstimate, 0, true);

        double uncertInvSqrd1 = 1 / pow(data1.PloidyUncertainty, 2);
        double uncertInvSqrd2 = 1 / pow(data2.PloidyUncertainty, 2);

        double sumUncertainty = uncertInvSqrd1 + uncertInvSqrd2;
        double sumObservedUncertainty = data1.PloidyEstimate * uncertInvSqrd1 + data2.PloidyEstimate * uncertInvSqrd2;

        double estPloidy = sumObservedUncertainty / sumUncertainty;

        double adjUncertainty = uncertInvSqrd1 * pow(max(data1.PloidyEstimate - estPloidy, data1.PloidyUncertainty/2),2);
        adjUncertainty += uncertInvSqrd2 * pow(max(data2.PloidyEstimate - estPloidy, data2.PloidyUncertainty/2),2);


        double estUncertainty = sqrt(2 * adjUncertainty / sumUncertainty);

        return new PloidyCalcData(estPloidy, estUncertainty, true);
    }

    public static boolean ploidyMatch(double ploidy1, double uncertainty1, double ploidy2, double uncertainty2)
    {
        return copyNumbersEqual(ploidy1, ploidy2) || ploidyOverlap(ploidy1, uncertainty1, ploidy2, uncertainty2);
    }

    public static boolean ploidyOverlap(double ploidy1, double uncertainty1, double ploidy2, double uncertainty2)
    {
        if(ploidy1 + uncertainty1 < ploidy2 - uncertainty2)
            return false;

        if(ploidy2 + uncertainty2 < ploidy1 - uncertainty1)
            return false;

        return true;
    }


}
