package com.hartwig.hmftools.linx.analysis;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.centromeres.Centromeres;
import com.hartwig.hmftools.common.chromosome.ChromosomeLengths;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

import java.util.List;
import java.util.Map;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;

// common utility methods for clustering logic

public class SvUtilities {

    private static final Map<String, GenomeRegion> CENTROMERES = Centromeres.grch37();
    public static final Map<String,Integer> CHROMOSOME_LENGTHS = ChromosomeLengths.getChromosomeLengths();

    public final static String CHROMOSOME_ARM_P = "P"; // short arm, and lower position
    public final static String CHROMOSOME_ARM_Q = "Q";
    public final static String CHROMOSOME_ARM_CENTROMERE = "C";

    public final static int NO_LENGTH = -1;

    public static final String getChromosomalArm(final String chromosome, final long position)
    {
        final GenomeRegion region = CENTROMERES.get(chromosome);

        if(region == null)
            return "INVALID";

        if(position <= region.start())
            return CHROMOSOME_ARM_P;
        else if(position >= region.end())
            return CHROMOSOME_ARM_Q;
        else
            return CHROMOSOME_ARM_CENTROMERE;
    }

    public static long getChromosomalArmLength(final String chromosome, final String armType)
    {
        final GenomeRegion region = CENTROMERES.get(chromosome);

        if(region == null)
            return 0;

        if(armType.equals(CHROMOSOME_ARM_P))
        {
            return region.start();
        }

        long chrLength = CHROMOSOME_LENGTHS.get(chromosome);

        return chrLength - region.end();
    }

    public static boolean isShortArmChromosome(final String chromosome)
    {
        return chromosome.equals("13") || chromosome.equals("14") || chromosome.equals("15")
                || chromosome.equals("20") || chromosome.equals("21");
    }

    public static void addSvToChrBreakendMap(final SvVarData var, Map<String, List<SvBreakend>> chrBreakendMap)
    {
        for(int be = SE_START; be <= SE_END; ++be)
        {
            if(be == SE_END && var.isSglBreakend())
                continue;

            SvBreakend breakend = var.getBreakend(be);

            long position = breakend.position();

            List<SvBreakend> breakendList = chrBreakendMap.get(breakend.chromosome());

            if (breakendList == null)
            {
                breakendList = Lists.newArrayList();
                chrBreakendMap.put(breakend.chromosome(), breakendList);
            }

            // add the variant in order by ascending position
            int index = 0;
            for (; index < breakendList.size(); ++index)
            {
                final SvBreakend otherBreakend = breakendList.get(index);

                if (position < otherBreakend.position())
                    break;
            }

            breakendList.add(index, breakend);
        }
    }

    public static int findCentromereBreakendIndex(final List<SvBreakend> breakendList, final String arm)
    {
        if(breakendList == null || breakendList.isEmpty())
            return -1;

        // return the last breakend list index prior to the centromere from either arm direction,
        // returning an index out of bounds if all breakends are on the other arm
        if(arm == CHROMOSOME_ARM_P)
        {
            int i = 0;
            for(; i < breakendList.size(); ++i)
            {
                if(breakendList.get(i).arm() == CHROMOSOME_ARM_Q)
                    break;
            }

            return i - 1;
        }
        else
        {
            int i = breakendList.size() - 1;
            for(; i >= 0; --i)
            {
                if(breakendList.get(i).arm() == CHROMOSOME_ARM_P)
                    break;
            }

            return i + 1 >= breakendList.size() ? -1 : i + 1;
        }
    }

    public static final String getSvTypesStr(final int[] typeCounts)
    {
        // the following map-based naming convention leads
        // to a predictable ordering of types: INV, CRS, BND, DEL and DUP
        String clusterTypeStr = "";

        for(int i = 0;i < typeCounts.length; ++i)
        {
            if(typeCounts[i] == 0)
                continue;

            if(!clusterTypeStr.isEmpty())
                clusterTypeStr += "_";

            clusterTypeStr += StructuralVariantType.values()[i] + "=" + typeCounts[i];
        }

        return clusterTypeStr;
    }

    public static String appendStr(final String dest, final String source, char delim)
    {
        return dest.isEmpty() ? source : dest + delim + source;
    }

    public static boolean isWithinRange(long pos1, long pos2, int permittedDistance)
    {
        if(pos1 < 0 || pos2 < 0)
            return false;

        return abs(pos1 - pos2) <= permittedDistance;
    }

    public static boolean isWithin(final SvVarData variant, final String chromosome, final long position)
    {
        if(!variant.chromosome(true).equals(chromosome) || !variant.chromosome(false).equals(chromosome))
            return false;

        if(variant.position(true) > position || variant.position(false) < position)
            return false;

        return true;
    }

    public static boolean isOverlapping(final SvVarData v1, final SvVarData v2)
    {
        // tests if either variant has an end within the other variant
        if(isWithin(v2, v1.chromosome(true), v1.position(true))
        || isWithin(v2, v1.chromosome(false), v1.position(false)))
        {
            return true;
        }

        if(isWithin(v1, v2.chromosome(true), v2.position(true))
        || isWithin(v1, v2.chromosome(false), v2.position(false)))
        {
            return true;
        }

        return false;
    }

    public static long getProximity(final SvVarData var1, final SvVarData var2)
    {
        long minDistance = -1;

        for(int se1 = SE_START; se1 <= SE_END; ++se1)
        {
            SvBreakend be1 = var1.getBreakend(se1);

            if(be1 == null)
                continue;

            for(int se2 = SE_START; se2 <= SE_END; ++se2)
            {
                SvBreakend be2 = var2.getBreakend(se2);

                if(be2 == null)
                    continue;

                if(be1.chromosome().equals(be2.chromosome()))
                {
                    long distance = abs(be1.position() - be2.position());
                    if(minDistance == -1 || distance < minDistance)
                        minDistance = distance;
                }
            }
        }

        return minDistance;
    }

    public static int calcConsistency(final List<SvVarData> svList)
    {
        int consistencyCount = 0;

        for(final SvVarData var : svList)
        {
            consistencyCount += calcConsistency(var);
        }

        return consistencyCount;
    }

    public static int calcConsistency(final SvVarData var)
    {
        int consistencyCount = 0;
        for(int be = SE_START; be <= SE_END; ++be)
        {
            if(be == SE_END && var.isSglBreakend())
                continue;

            consistencyCount += calcConsistency(var.getBreakend(be));
        }

        return consistencyCount;
    }

    public static int calcConsistency(final SvBreakend breakend)
    {
        return (breakend.arm() == CHROMOSOME_ARM_P ? 1 : -1) * breakend.orientation() * 1;
    }

    public static double MAX_COPY_NUM_DIFF = 0.5;
    public static double MAX_COPY_NUM_DIFF_PERC = 0.15;

    public static boolean copyNumbersEqual(double cn1, double cn2)
    {
        return copyNumbersEqual(cn1, cn2, MAX_COPY_NUM_DIFF, MAX_COPY_NUM_DIFF_PERC);
    }

    public static boolean copyNumbersEqual(double cn1, double cn2, double maxDiff, double maxDiffPerc)
    {
        double copyNumDiff = abs(cn2 - cn1);
        double copyNumDiffPerc = copyNumDiff / max(abs(cn1), abs(cn2));

        if (copyNumDiff > maxDiff && copyNumDiffPerc > maxDiffPerc)
            return false;

        return true;
    }

    public static String formatPloidy(double ploidy)
    {
        if(ploidy > 10)
            return String.format("%.0f", ploidy);
        else if(ploidy < 0.5)
            return String.format("%.2f", ploidy);
        else
            return String.format("%.1f", ploidy);
    }

    public static final String makeChrArmStr(final SvVarData var, boolean useStart)
    {
        return makeChrArmStr(var.chromosome(useStart), var.arm(useStart));
    }

    public static final String makeChrArmStr(final String chr, final String arm) { return chr + "_" + arm; }

}
