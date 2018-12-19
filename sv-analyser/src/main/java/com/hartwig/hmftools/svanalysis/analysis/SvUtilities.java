package com.hartwig.hmftools.svanalysis.analysis;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.centromeres.Centromeres;
import com.hartwig.hmftools.common.chromosome.ChromosomeLengths;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import java.util.List;
import java.util.Map;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;

// common utility methods for clustering logic

public class SvUtilities {

    private static final Map<String, GenomeRegion> CENTROMERES = Centromeres.grch37();
    public static final Map<String,Integer> CHROMOSOME_LENGTHS = ChromosomeLengths.getChromosomeLengths();

    public static String CHROMOSOME_ARM_P = "P"; // short arm, and lower position
    public static String CHROMOSOME_ARM_Q = "Q";
    public static String CHROMOSOME_ARM_CENTROMERE = "C";

    // list of common annotations
    public static String REPLICATION_EVENT = "REP";
    public static String DOUBLE_STRANDED_BREAK = "DSB";
    public static String TEMPLATED_INSERTION = "TI";
    public static String DOUBLE_MINUTE = "DM";
    public static String BREAKAGE_FUSION_BRIDGE = "BFB";
    public static String FRAGILE_SITE = "FRAG";
    public static String UNPHASED_EVENTS = "UNPHAS";

    public static String SV_GROUP_ENCLOSED = "ENCLD";
    public static String SV_GROUP_ENCLOSING = "ENCLG";
    public static String SV_GROUP_OVERLAP = "OVLP";
    public static String SV_GROUP_NEIGHBOURS = "NHRB";

    public static int PERMITED_DUP_BE_DISTANCE = 1;

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

    public static void addSvToChrBreakendMap(final SvVarData var, Map<String, List<SvBreakend>> chrBreakendMap)
    {
        for(int be = SVI_START; be <= SVI_END; ++be)
        {
            if(be == SVI_END && var.isNullBreakend())
                continue;

            boolean useStart = isStart(be);

            final String chr = var.chromosome(useStart);
            long position = var.position(useStart);

            if (!chrBreakendMap.containsKey(chr))
            {
                List<SvBreakend> breakendList = Lists.newArrayList();
                breakendList.add(var.getBreakend(useStart));
                chrBreakendMap.put(chr, breakendList);
                continue;
            }

            // otherwise add the variant in order by ascending position
            List<SvBreakend> breakendList = chrBreakendMap.get(chr);

            int index = 0;
            for (; index < breakendList.size(); ++index)
            {
                final SvBreakend breakend = breakendList.get(index);

                if (position < breakend.position())
                    break;
            }

            breakendList.add(index, var.getBreakend(useStart));
        }
    }

    public static boolean areVariantsLinkedByDistance(final SvVarData v1, final SvVarData v2, int permittedDistance)
    {
        if(v1.id().equals(v2.id()))
            return false;

        for(int be1 = SVI_START; be1 <= SVI_END; ++be1)
        {
            if(be1 == SVI_END && v1.isNullBreakend())
                continue;

            boolean v1Start = isStart(be1);

            for(int be2 = SVI_START; be2 <= SVI_END; ++be2)
            {
                if(be2 == SVI_END && v2.isNullBreakend())
                    continue;

                boolean v2Start = isStart(be2);

                if (areVariantsLinkedByDistance(v1, v1Start, v2, v2Start, permittedDistance))
                    return true;
            }
        }

        return false;
    }

    public static boolean areVariantsLinkedByDistance(final SvVarData v1, final boolean v1UseStart, final SvVarData v2, final boolean v2UseStart, int permittedDistance)
    {
        // search all remaining SVs for proximity
        if(v1 == v2)
            return false;

        if(v1.position(v1UseStart) < 0 || v2.position(v2UseStart) < 0) // for single breakends
            return false;

        if(!v1.chromosome(v1UseStart).equals(v2.chromosome(v2UseStart)))
            return false;

        if (!isWithinRange(v1.position(v1UseStart), v2.position(v2UseStart), permittedDistance))
            return false;

        return true;
    }

    public static boolean isWithinRange(long pos1, long pos2, int permittedDistance)
    {
        if(pos1 < 0 || pos2 < 0)
            return false;

        return abs(pos1 - pos2) <= permittedDistance;
    }

    public static boolean isWithin(final SvVarData outer, final SvVarData inner)
    {
        // tests if the inner variant is wholly contained within the outer variant
        if(!outer.isLocal() || !inner.isLocal())
            return false;

        if(!outer.chromosome(true).equals(inner.chromosome(true)))
            return false;

        if(inner.position(false) < 0 || outer.position(false) < 0)
            return false;

        if(inner.position(true) < outer.position(true))
            return false;

        if(inner.position(false) > outer.position(false))
            return false;

        return true;
    }

    public static boolean isWithin(final SvVarData variant, final String chromosome, final long position)
    {
        if(!variant.chromosome(true).equals(chromosome) || !variant.chromosome(false).equals(chromosome))
            return false;

        if(variant.position(true) > position || variant.position(false) < position)
            return false;

        return true;
    }

    public static boolean isLocalOverlap(final SvVarData v1, final SvVarData v2)
    {
        // tests if the inner variant is wholy contained within the outer variant
        if(!v1.isLocal() || !v2.isLocal())
            return false;

        if(!v1.chromosome(true).equals(v2.chromosome(true)))
            return false;

        // v1's right end overlaps with v2's left
        if(v1.position(true) < v2.position(true)
        && v1.position(false) < v2.position(false)
        && v1.position(false) >= v2.position(true))
        {
            return true;
        }

        // v1's left end overlaps with v2's right
        if(v1.position(false) > v2.position(false)
        && v1.position(true) > v2.position(true)
        && v1.position(true) <= v2.position(false))
        {
            return true;
        }

        return false;
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

    public static int getShortestProximity(final SvVarData v1, final SvVarData v2)
    {
        int minLength = -1;

        if(v1.chromosome(true).equals(v2.chromosome(true)))
        {
            int length = Math.abs((int)v1.position(true) - (int)v2.position(true));
            minLength = (minLength >= 0) ? Math.min(minLength, length) : length;
        }

        if(v1.chromosome(true).equals(v2.chromosome(false)) && !v2.isNullBreakend())
        {
            int length = Math.abs((int)v1.position(true) - (int)v2.position(false));
            minLength = (minLength >= 0) ? Math.min(minLength, length) : length;
        }

        if(v1.chromosome(false).equals( v2.chromosome(true)) && !v1.isNullBreakend())
        {
            int length = Math.abs((int)v1.position(false) - (int)v2.position(true));
            minLength = (minLength >= 0) ? Math.min(minLength, length) : length;
        }

        if(v1.chromosome(false).equals(v2.chromosome(false)) && !v1.isNullBreakend() && !v2.isNullBreakend())
        {
            int length = Math.abs((int)v1.position(false) - (int)v2.position(false));
            minLength = (minLength >= 0) ? Math.min(minLength, length) : length;
        }

        return minLength;
    }

    public static int getProximity(final SvVarData v1, final SvVarData v2, boolean v1Start, boolean v2Start)
    {
        // warning: no check for chromosome or arm
//        if(!v1.chromosome(v1Start).equals(v2.chromosome(v2Start)))
//            return -1;

        if(v1.position(v1Start) < 0 || v2.position(v1Start) < 0)
            return 0;

        return Math.abs((int)v1.position(v1Start) - (int)v2.position(v2Start));
    }

    public static boolean sameChrArm(final SvVarData v1, final SvVarData v2, boolean v1Start, boolean v2Start)
    {
        if(v1.position(v1Start) < 0 || v2.position(v1Start) < 0)
            return false;

        return v1.chromosome(v1Start).equals(v2.chromosome(v2Start)) && v1.arm(v1Start) == v2.arm(v2Start);
    }

    public static boolean variantMatchesBreakend(final SvVarData var, final SvBreakend breakend, boolean useStart, int permittedDist)
    {
        return breakend.chromosome().equals(var.chromosome(useStart))
                && abs(breakend.position() - var.position(useStart)) <= permittedDist
                && breakend.orientation() == var.orientation(useStart);
    }

    public static boolean breakendsMatch(final SvVarData var1, final SvVarData var2, boolean v1Start, boolean v2Start, int permittedDist)
    {
        if(var1.position(v1Start) < 0 || var2.position(v2Start)< 0)
            return false;

        return var1.chromosome(v1Start).equals(var2.chromosome(v2Start))
                && abs(var1.position(v1Start) - var2.position(v2Start)) <= permittedDist
                && var1.orientation(v1Start) == var2.orientation(v2Start);
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
        for(int be = SVI_START; be <= SVI_END; ++be)
        {
            boolean useStart = isStart(be);

            consistencyCount += calcConsistency(var, useStart);
        }

        return consistencyCount;
    }

    public static int calcConsistency(final SvVarData var, boolean useStart)
    {
        if(!useStart && var.isNullBreakend())
            return 0;

        return (var.arm(useStart) == CHROMOSOME_ARM_P ? 1 : -1) * var.orientation(useStart) * 1;
    }

    public static double DEFAULT_MAX_COPY_NUM_DIFF = 0.25;
    public static double DEFAULT_MAX_COPY_NUM_DIFF_PERC = 0.1;

    public static boolean copyNumbersEqual(double cn1, double cn2)
    {
        return copyNumbersEqual(cn1, cn2, DEFAULT_MAX_COPY_NUM_DIFF, DEFAULT_MAX_COPY_NUM_DIFF_PERC);
    }

    public static boolean copyNumbersEqual(double cn1, double cn2, double maxDiff, double maxDiffPerc)
    {
        if(round(cn1) == round(cn1))
            return true;

        double copyNumDiff = abs(cn2 - cn1);
        double copyNumDiffPerc = copyNumDiff / max(abs(cn1), abs(cn2));

        if (copyNumDiff > maxDiff && copyNumDiffPerc > maxDiffPerc)
            return false;

        return true;
    }

    public static int calcTypeCount(final List<SvVarData> svList, StructuralVariantType type)
    {
        int count = 0;

        for(final SvVarData var : svList)
        {
            if(var.type() == type)
                ++count;
        }

        return count;
    }

    public static final String getVariantChrArm(final SvVarData var, boolean isStart)
    {
        return makeChrArmStr(var.chromosome(isStart), var.arm(isStart));
    }

    public static final String makeChrArmStr(final SvVarData var, boolean useStart)
    {
        return makeChrArmStr(var.chromosome(useStart), var.arm(useStart));
    }

    public static final String makeChrArmStr(final String chr, final String arm) { return chr + "_" + arm; }

    public static final String getChrFromChrArm(final String chrArm)
    {
        return chrArm.split("_")[0];
    }

    public static final String getArmFromChrArm(final String chrArm)
    {
        return chrArm.split("_")[1];
    }

    public static int MAX_FACTORIAL_VALUE = 40;

    public static long factorial(long i)
    {
        if(i <= 1 || i > MAX_FACTORIAL_VALUE)
            return 1;

        return i * factorial(i - 1);
    }

    public static long combination(long n, long k)
    {
        long nFact = factorial(n);
        long kFact = factorial(k);
        long nKFact = factorial(n- k);
        long calc = nFact / nKFact;
        calc /= kFact;
        return calc;
    }


}
