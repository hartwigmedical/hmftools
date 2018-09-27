package com.hartwig.hmftools.svanalysis.analysis;

import com.hartwig.hmftools.common.centromeres.Centromeres;
import com.hartwig.hmftools.common.chromosome.ChromosomeLengths;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;

import java.util.List;
import java.util.Map;

import static java.lang.Math.abs;


// common utility methods for clustering logic

public class SvUtilities {

    private int mClusterBaseDistance;

    private static final Map<String, GenomeRegion> CENTROMERES = Centromeres.grch37();
    final Map<String,Integer> CHROMOSOME_LENGTHS = ChromosomeLengths.getChromosomeLengths();

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

    public SvUtilities(int baseDistance)
    {
        mClusterBaseDistance = baseDistance;
    }

    public int getBaseDistance() { return mClusterBaseDistance; }

    public final String getChromosomalArm(final String chromosome, final long position)
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

    public long getChromosomalArmLength(final String chromosome, final String armType)
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

    public boolean areVariantsLinkedByDistance(final SvClusterData v1, final SvClusterData v2)
    {
        if(v1.id().equals(v2.id()))
            return false;

        if (this.areVariantsLinkedByDistance(v1, true, v2, true)
        || this.areVariantsLinkedByDistance(v1, false, v2, true)
        || this.areVariantsLinkedByDistance(v1, true, v2, false)
        || this.areVariantsLinkedByDistance(v1, false, v2, false))
        {
            return true;
        }

        return false;
    }

    public boolean areVariantsLinkedByDistance(final SvClusterData v1, final boolean v1UseStart, final SvClusterData v2, final boolean v2UseStart)
    {
        // search all remaining SVs for proximity
        if(v1.id().equals(v2.id()))
            return false;

        if(v1.position(v1UseStart) < 0 || v2.position(v2UseStart) < 0) // for single breakends
            return false;

        if(!v1.chromosome(v1UseStart).equals(v2.chromosome(v2UseStart)))
            return false;

        if (!isWithinRange(v1.position(v1UseStart), v2.position(v2UseStart)))
            return false;

        return true;
    }

    public boolean areAllWithinRange(long start1, long end1, long start2, long end2, int permittedDistance)
    {
        if(start1 > start2 - permittedDistance && end1 < end2 + permittedDistance)
            return true;

        if(start2 > start1 - permittedDistance && end2 < end1 + permittedDistance)
            return true;

        return false;
    }

    public boolean areAnyWithinRange(long start1, long end1, long start2, long end2, int permittedDistance)
    {
        return isWithinRange(start1, start2) || isWithinRange(end1, end2)
            || isWithinRange(end1, start2) || isWithinRange(start1, end2);
    }

    public static boolean areTypePair(final SvClusterData v1, final SvClusterData v2, StructuralVariantType type1, StructuralVariantType type2)
    {
        return (v1.type() == type1 && v2.type() == type2) || (v1.type() == type2 && v2.type() == type1);
    }

    public boolean isWithinRange(long pos1, long pos2)
    {
        if(pos1 < 0 || pos2 < 0)
            return false;

        return abs(pos1 - pos2) <= mClusterBaseDistance;
    }

    public static boolean isWithin(final SvClusterData outer, final SvClusterData inner)
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

    public static boolean isWithin(final SvClusterData variant, final String chromosome, final long position)
    {
        if(!variant.chromosome(true).equals(chromosome) || !variant.chromosome(false).equals(chromosome))
            return false;

        if(variant.position(true) > position || variant.position(false) < position)
            return false;

        return true;
    }

    public static boolean isLocalOverlap(final SvClusterData v1, final SvClusterData v2)
    {
        // tests if the inner variant is wholy contained within the outer variant
        if(!v1.isLocal() || !v2.isLocal())
            return false;

        if(!v1.chromosome(true).equals(v2.chromosome(true)))
            return false;

        // v1's right end overlaps with v2's left
        if(v1.position(true) < v2.position(true)
        && v1.position(false) < v2.position(false)
        && v1.position(false) >= v2.position(true)) {
            return true;
        }

        // v1's left end overlaps with v2's right
        if(v1.position(false) > v2.position(false)
        && v1.position(true) > v2.position(true)
        && v1.position(true) <= v2.position(false)) {
            return true;
        }

        return false;
    }

    public static boolean isOverlapping(final SvClusterData v1, final SvClusterData v2)
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

    public static int getShortestProximity(final SvClusterData v1, final SvClusterData v2)
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

    public static int getProximity(final SvClusterData v1, final SvClusterData v2, boolean v1Start, boolean v2Start)
    {
        // warning: no check for chromosome or arm
//        if(!v1.chromosome(v1Start).equals(v2.chromosome(v2Start)))
//            return -1;

        if(v1.position(v1Start) < 0 || v2.position(v1Start) < 0)
            return 0;

        return Math.abs((int)v1.position(v1Start) - (int)v2.position(v2Start));
    }

    public static boolean isFromCentromere(final SvClusterData var)
    {
        return (var.orientation(true) == 1) == (var.getStartArm() == CHROMOSOME_ARM_Q);
    }

    public static boolean isFromTelomere(final SvClusterData var)
    {
        return (var.orientation(true) == 1) == (var.getStartArm() == CHROMOSOME_ARM_P);
    }

    public static boolean isToCentromere(final SvClusterData var)
    {
        return (var.orientation(false) == -1) == (var.getStartArm() == CHROMOSOME_ARM_P);
    }

    public static boolean isToTelomere(final SvClusterData var)
    {
        return (var.orientation(false) == -1) == (var.getStartArm() == CHROMOSOME_ARM_Q);
    }

    public static boolean sameChrArm(final SvClusterData v1, final SvClusterData v2, boolean v1Start, boolean v2Start)
    {
        if(v1.position(v1Start) < 0 || v2.position(v1Start) < 0)
            return false;

        return v1.chromosome(v1Start).equals(v2.chromosome(v2Start)) && v1.arm(v1Start) == v2.arm(v2Start);
    }

    public static boolean areLinkedSection(final SvClusterData v1, final SvClusterData v2, boolean v1Start, boolean v2Start)
    {
        // must be same chromosomal arm to be considered facing
        if(!sameChrArm(v1, v2, v1Start, v2Start))
            return false;

        // start apart and heading towards each other
        long pos1 = v1.position(v1Start);
        boolean headsLeft1 = (v1.orientation(v1Start) == 1);
        long pos2 = v2.position(v2Start);
        boolean headsLeft2 = (v2.orientation(v2Start) == 1);

        if(pos1 < pos2 && !headsLeft1 && headsLeft2)
            return true;

        if(pos2 < pos1 && headsLeft1 && !headsLeft2)
            return true;

        return false;
    }

    public static boolean areSectionBreak(final SvClusterData v1, final SvClusterData v2, boolean v1Start, boolean v2Start)
    {
        // only relevant if on same chromosomal arm
        if(!sameChrArm(v1, v2, v1Start, v2Start))
            return false;

        // start apart or equal and heading same direction
        long pos1 = v1.position(v1Start);
        boolean headsLeft1 = (v1.orientation(v1Start) == 1);
        long pos2 = v2.position(v2Start);
        boolean headsLeft2 = (v2.orientation(v2Start) == 1);

        if(pos1 <= pos2 && headsLeft1 && !headsLeft2)
            return true;

        if(pos2 <= pos1 && !headsLeft1 && headsLeft2)
            return true;

        return false;
    }

    public static boolean variantMatchesBreakend(final SvClusterData var, final SvBreakend breakend, boolean useStart, int permittedDist)
    {
        return breakend.chromosome().equals(var.chromosome(useStart))
                && abs(breakend.position() - var.position(useStart)) <= permittedDist
                && breakend.orientation() == var.orientation(useStart);
    }

    public static boolean breakendsMatch(final SvClusterData var1, final SvClusterData var2, boolean v1Start, boolean v2Start, int permittedDist)
    {
        if(var1.position(v1Start) < 0 || var2.position(v2Start)< 0)
            return false;

        return var1.chromosome(v1Start).equals(var2.chromosome(v2Start))
                && abs(var1.position(v1Start) - var2.position(v2Start)) <= permittedDist
                && var1.orientation(v1Start) == var2.orientation(v2Start);
    }

    public static int calcConsistency(final List<SvClusterData> svList)
    {
        int consistencyCount = 0;

        for(final SvClusterData var : svList)
        {
            consistencyCount += (var.getStartArm() == CHROMOSOME_ARM_P ? 1 : -1) * var.orientation(true);
            consistencyCount += (var.getEndArm() == CHROMOSOME_ARM_P ? 1 : -1) * var.orientation(false);
        }

        return consistencyCount;
    }

    public static final String getVariantChrArm(final SvClusterData var, boolean isStart)
    {
        return var.chromosome(isStart) + "_" + var.arm(isStart);
    }

    public static final String getChrFromChrArm(final String chrArm)
    {
        return chrArm.split("_")[0];
    }

    public static final String getArmFromChrArm(final String chrArm)
    {
        return chrArm.split("_")[1];
    }

}
