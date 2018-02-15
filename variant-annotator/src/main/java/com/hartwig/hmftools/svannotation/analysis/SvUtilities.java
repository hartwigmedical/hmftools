package com.hartwig.hmftools.svannotation.analysis;

import com.hartwig.hmftools.common.centromeres.Centromeres;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.svannotation.analysis.SvClusterData;

import java.util.Map;

import static java.lang.Math.abs;


// common utility methods for clustering logic

public class SvUtilities {

    private int mClusterBaseDistance;

    private static final Map<String, GenomeRegion> CENTROMERES = Centromeres.grch37();

    public static String CHROMOSOME_ARM_P = "P"; // short arm, and lower position
    public static String CHROMOSOME_ARM_Q = "Q";
    public static String CHROMOSOME_ARM_CENTROMERE = "C";

    public SvUtilities(int baseDistance)
    {
        mClusterBaseDistance = baseDistance;
    }

    public final String getChromosomalArm(final String chromosome, final long position)
    {
        final GenomeRegion region = CENTROMERES.get(chromosome);

        if(region == null)
            return "INVALID";

        if(position <= region.start())
            return CHROMOSOME_ARM_P;
        if(position >= region.end())
            return CHROMOSOME_ARM_P;
        else
            return CHROMOSOME_ARM_CENTROMERE;
    }

    public boolean areVariantsLinked(final SvClusterData v1, final SvClusterData v2)
    {
        if(v1.id().equals(v2.id()))
            return false;

        if (areVariantsLinked(v1, true, v2, true)
        || areVariantsLinked(v1, false, v2, true)
        || areVariantsLinked(v1, true, v2, false)
        || areVariantsLinked(v1, false, v2, false))
        {
            return true;
        }

        return false;
    }

    public boolean areVariantsLinked(final SvClusterData v1, final boolean v1UseStart, final SvClusterData v2, final boolean v2UseStart)
    {
        // search all remaining SVs for proximity
        if(v1.id().equals(v2.id()))
            return false;

        if(!v1.chromosome(v1UseStart).equals(v2.chromosome(v2UseStart)))
            return false;

        if (!isWithinRange(v1.position(v1UseStart), v2.position(v2UseStart)))
            return false;

        return true;
    }

    public boolean isWithinRange(long pos1, long pos2)
    {
        return abs(pos1 - pos2) <= mClusterBaseDistance;
    }

    public boolean isWithin(final SvClusterData outer, final SvClusterData inner)
    {
        // tests if the inner variant is wholy contained within the outer variant
        if(!outer.isLocal() || !inner.isLocal())
            return false;

        if(!outer.chromosome(true).equals(inner.chromosome(true)))
            return false;

        if(inner.position(true) < outer.position(true))
            return false;

        if(inner.position(false) > outer.position(false))
            return false;

        return true;
    }

    public boolean isWithin(final SvClusterData variant, final String chromosome, final long position)
    {
        if(!variant.chromosome(true).equals(chromosome))
            return false;

        if(variant.position(true) > position || variant.position(false) < position)
            return false;

        return true;
    }

    public boolean isLocalOverlap(final SvClusterData v1, final SvClusterData v2)
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
}
