package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.LinxConfig.RG_VERSION;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.Q_ARM;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.asStr;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.purple.ChromosomeArm;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

// common utility methods for SVs

public class SvUtilities {

    public static final int NO_LENGTH = -1;

    public static RefGenomeCoordinates refGenomeLengths()
    {
        return RG_VERSION == RefGenomeVersion.V38 ? RefGenomeCoordinates.COORDS_38 : RefGenomeCoordinates.COORDS_37;
    }

    public static int getChromosomeLength(final String chromosome)
    {
        return refGenomeLengths().length(chromosome);
    }

    public static ChromosomeArm getChromosomalArm(final String chromosome, final int position)
    {
        int centromerePos = refGenomeLengths().centromere(chromosome);
        return position < centromerePos ? P_ARM : Q_ARM;
    }

    public static int getChromosomalArmLength(final String chromosome, final ChromosomeArm armType)
    {
        final RefGenomeCoordinates refGenome = refGenomeLengths();
        final HumanChromosome chr = HumanChromosome.fromString(chromosome);

        final Integer centromerePos = refGenome.centromeres().get(chr);

        if(centromerePos == null)
            return 0;

        if(armType == P_ARM)
        {
            return centromerePos.intValue();
        }

        int chrLength = refGenome.lengths().get(chr).intValue();

        return chrLength - centromerePos.intValue();
    }

    public static boolean isShortArmChromosome(final String chromosome)
    {
        return chromosome.equals("13") || chromosome.equals("14") || chromosome.equals("15")
                || chromosome.equals("21") || chromosome.equals("22");
    }

    public static void addSvToChrBreakendMap(final SvVarData var, Map<String, List<SvBreakend>> chrBreakendMap)
    {
        for(int be = SE_START; be <= SE_END; ++be)
        {
            if(be == SE_END && var.isSglBreakend())
                continue;

            SvBreakend breakend = var.getBreakend(be);

            int position = breakend.position();

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

                // special case of inferred placed at another SV's location to explain a CN change
                // ensure the inferred is placed such that it faces its partner breakend
                if(position == otherBreakend.position() && breakend.orientation() != otherBreakend.orientation())
                {
                    if(var.isInferredSgl() && breakend.orientation() == -1)
                        break;
                    else if(otherBreakend.getSV().isInferredSgl() &&  otherBreakend.orientation() == 1)
                        break;
                }
            }

            breakendList.add(index, breakend);
        }
    }

    public static int findCentromereBreakendIndex(final List<SvBreakend> breakendList, final ChromosomeArm arm)
    {
        if(breakendList == null || breakendList.isEmpty())
            return -1;

        // return the last breakend list index prior to the centromere from either arm direction,
        // returning an index out of bounds if all breakends are on the other arm
        if(arm == P_ARM)
        {
            int i = 0;
            for(; i < breakendList.size(); ++i)
            {
                if(breakendList.get(i).arm() == Q_ARM)
                    break;
            }

            return i - 1;
        }
        else
        {
            int i = breakendList.size() - 1;
            for(; i >= 0; --i)
            {
                if(breakendList.get(i).arm() == P_ARM)
                    break;
            }

            return i + 1 >= breakendList.size() ? -1 : i + 1;
        }
    }

    public static String getSvTypesStr(final int[] typeCounts)
    {
        // the following map-based naming convention leads
        // to a predictable ordering of types: INV, CRS, BND, DEL and DUP
        String clusterTypeStr = "";

        for(int i = 0;i < typeCounts.length; ++i)
        {
            if(typeCounts[i] == 0)
                continue;

            clusterTypeStr = appendStr(clusterTypeStr, StructuralVariantType.values()[i] + "=" + typeCounts[i], '_');
        }

        return clusterTypeStr;
    }

    public static boolean isWithin(final SvVarData variant, final String chromosome, final int position)
    {
        if(!variant.chromosome(true).equals(chromosome) || !variant.chromosome(false).equals(chromosome))
            return false;

        return positionWithin(position, variant.position(true), variant.position(false));
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

    public static int getProximity(final SvVarData var1, final SvVarData var2)
    {
        int minDistance = -1;

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
                    int distance = abs(be1.position() - be2.position());
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
        return (breakend.arm() == P_ARM ? 1 : -1) * breakend.orientation();
    }

    public static final double MAX_COPY_NUM_DIFF = 0.5;
    public static final double MAX_COPY_NUM_DIFF_PERC = 0.15;

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
        return String.format("%.3f", ploidy);
    }

    public static String formatJcn(double jcn)
    {
        if(jcn > 10)
            return String.format("%.0f", jcn);
        else if(jcn < 0.5)
            return String.format("%.2f", jcn);
        else
            return String.format("%.1f", jcn);
    }

    public static String makeChrArmStr(final SvVarData var, boolean useStart)
    {
        return makeChrArmStr(var.chromosome(useStart), var.arm(useStart));
    }

    public static String makeChrArmStr(final String chr, final String arm) { return chr + "_" + arm; }
    public static String makeChrArmStr(final String chr, final ChromosomeArm arm) { return makeChrArmStr(chr, asStr(arm)); }

}
