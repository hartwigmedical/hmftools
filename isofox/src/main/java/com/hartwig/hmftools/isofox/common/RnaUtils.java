package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class RnaUtils
{
    public static List<int[]> deriveCommonRegions(final List<int[]> regions1, final List<int[]> regions2)
    {
        // merges any overlapping regoins to create a combined set without overlaps
        if(regions1.isEmpty() || regions2.isEmpty())
            return regions1.isEmpty() ? regions2 : regions1;

        List<int[]> newRegions = Lists.newArrayList();

        // early exit for non-overlapping regions
        if(regions1.get(regions1.size() - 1)[SE_END] < regions2.get(0)[SE_START])
        {
            newRegions.addAll(regions1);
            newRegions.addAll(regions2);
            return newRegions;
        }
        else if(regions2.get(regions2.size() - 1)[SE_END] < regions1.get(0)[SE_START])
        {
            newRegions.addAll(regions2);
            newRegions.addAll(regions1);
            return newRegions;
        }

        int index1 = 0;
        int index2 = 0;

        while(index1 < regions1.size() || index2 < regions2.size())
        {
            int[] region1 = index1 < regions1.size() ? regions1.get(index1) : null;
            int[] region2 = index2 < regions2.size() ? regions2.get(index2) : null;

            if(region1 != null && region2 != null)
            {
                // add the earlier region if not overlapping
                if(region1[SE_END] < region2[SE_START] - 1)
                {
                    newRegions.add(region1);
                    ++index1;
                    continue;
                }
                else if(region2[SE_END] < region1[SE_START] - 1)
                {
                    newRegions.add(region2);
                    ++index2;
                    continue;
                }

                // merge the overlapping regions
                int[] newRegion = new int[] { min(region1[SE_START], region2[SE_START]), max(region1[SE_END], region2[SE_END]) };
                newRegions.add(newRegion);

                ++index1;
                ++index2;

                // continuing merging subsequent overlapping regions
                while(index1 < regions1.size() || index2 < regions2.size())
                {
                    region1 = index1 < regions1.size() ? regions1.get(index1) : null;
                    region2 = index2 < regions2.size() ? regions2.get(index2) : null;

                    boolean merged = false;

                    if(region1 != null && region1[SE_START] <= newRegion[SE_END] + 1)
                    {
                        newRegion[SE_END] = max(region1[SE_END], newRegion[SE_END]);
                        ++index1;
                        merged = true;
                    }

                    if(region2 != null && region2[SE_START] <= newRegion[SE_END] + 1)
                    {
                        newRegion[SE_END] = max(region2[SE_END], newRegion[SE_END]);
                        ++index2;
                        merged = true;
                    }

                    if(!merged)
                        break;
                }
            }
            else if(region1 != null)
            {
                newRegions.add(region1);
                ++index1;
            }
            else
            {
                newRegions.add(region2);
                ++index2;
            }
        }

        return newRegions;
    }

    private static final double MIN_BASE_MATCH_PERC = 0.9;

    public static int findStringOverlaps(final String str1, final String str2)
    {
        if(str1.length() == 0 || str2.length() == 0)
            return 0;

        int matched = 0;
        int i = 0;
        int j = 0;
        int mismatchIndex = -1;

        // first compare bases at same indices, making note of the first difference if there is one
        while(i < str1.length() && j < str2.length())
        {
            if (str1.charAt(i) == str2.charAt(j))
                ++matched;
            else if(mismatchIndex == -1)
                mismatchIndex = i;

            ++i;
            ++j;
        }

        if(matched > MIN_BASE_MATCH_PERC * min(str1.length(), str2.length()))
            return matched;

        i = j = mismatchIndex;
        matched = mismatchIndex;

        while(i < str1.length() && j < str2.length())
        {
            if(str1.charAt(i) == str2.charAt(j))
            {
                ++i;
                ++j;
                ++matched;
                continue;
            }

            // search ahead in each string in turn for the next short matching sequence
            int startI = i;
            boolean seqFound = false;
            for(; i < str1.length() - 2 && j < str2.length() - 2; ++i)
            {
                if(str1.charAt(i) == str2.charAt(j) && str1.charAt(i+1) == str2.charAt(j+1) && str1.charAt(i+2) == str2.charAt(j+2))
                {
                    seqFound = true;
                    break;
                }
            }

            if(seqFound)
                continue;

            i = startI;

            for(; i < str1.length() - 2 && j < str2.length() - 2; ++j)
            {
                if(str1.charAt(i) == str2.charAt(j) && str1.charAt(i+1) == str2.charAt(j+1) && str1.charAt(i+2) == str2.charAt(j+2))
                {
                    seqFound = true;
                    break;
                }
            }

            if(!seqFound)
                break;
        }

        return matched;
    }

    public static StructuralVariantType impliedSvType(final String[] chromosomes, final byte[] orientations)
    {
        if(chromosomes[SE_START].equals(chromosomes[SE_END]))
        {
            if(orientations[SE_START] == orientations[SE_END])
            {
                return INV;
            }
            else if(orientations[SE_START] == POS_ORIENT)
            {
                return DEL;
            }
            else
            {
                return DUP;
            }
        }
        else
        {
            return BND;
        }
    }

    public static final String SP_SEQ_DONOR_1 = "GT";
    public static final String SP_SEQ_DONOR_2 = "GC";
    public static final String SP_SEQ_ACCEPTOR = "AG";

    public static final String SP_SEQ_NEG_STRAND_DONOR_1 = "AC";
    public static final String SP_SEQ_NEG_STRAND_DONOR_2 = "GC";
    public static final String SP_SEQ_NEG_STRAND_ACCEPTOR = "CT";

    public static boolean canonicalAcceptor(final String context, byte strand)
    {
        if(strand == 1)
            return context.equals(SP_SEQ_ACCEPTOR);
        else
            return context.equals(SP_SEQ_NEG_STRAND_ACCEPTOR);
    }

    public static boolean canonicalDonor(final String context, byte strand)
    {
        if(strand == 1)
            return context.equals(SP_SEQ_DONOR_1) || context.equals(SP_SEQ_DONOR_2);
        else
            return context.equals(SP_SEQ_NEG_STRAND_DONOR_1) || context.equals(SP_SEQ_NEG_STRAND_DONOR_2);
    }

    public static long getChromosomeLength(final String chromosome, final RefGenomeVersion version)
    {
        final RefGenomeCoordinates coords = version == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;
        return coords.length(chromosome);
    }

    public static Cigar cigarFromStr(final String cigarStr)
    {
        List<CigarElement> cigarElements = org.apache.commons.compress.utils.Lists.newArrayList();

        int index = 0;
        String basesStr = "";
        while(index < cigarStr.length())
        {
            char c = cigarStr.charAt(index);

            try
            {
                CigarOperator operator = CigarOperator.valueOf(String.valueOf(c));
                cigarElements.add(new CigarElement(Integer.parseInt(basesStr), operator));
                basesStr = "";
            }
            catch (Exception e)
            {
                basesStr += c;
            }
            ++index;
        }

        return new Cigar(cigarElements);
    }
}
