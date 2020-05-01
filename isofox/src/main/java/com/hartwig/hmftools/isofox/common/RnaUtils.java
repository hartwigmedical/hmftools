package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RnaUtils
{
    public static boolean positionsOverlap(int posStart1, int posEnd1, int posStart2, int posEnd2)
    {
        return !(posStart1 > posEnd2 || posEnd1 < posStart2);
    }

    public static boolean positionWithin(int position, int otherPosStart, int otherPosEnd)
    {
        return (position >= otherPosStart && position <= otherPosEnd);
    }

    public static boolean positionsWithin(int innerStart, int innerEnd, int outerStart, int outerEnd)
    {
        return (innerStart < innerEnd && innerStart >= outerStart && innerEnd <= outerEnd);
    }

    public static List<int[]> deriveCommonRegions(final List<int[]> regions1, final List<int[]> regions2)
    {
        // merges any overlapping regoins to create a combined set without overlaps
        List<int[]> newRegions = Lists.newArrayList();

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

    public static Map<String,Integer> createFieldsIndexMap(final String fieldsHeader, final String delimiter)
    {
        final String[] items = fieldsHeader.split(delimiter,-1);
        final Map<String,Integer> fieldsIndexMap = Maps.newHashMap();

        for(int i = 0; i < items.length; ++i)
        {
            fieldsIndexMap.put(items[i], i);
        }

        return fieldsIndexMap;
    }

    public static void calcPercentileValues(final List<Double> values, final double[] percentileValues)
    {
        int sampleCount = values.size();

        // populate the upper and lower bounds
        double percSlots = percentileValues.length;

        double samplesPerPercentile = sampleCount/percSlots;

        for(int i = 0; i < percentileValues.length; ++i)
        {
            double lowerIndex = i * samplesPerPercentile;
            double upperIndex = lowerIndex + samplesPerPercentile * 0.9999;

            int lowerBound = (int)floor(lowerIndex);
            int upperBound = (int)ceil(upperIndex) - 1;
            upperBound = min(upperBound, sampleCount);

            if(lowerBound == upperBound)
            {
                percentileValues[i] = values.get(lowerBound);
                continue;
            }

            double tpmTotal = 0;
            double sampleTotal = 0;

            for(int s = lowerBound; s <= upperBound; ++s)
            {
                double tpm = values.get(s);

                double fractionOfTpm;

                if(s == lowerBound)
                {
                    fractionOfTpm = 1 - (lowerIndex - lowerBound);
                    sampleTotal += fractionOfTpm;
                    tpmTotal += fractionOfTpm * tpm;
                }
                else if(s == upperBound)
                {
                    fractionOfTpm = upperIndex - upperBound;
                    sampleTotal += fractionOfTpm;
                    tpmTotal += fractionOfTpm * tpm;
                }
                else
                {
                    ++sampleTotal;
                    tpmTotal += tpm;
                }
            }

            percentileValues[i] = tpmTotal / sampleTotal;
        }
    }

    public static StructuralVariantType impliedSvType(final String[] chromosomes, final byte[] orientations)
    {
        if(chromosomes[SE_START].equals(chromosomes[SE_END]))
        {
            if(orientations[SE_START] == orientations[SE_END])
            {
                return INV;
            }
            else if(orientations[SE_START] == 1)
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

    public static void setJunctionBaseContext(
            final IndexedFastaSequenceFile refGenome, final String[] chromosomes, final int[] junctions, final String[] baseContext)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            int startOffset = (se == SE_START) ? 1 : 10;
            int endOffset = startOffset == 1 ? 10: 1;

            baseContext[se] = refGenome.getSubsequenceAt(
                    chromosomes[se], junctions[se] - startOffset, junctions[se] + endOffset).getBaseString();
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

}
