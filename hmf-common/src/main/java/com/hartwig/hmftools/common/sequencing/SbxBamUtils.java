package com.hartwig.hmftools.common.sequencing;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class SbxBamUtils
{
    public static final byte RAW_DUPLEX_QUAL = 93;
    public static final byte RAW_SIMPLEX_QUAL = 18;
    public static final byte RAW_DUPLEX_MISMATCH_QUAL = 0;

    // values assigned by Redux, used in BQR and all downstream tools
    public static final byte SBX_SIMPLEX_QUAL = 27;
    public static final byte SBX_DUPLEX_QUAL = 40;
    public static final byte SBX_DUPLEX_MISMATCH_QUAL = 1;

    public static final byte SBX_DUPLEX_ADJACENT_1_QUAL_REF_MATCH = 15;
    public static final byte SBX_DUPLEX_ADJACENT_1_QUAL = 20;
    public static final byte SBX_DUPLEX_ADJACENT_2_3_QUAL = 25;

    private static final byte SBX_MEDIUM_QUAL_LOWER = SBX_SIMPLEX_QUAL;
    private static final byte SBX_MEDIUM_QUAL_UPPER = 29;

    public static final String SBX_YC_TAG = "YC";
    public static final String SBX_DUPLEX_READ_INDEX_TAG = "YX";

    public static boolean isHighBaseQual(final byte qual)
    {
        return qual > SBX_MEDIUM_QUAL_UPPER;
    }

    public static boolean isMediumBaseQual(final byte qual)
    {
        return qual >= SBX_MEDIUM_QUAL_LOWER && qual <= SBX_MEDIUM_QUAL_UPPER;
    }

    public static int extractDuplexBaseIndex(final SAMRecord record)
    {
        Integer baseIndex = record.getIntegerAttribute(SBX_DUPLEX_READ_INDEX_TAG);
        return baseIndex != null ? baseIndex : -1;
    }

    public static boolean inDuplexRegion(final SAMRecord record, int baseIndex)
    {
        return inDuplexRegion(!record.getReadNegativeStrandFlag(), extractDuplexBaseIndex(record), baseIndex);
    }

    public static boolean inDuplexRegion(final boolean posOrientationRead, int duplexBaseIndex, int baseIndex)
    {
        if(duplexBaseIndex < 0)
            return false;

        if(posOrientationRead)
            return baseIndex >= duplexBaseIndex;
        else
            return baseIndex <= duplexBaseIndex;
    }

    public static List<Integer> getDuplexIndelIndices(final String ycTagStr)
    {
        List<Integer> duplexIndels = null;
        String[] ycTagComponents = ycTagStr.split("-");

        int simplexHeadLength = Integer.parseInt(ycTagComponents[0]);
        int readIndex = simplexHeadLength;

        String duplexRegion = ycTagComponents[1];

        for(int i = 0; i < duplexRegion.length();)
        {
            String intString = parseInt(duplexRegion, i);
            if(intString != null)
            {
                int duplexMatchLength = Integer.parseInt(intString);

                readIndex += duplexMatchLength;
                i += intString.length();
                continue;
            }

            char code = duplexRegion.charAt(i);
            i++;

            switch(code)
            {
                // indel related (ie not SNV)
                case 'I':
                case 'L':
                case 'P':
                case 'Q':
                case 'J':
                case 'O':
                case 'X':
                case 'Z':

                    if(duplexIndels == null)
                        duplexIndels = Lists.newArrayList();

                    duplexIndels.add(readIndex);
                    break;

                default:
                    break;
            }

            ++readIndex;
        }

        return duplexIndels;
    }

    @Nullable
    private static String parseInt(final String s, int start)
    {
        if(start < 0 || start >= s.length())
            return null;

        if(s.charAt(start) < '0' || s.charAt(start) > '9')
            return null;

        StringBuilder intString = new StringBuilder();
        for(int i = start; i < s.length(); i++)
        {
            if(s.charAt(i) < '0' || s.charAt(i) > '9')
                break;

            intString.append(s.charAt(i));
        }

        return intString.toString();
    }

    public static Boolean isHomopolymerLowBaseQualAtStart(final SAMRecord record)
    {
        byte previousBase = record.getReadBases()[0];
        byte previousQual = record.getBaseQualities()[0];
        byte hpStartQual = 0;
        boolean inHomopolymer = false;

        for(int i = 1; i < record.getReadBases().length; ++i)
        {
            byte base = record.getReadBases()[i];
            byte qual = record.getBaseQualities()[i];

            if(inHomopolymer)
            {
                if(base != previousBase)
                {
                    inHomopolymer = false;

                    if(hpStartQual != previousQual)
                    {
                        if(hpStartQual <= SBX_DUPLEX_MISMATCH_QUAL)
                            return Boolean.TRUE;
                        else if(previousQual <= SBX_DUPLEX_MISMATCH_QUAL)
                            return Boolean.FALSE;
                    }
                }
            }
            else
            {
                if(base == previousBase)
                {
                    inHomopolymer = true;
                    hpStartQual = previousQual;
                }
            }

            previousQual = qual;
            previousBase = base;
        }

        return null;
    }
}
