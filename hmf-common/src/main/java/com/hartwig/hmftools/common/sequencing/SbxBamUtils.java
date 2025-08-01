package com.hartwig.hmftools.common.sequencing;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

public class SbxBamUtils
{
    public static final byte DUPLEX_QUAL = 93;
    public static final byte SIMPLEX_QUAL = 18;

    public static final String SBX_YC_TAG = "YC";

    public static List<Boolean> getDuplexIndels(final String ycTagStr)
    {
        List<Boolean> duplexIndels = Lists.newArrayList();
        String[] ycTagComponents = ycTagStr.split("-");

        int simplexHeadLength = Integer.parseInt(ycTagComponents[0]);
        String duplexRegion = ycTagComponents[1];
        for(int i = 0; i < simplexHeadLength; i++)
        {
            duplexIndels.add(false);
        }

        for(int i = 0; i < duplexRegion.length();)
        {
            String intString = parseInt(duplexRegion, i);
            if(intString != null)
            {
                int duplexMatchLength = Integer.parseInt(intString);
                for(int j = 0; j < duplexMatchLength; j++)
                {
                    duplexIndels.add(false);
                }
                i += intString.length();
                continue;
            }

            char code = duplexRegion.charAt(i);
            i++;
            switch(code)
            {
                case 'I':
                case 'L':
                case 'P':
                case 'Q':
                case 'J':
                case 'O':
                case 'X':
                case 'Z':
                    duplexIndels.add(true);
                    break;
                default:
                    duplexIndels.add(false);
            }
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
}
