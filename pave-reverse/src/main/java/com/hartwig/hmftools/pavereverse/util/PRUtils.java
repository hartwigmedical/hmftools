package com.hartwig.hmftools.pavereverse.util;

import com.google.common.base.Preconditions;

public class PRUtils
{
    public static int substitutionDistance(String a, String b)
    {
        Preconditions.checkArgument(a.length() == b.length());
        int distance = 0;
        for(int i = 0; i < a.length(); ++i)
        {
            if(a.charAt(i) != b.charAt(i))
            {
                distance += 1;
            }
        }
        return distance;
    }
}
