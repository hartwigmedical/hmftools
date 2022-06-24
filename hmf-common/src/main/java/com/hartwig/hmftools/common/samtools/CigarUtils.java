package com.hartwig.hmftools.common.samtools;

import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public final class CigarUtils
{
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

    public static int calcCigarLength(final String cigarStr)
    {
        int index = 0;
        int baseLength = 0;
        String basesStr = "";
        while(index < cigarStr.length())
        {
            char c = cigarStr.charAt(index);
            boolean isAddItem = (c == 'D' || c == 'M');
            boolean isIgnoreItem = (c == 'I' || c == 'N' || c == 'S' || c == 'H' || c == 'P' || c == '=' || c == 'X');

            if(isAddItem)
            {
                try { baseLength += Integer.parseInt(basesStr); } catch (Exception e) {}
                basesStr = "";
            }
            else if(isIgnoreItem)
            {
                basesStr = "";
            }
            else
            {
                basesStr += c;
            }

            ++index;
        }

        return baseLength;
    }

}
