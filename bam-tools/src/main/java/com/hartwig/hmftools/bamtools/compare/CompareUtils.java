package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.common.codon.Nucleotides.complement;

import java.io.File;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class CompareUtils
{
    public static SamReaderFactory makeSamReaderFactory(CompareConfig config)
    {
        SamReaderFactory readerFactory = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT);
        if(config.RefGenomeFile != null && !config.RefGenomeFile.isEmpty())
        {
            readerFactory.referenceSequence(new File(config.RefGenomeFile));
        }
        return readerFactory;
    }

    // check if the string match, with the extra caveat that they could be reversed
    public static boolean stringsMatch(String str1, boolean str1Reversed, String str2, boolean str2Reversed)
    {
        if(str1Reversed == str2Reversed)
        {
            return str1.equals(str2);
        }
        // they are not the same orientation, check str1 forward and str2 backwards
        if (str1.length() == str2.length())
        {
            for(int i = 0; i < str1.length(); i++)
            {
                if(str1.charAt(i) != str2.charAt(str2.length() - i - 1))
                {
                    return false;
                }
            }
        }
        return true;
    }

    // check if the bases match, with the extra caveat that they could be reverse complemented
    public static boolean basesMatch(String bases1, boolean bases1Reversed, String bases2, boolean bases2Reversed)
    {
        if(bases1Reversed == bases2Reversed)
        {
            return bases1.equals(bases2);
        }
        // they are not the same orientation, check bases1 forward and bases2 backwards and apply complement
        if (bases1.length() == bases2.length())
        {
            for(int i = 0; i < bases1.length(); i++)
            {
                // need to be complemented
                if(complement(bases1.charAt(i)) != bases2.charAt(bases2.length() - i - 1))
                {
                    return false;
                }
            }
        }
        return true;
    }
}
