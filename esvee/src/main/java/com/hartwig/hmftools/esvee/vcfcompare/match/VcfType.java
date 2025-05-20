package com.hartwig.hmftools.esvee.vcfcompare.match;

import java.io.File;

public enum VcfType
{
    SOMATIC,
    GERMLINE,
    UNFILTERED,
    TRUTH,
    OTHER;

    public static VcfType fromVcfPath(String path)
    {
        String basename = new File(path).getName();

        if(basename.contains("unfiltered"))
        {
            return UNFILTERED;
        }

        if(basename.contains("somatic"))
        {
            return SOMATIC;
        }

        if(basename.contains("germline"))
        {
            return GERMLINE;
        }

        if(basename.contains("truth"))
        {
            return TRUTH;
        }

        return OTHER;
    }
}
