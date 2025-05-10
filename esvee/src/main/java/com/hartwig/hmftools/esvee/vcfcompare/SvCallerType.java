package com.hartwig.hmftools.esvee.vcfcompare;

import java.io.File;

public enum SvCallerType
{
    ESVEE,
    GRIDSS,
    TRUTH,
    OTHER;

    public static SvCallerType fromVcfPath(String path)
    {
        String basename = new File(path).getName();

        if(basename.contains("esvee") || basename.contains("purple"))
            return ESVEE;

        if(basename.contains("gridss") || basename.contains("gripss"))
            return GRIDSS;

        if(basename.contains("truth"))
            return TRUTH;

        return OTHER;
    }
}
