package com.hartwig.hmftools.esvee.vcfcompare.common;

import java.io.File;

enum SvCaller
{
    ESVEE,
    GRIDSS,
    TRUTH,
    OTHER;

    public static SvCaller fromVcfPath(String path)
    {
        String basename = new File(path).getName();

        if(basename.contains("esvee"))
            return ESVEE;

        if(basename.contains("gridss") || basename.contains("gripss"))
            return GRIDSS;

        if(basename.contains("truth"))
            return TRUTH;

        return OTHER;
    }
}
