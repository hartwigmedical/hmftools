package com.hartwig.hmftools.esvee.utils.vcfcompare;

import java.io.File;

enum SvCaller
{
    ESVEE,
    GRIDSS,
    OTHER;

    public static SvCaller fromVcfPath(String path)
    {
        String basename = new File(path).getName();

        if(basename.contains("esvee"))
        {
            return ESVEE;
        }

        if(basename.contains("gridss") || basename.contains("gripss"))
        {
            return GRIDSS;
        }

        return OTHER;
    }
}
