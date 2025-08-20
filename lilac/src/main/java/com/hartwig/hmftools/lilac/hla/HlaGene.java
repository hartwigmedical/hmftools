package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.LilacConstants.HLA_PREFIX;
import static com.hartwig.hmftools.lilac.MhcClass_.CLASS_1;
import static com.hartwig.hmftools.lilac.MhcClass_.CLASS_2;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.MhcClass_;

public enum HlaGene
{
    HLA_A(CLASS_1, false),
    HLA_B(CLASS_1, false),
    HLA_C(CLASS_1, false),

    HLA_Y(CLASS_1, true),
    HLA_H(CLASS_1, true),

    HLA_DQB1(CLASS_2, false);

    // TODO: do not include in headers...
//    NONE(CLASS_1, false); // used for debugging

    private final MhcClass_ mMhcClass;
    private final boolean mIsPseudo;

    HlaGene(final MhcClass_ mhcClass, boolean isPseudo)
    {
        mMhcClass = mhcClass;
        mIsPseudo = isPseudo;
    }

    public MhcClass_ mhcClass()
    {
        return mMhcClass;
    }

    public boolean isPseudo()
    {
        return mIsPseudo;
    }

    public static HlaGene fromString(final String s)
    {
        String gene = s.startsWith(HLA_PREFIX) ? s.substring(HLA_PREFIX.length()) : s;
        return valueOf("HLA_" + gene);
    }

    @Override
    public String toString()
    {
        return longName();
    }

    public String longName()
    {
        // TODO:
//        if(this == NONE)
//            return "";

        return super.toString().replace("_", "-");
    }

    public String shortName()
    {
        // TODO:
//        if(this == NONE)
//            return "";

        return super.toString().substring(HLA_PREFIX.length());
    }

    public static List<String> getShortNames(final LilacConfig config)
    {
        List<String> geneStrings = Lists.newArrayList();
        for(HlaGene gene : values())
        {
            if(gene.isPseudo())
                continue;

            if(config.ClassType == null || config.ClassType == gene.mhcClass())
                geneStrings.add(gene.shortName());
        }

        return geneStrings;
    }
}
