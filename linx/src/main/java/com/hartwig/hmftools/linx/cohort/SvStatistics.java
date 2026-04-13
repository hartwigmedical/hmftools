package com.hartwig.hmftools.linx.cohort;

import com.hartwig.hmftools.common.sv.StructuralVariantType;

public class SvStatistics
{
    public int Total;
    public int NonPON;
    public int[] TypeCounts;
    public int MeiLINE;
    public int MeiSINE;
    public int MeiOther;
    public int Genic;
    public int GeneDisruptive;

    public SvStatistics()
    {
        Total = 0;
        NonPON = 0;
        MeiSINE = 0;
        MeiLINE = 0;
        MeiOther = 0;
        Genic = 0;
        GeneDisruptive = 0;

        TypeCounts = new int[SvCategory.values().length];
    }

    public void addSvData(final SvData sv)
    {
        ++Total;

        if(sv.hasPonMatch())
            ++NonPON;

        if(sv.category() != null)
        {
            ++TypeCounts[sv.category().ordinal()];
        }

        if(sv.IsLine || sv.RepeatClass.contains("LINE"))
            ++MeiLINE;
        else if(sv.RepeatClass.contains("SINE"))
            ++MeiSINE;
        else if(!sv.RepeatClass.isEmpty())
            ++MeiOther;

        if(sv.genic())
        {
            ++Genic;

            if(sv.geneDisruptive() && sv.SvType != StructuralVariantType.SGL)
                ++GeneDisruptive;
        }
    }
}
