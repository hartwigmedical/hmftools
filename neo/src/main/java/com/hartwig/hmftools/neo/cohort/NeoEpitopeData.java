package com.hartwig.hmftools.neo.cohort;

import com.hartwig.hmftools.common.neo.NeoEpitopeType;

public class NeoEpitopeData
{
    public final int Id;
    public final NeoEpitopeType VariantType;
    public final String VariantInfo;
    public final String GeneName;

    public NeoEpitopeData(final int id, final NeoEpitopeType variantType, final String variantInfo, final String geneName)
    {
        Id = id;
        VariantType = variantType;
        VariantInfo = variantInfo;
        GeneName = geneName;
    }
}
