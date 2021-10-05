package com.hartwig.hmftools.common.ensemblcache;

public class GeneMappingData
{
    public final String GeneId;
    public final String GeneNameNew;
    public final String GeneNameOld;

    public GeneMappingData(final String geneId, final String geneNameNew, final String geneNameOld)
    {
        GeneId = geneId;
        GeneNameNew = geneNameNew;
        GeneNameOld = geneNameOld;
    }
}
