package com.hartwig.hmftools.linx.visualiser.file;

import java.util.List;

import com.google.common.collect.Lists;

public class VisGeneData
{
    public final int ClusterId;
    public final String GeneId;
    public final String GeneName;
    public final String TransName;
    public final int TransId;
    public final String Chromosome;
    public final String AnnotationType;

    public final List<Integer> SpecificExons;

    public VisGeneData(int clusterId, final String geneId, final String geneName, final String transName, int transId,
            final String chromosome, final String annotationType)
    {
        ClusterId = clusterId;
        GeneId = geneId;
        GeneName = geneName;
        TransName = transName;
        TransId = transId;
        Chromosome = chromosome;
        AnnotationType = annotationType;

        SpecificExons = Lists.newArrayList();
    }
}
