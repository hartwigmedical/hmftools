package com.hartwig.hmftools.linx.visualiser.file;

import java.util.Map;

import com.google.common.collect.Maps;

public class VisGeneData
{
    public final int ClusterId;
    public final String GeneId;
    public final String GeneName;
    public final String TransName;
    public final int TransId;
    public final String Chromosome;
    public final String AnnotationType;

    // optional list of exons by rank and any adjustments to their positions for display
    public final Map<Integer,int[]> ExonPositionOffsets;

    // exons by rank and their adjusted positions for exons not retained in this cluster
    public final Map<Integer,int[]> ExonsLostOffsets;

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

        ExonPositionOffsets = Maps.newHashMap();
        ExonsLostOffsets = Maps.newHashMap();
    }
}
