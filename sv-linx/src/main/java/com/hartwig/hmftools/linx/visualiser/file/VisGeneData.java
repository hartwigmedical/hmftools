package com.hartwig.hmftools.linx.visualiser.file;

public class VisGeneData
{
    public final int ClusterId;
    public final String GeneId;
    public final String GeneName;
    public final String TranscriptId;
    public final String Chromosome;
    public final String AnnotationType;


    public VisGeneData(int clusterId, final String geneId, final String geneName, final String transcriptId,
            final String chromosome, final String annotationType)
    {
        ClusterId = clusterId;
        GeneId = geneId;
        GeneName = geneName;
        TranscriptId = transcriptId;
        Chromosome = chromosome;
        AnnotationType = annotationType;
    }
}
