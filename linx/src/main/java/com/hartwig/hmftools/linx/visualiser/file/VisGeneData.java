package com.hartwig.hmftools.linx.visualiser.file;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

public class VisGeneData
{
    public final int ClusterId;
    public final String GeneId;
    public final String GeneName;
    public final String TransName;
    public final int TransId;
    public final String Chromosome;
    public final VisGeneAnnotationType AnnotationType;

    // both sets of exon position adjustments below are used to rectify small differences in pseudogene deletions

    // optional list of exons by rank and any adjustments to their positions for display
    public final Map<Integer,int[]> ExonPositionOffsets;

    // exons by rank and their adjusted positions for exons not retained in this cluster
    public final Map<Integer,int[]> ExonsLostOffsets;

    public final ChrBaseRegion RegionOverride; // used for enhancer regions with matched data in Ensembl

    public VisGeneData(
            int clusterId, final String geneId, final String geneName, final String transName, int transId,
            final String chromosome, final VisGeneAnnotationType annotationType, @Nullable final ChrBaseRegion regionOverride)
    {
        ClusterId = clusterId;
        GeneId = geneId;
        GeneName = geneName;
        TransName = transName;
        TransId = transId;
        Chromosome = chromosome;
        AnnotationType = annotationType;
        RegionOverride = regionOverride;

        ExonPositionOffsets = Maps.newHashMap();
        ExonsLostOffsets = Maps.newHashMap();
    }
}
