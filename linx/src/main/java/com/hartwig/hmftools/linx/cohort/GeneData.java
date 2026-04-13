package com.hartwig.hmftools.linx.cohort;

import static com.hartwig.hmftools.linx.types.SvVarData.GENE_DATA_ITEM_DELIM;
import static com.hartwig.hmftools.linx.types.SvVarData.SV_DISRUPTIVE_STR;

import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;

public class GeneData
{
    public final String GeneName;
    public final TranscriptCodingType CodingType;
    public final TranscriptRegionType RegionType;
    public final boolean Disruptive;
    public final int Exon;

    private boolean mPaired;
    private boolean mMarkedNonDisruptive;

    public GeneData(
            final String geneName, final TranscriptCodingType codingType, final TranscriptRegionType regionType,
            final boolean disruptive, final int exon)
    {
        GeneName = geneName;
        CodingType = codingType;
        RegionType = regionType;
        Disruptive = disruptive;
        Exon = exon;

        mMarkedNonDisruptive = false;
        mPaired = false;
    }

    public boolean markedNonDisruptive() { return mMarkedNonDisruptive; }
    public void markNonDisruptive() { mMarkedNonDisruptive = true; }
    public boolean disruptive() { return Disruptive && !mMarkedNonDisruptive; }

    public boolean isPaired() { return mPaired; }
    public void markPaired() { mPaired = true; }

    public static GeneData fromValues(final String geneDataStr)
    {
        String[] geneParts = geneDataStr.split("\\" + GENE_DATA_ITEM_DELIM, -1);

        if(geneParts.length < 3)
            return null;

        // format: geneName|regionType|codingType|disruption or not present|exon=123 or not present

        int partIndex = 0;
        String geneName = geneParts[partIndex++];
        TranscriptRegionType regionType = TranscriptRegionType.valueOf(geneParts[partIndex++]);
        TranscriptCodingType codingType = TranscriptCodingType.valueOf(geneParts[partIndex++]);

        boolean isDisruptive = false;

        if(partIndex < geneParts.length && geneParts[partIndex].equals(SV_DISRUPTIVE_STR))
        {
            ++partIndex;
            isDisruptive = true;
        }

        int exon = -1;
        if(partIndex < geneParts.length && geneParts[partIndex].contains("exon"))
        {
            String[] exonData = geneParts[partIndex].split("=", 2);
            exon = Integer.parseInt(exonData[1]);
        }

        return new GeneData(geneName, codingType, regionType, isDisruptive, exon);
    }
}
