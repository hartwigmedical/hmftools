package com.hartwig.hmftools.linx.cohort;

import static com.hartwig.hmftools.linx.types.SvVarData.SV_DISRUPTIVE_STR;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.linx.types.ResolvedType;

public class SvData
{
    public final int Id;
    public final String ChrStart;
    public final int PosStart;
    public final Orientation OrientStart;
    public final String ChrEnd;
    public final int PosEnd;
    public final Orientation OrientEnd;

    public final StructuralVariantType SvType;
    public final ResolvedType Resolved;

    public final String GeneStart;
    public final String GeneEnd;

    public final int ClusterId;
    public final int ClusterCount;
    public final boolean IsLine;

    public final int ChainId;
    public final int ChainCount;
    public final String ChainIndex;

    public final double AF;
    public final String RepeatClass;

    private SvType mType;
    private int mLength;

    private boolean mGeneDisruptive;
    private PonMatchType mPonMatch;

    public SvData(
            final int id, final String chrStart, final int posStart, final Orientation orientStart, final String chrEnd,
            final int posEnd, final Orientation orientEnd, final StructuralVariantType svType, final ResolvedType resolved,
            final String geneStart, final String geneEnd, final int clusterId, final int clusterCount, final boolean isLine,
            final int chainId, final int chainCount, final String chainIndex, final double af, final String repeatClass)
    {
        Id = id;
        ChrStart = chrStart;
        PosStart = posStart;
        OrientStart = orientStart;
        ChrEnd = chrEnd;
        PosEnd = posEnd;
        OrientEnd = orientEnd;
        SvType = svType;
        Resolved = resolved;
        GeneStart = geneStart;
        GeneEnd = geneEnd;
        ClusterId = clusterId;
        ClusterCount = clusterCount;
        IsLine = isLine;
        ChainId = chainId;
        ChainCount = chainCount;
        ChainIndex = chainIndex;
        AF = af;
        RepeatClass = repeatClass;

        mType = null;

        if(svType == StructuralVariantType.DEL ||svType == StructuralVariantType.DUP || svType == StructuralVariantType.INV)
            mLength = posEnd - posStart;
        else if(svType == StructuralVariantType.INS)
            mLength = 1;
        else
            mLength = 0;

        mGeneDisruptive = false;
        mPonMatch = PonMatchType.NONE;
    }

    public int length() { return mLength; }

    public void setSvType(final SvType type) { mType = type; }
    public SvType svType() { return mType; }

    public boolean genic() { return !GeneStart.isEmpty() || !GeneEnd.isEmpty(); }
    public boolean geneDisruptive() { return mGeneDisruptive; }

    public void setGeneDisruptive()
    {
        if(!genic())
            return;

        if(GeneStart.contains(SV_DISRUPTIVE_STR) || GeneEnd.contains(SV_DISRUPTIVE_STR))
        {
            mGeneDisruptive = true;
        }

        // any other exceptions?
    }

    public PonMatchType ponMatch() { return mPonMatch; }
    public void setPonMatch(final PonMatchType ponMatch) { mPonMatch = ponMatch; }

}
