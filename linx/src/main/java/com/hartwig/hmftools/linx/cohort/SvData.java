package com.hartwig.hmftools.linx.cohort;

import static com.hartwig.hmftools.common.gene.TranscriptCodingType.NON_CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StartEndPair;
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

    private SvCategory mSvCategory;
    private int mLength;

    private Map<String,PonMatch> mPonMatches;
    private final StartEndPair<List<GeneData>> mGeneDataList;

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

        mSvCategory = null;

        if(svType == DEL ||svType == DUP || svType == StructuralVariantType.INV)
            mLength = posEnd - posStart;
        else if(svType == StructuralVariantType.INS)
            mLength = 1;
        else
            mLength = 0;

        mPonMatches = Maps.newHashMap();

        mGeneDataList = new StartEndPair<>(Lists.newArrayList(), Lists.newArrayList());
        parseGeneData();
    }

    public int length() { return mLength; }

    public void setSvType(final SvCategory type) { mSvCategory = type; }
    public SvCategory category() { return mSvCategory; }

    public boolean genic() { return !GeneStart.isEmpty() || !GeneEnd.isEmpty(); }

    public boolean geneDisruptive()
    {
        return mGeneDataList.start().stream().anyMatch(x -> x.disruptive()) || mGeneDataList.end().stream().anyMatch(x -> x.disruptive());
    }

    private void parseGeneData()
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            String genesDataStr = se == SE_START ? GeneStart : GeneEnd;

            if((se == SE_END && SvType == SGL) || genesDataStr.isEmpty())
                continue;

            List<GeneData> geneDataList = mGeneDataList.get(se);

            String[] geneDataItems = genesDataStr.split(ITEM_DELIM, -1);

            for(String geneDataStr : geneDataItems)
            {
                GeneData geneData = GeneData.fromValues(geneDataStr);
                if(geneData != null)
                    geneDataList.add(geneData);
            }
        }
    }

    public List<GeneData> genesStart() { return mGeneDataList.start(); }
    public List<GeneData> genesEnd() { return mGeneDataList.end(); }

    public void setGeneDisruptive()
    {
        if(mGeneDataList.start().isEmpty() && mGeneDataList.end().isEmpty())
            return;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            for(GeneData geneData : mGeneDataList.get(se))
            {
                if(geneData.CodingType == NON_CODING)
                    geneData.setNonDisruptive();
            }
        }

        // tailored treatment of SGLs
        if(SvType == SGL)
        {
            boolean isUnclusteredSgl = ClusterCount == 1 && SvType == SGL;

            for(GeneData startGeneData : mGeneDataList.start())
            {
                if(!startGeneData.Disruptive)
                    continue;

                if(isUnclusteredSgl && startGeneData.CodingType == UTR_3P)
                    startGeneData.setNonDisruptive();

                // could assume also true for intronic
            }

            return;
        }

        boolean simpleDelDup = ClusterCount == 1 && (SvType == DEL || SvType == DUP);

        // paired breakend logic
        for(GeneData startGeneData : mGeneDataList.start())
        {
            GeneData endGeneData = mGeneDataList.end().stream().filter(x -> x.GeneName.equals(startGeneData.GeneName)).findFirst().orElse(null);

            if(endGeneData == null)
                continue;

            startGeneData.markPaired();
            endGeneData.markPaired();

            if(!startGeneData.Disruptive && !endGeneData.Disruptive)
                continue;

            if(simpleDelDup && (startGeneData.CodingType == UTR_3P && endGeneData.CodingType == UTR_3P))
            {
                startGeneData.setNonDisruptive();
                endGeneData.setNonDisruptive();
                continue;
            }
        }

        // any other exceptions?
    }

    public Map<String,PonMatch> ponMatches() { return mPonMatches; }
    public boolean hasPonMatch() { return !mPonMatches.isEmpty(); }

    protected final String UNNAMED_PON_MATCH = "UNNAMED";
    public void addPonMatch(final PonMatch ponMatch) { mPonMatches.put(UNNAMED_PON_MATCH, ponMatch); }
    public void addPonMatch(final String ponName, final PonMatch ponMatch) { mPonMatches.put(ponName, ponMatch); }

    public PonMatch getUnnamedPonMatch() { return mPonMatches.get(UNNAMED_PON_MATCH); }
}
