package com.hartwig.hmftools.svtools.neo;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.UNKNOWN;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.TranscriptCodingType;
import com.hartwig.hmftools.common.fusion.TranscriptRegionType;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;

public class NeoEpitopeData
{
    private final PointMutationData mPointMutation;
    private final NeoEpitopeFusion mSvFusion;

    // transcript context
    public final TranscriptData[] TransData; // only start is populated for same-gene NEs

    public final int[] Phases;
    public final int[] ExonRank;

    public TranscriptCodingType[] CodingType;
    public TranscriptRegionType[] RegionType;

    public String UpstreamAcids;
    public String DownstreamAcids;
    public String NovelAcid;
    public int DownstreamNmdBases;

    public NeoEpitopeData(final PointMutationData pointMutation, final NeoEpitopeFusion fusion)
    {
        mPointMutation = pointMutation;
        mSvFusion = fusion;
        TransData = new TranscriptData[] {null, null};
        Phases = new int[] {-1, -1};
        ExonRank = new int[FS_PAIR];
        CodingType = new TranscriptCodingType[] {TranscriptCodingType.UNKNOWN, TranscriptCodingType.UNKNOWN};
        RegionType = new TranscriptRegionType[] {TranscriptRegionType.UNKNOWN, TranscriptRegionType.UNKNOWN};

        UpstreamAcids = "";
        DownstreamAcids = "";
        NovelAcid = "";
        DownstreamNmdBases = 0;
    }

    public final PointMutationData pointMutation() { return mPointMutation; }
    private final NeoEpitopeFusion fusion() { return mSvFusion; }

    public byte orientation(int fs)
    {
        return (fs == FS_UPSTREAM) == (TransData[fs].Strand == POS_STRAND) ? POS_ORIENT : NEG_ORIENT;
    }

    public int positon(int stream)
    {
        if(mSvFusion != null)
            return mSvFusion.Positions[stream];

        int indelBaseDiff = mPointMutation.Alt.length() - mPointMutation.Ref.length();

        int pmPosition = mPointMutation.Position;

        if(indelBaseDiff >= 0)
            return pmPosition;

        if((TransData[FS_UPSTREAM].Strand == POS_STRAND) == (stream == FS_UPSTREAM))
            return pmPosition;
        else
            return pmPosition + abs(indelBaseDiff);
    }

    public String chromosome(int stream)
    {
        if(mSvFusion != null)
            return mSvFusion.Chromosomes[stream];
        else
            return mPointMutation.Chromosome;
    }

    public byte strand(int stream) { return TransData[stream].Strand; }

    public boolean phaseMatched() { return Phases[FS_UPSTREAM] == Phases[FS_DOWNSTREAM]; }

    public String toString()
    {
        if(mSvFusion != null)
        {
            return String.format("fusion up(%s: %s:%d:%d) down(%s: %s:%d:%d)",
                    mSvFusion.GeneNames[FS_UPSTREAM], mSvFusion.Chromosomes[FS_UPSTREAM],
                    mSvFusion.Positions[FS_UPSTREAM], mSvFusion.Orientations[FS_UPSTREAM],
                    mSvFusion.GeneNames[FS_DOWNSTREAM], mSvFusion.Chromosomes[FS_DOWNSTREAM],
                    mSvFusion.Positions[FS_DOWNSTREAM], mSvFusion.Orientations[FS_DOWNSTREAM]);
        }
        else
        {
            return String.format("pointMut (%s: %s:%d %s -> %s)",
                    mPointMutation.Gene, mPointMutation.Chromosome,  mPointMutation.Position,
                    mPointMutation.Ref, mPointMutation.Alt);
        }
    }
}
