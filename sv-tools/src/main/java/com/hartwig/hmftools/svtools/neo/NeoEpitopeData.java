package com.hartwig.hmftools.svtools.neo;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;

import com.hartwig.hmftools.common.fusion.TranscriptCodingType;
import com.hartwig.hmftools.common.fusion.TranscriptRegionType;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;

public class NeoEpitopeData
{
    private final PointMutationData mPointMutation;
    private final NeoEpitopeFusion mSvFusion;

    // transcript context
    public final String[] TransNames; // only start is populated for same-gene NEs

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
        TransNames = new String[FS_PAIR];
        Phases = new int[FS_PAIR];
        ExonRank = new int[FS_PAIR];
        CodingType = new TranscriptCodingType[FS_PAIR];
        RegionType = new TranscriptRegionType[FS_PAIR];

        UpstreamAcids = "";
        DownstreamAcids = "";
        NovelAcid = "";
        DownstreamNmdBases = 0;
    }

    public final PointMutationData pointMutation() { return mPointMutation; }
    private final NeoEpitopeFusion fusion() { return mSvFusion; }

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
