package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class SvDisruptionData
{
    public final SvVarData Var;
    public final boolean IsStart;

    public final GeneData Gene;
    public final TranscriptData Transcript;
    public final int[] Exons;
    public final TranscriptCodingType CodingType;
    public final TranscriptRegionType RegionType;

    public final double UndisruptedCopyNumber;

    private boolean mReportable;
    private boolean mPseudogeneDeletion;

    public SvDisruptionData(
            final SvVarData var, boolean isStart, final GeneData gene, final TranscriptData transcript, final int[] exons,
            final TranscriptCodingType codingType, final TranscriptRegionType regionType, final double undisruptedCopyNumber)
    {
        Var = var;
        IsStart = isStart;
        Gene = gene;
        Transcript = transcript;
        Exons = exons;
        CodingType = codingType;
        RegionType = regionType;
        UndisruptedCopyNumber = undisruptedCopyNumber;

        mReportable = true;
        mPseudogeneDeletion = false;
    }

    public boolean reportable() { return mReportable; }
    public void setReportable(boolean reportable) { mReportable = reportable; }

    public boolean isPseudogeneDeletion() { return mPseudogeneDeletion; }
    public void markPseudogeneDeletion() { mPseudogeneDeletion = true; }

    public String toString() { return String.format("%s: reportable(%s) gene(%s)",
            Var.toString(), mReportable, Gene.GeneName); }

    public String asCsv()
    {
        StringBuilder sb = new StringBuilder();

        sb.append(String.format("%s,%d,%s,%s,%s,%d,%d",
                mReportable, Var.id(), IsStart, Var.type(),
                Gene.Chromosome, Var.position(IsStart), Var.orientation(IsStart)));

        sb.append(String.format(",%s,%s,%d,%s,%d,%d,%s,%s",
                Gene.GeneId, Gene.GeneName, Gene.Strand, Transcript.TransName,
                Exons[FS_UP], Exons[FS_DOWN], CodingType, RegionType));

        final SvCluster cluster = Var.getCluster();

        sb.append(String.format(",%d,%s,%d",
                cluster.id(), cluster.getResolvedType(), cluster.getSvCount()));

        return sb.toString();
    }
}
