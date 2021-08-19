package com.hartwig.hmftools.linx.fusion;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.linx.types.SvVarData;

public class SvDisruptionData
{
    public final SvVarData Var;
    public final boolean IsStart;

    public final boolean Reportable;

    public final GeneData Gene;
    public final TranscriptData Transcript;
    public final int[] Exons;
    public final TranscriptCodingType CodingType;
    public final TranscriptRegionType RegionType;

    public final double UndisruptedCopyNumber;

    public SvDisruptionData(
            final SvVarData var, boolean isStart, boolean reportable, final GeneData gene, final TranscriptData transcript, final int[] exons,
            final TranscriptCodingType codingType, final TranscriptRegionType regionType, final double undisruptedCopyNumber)
    {
        Var = var;
        IsStart = isStart;
        Reportable = reportable;
        Gene = gene;
        Transcript = transcript;
        Exons = exons;
        CodingType = codingType;
        RegionType = regionType;
        UndisruptedCopyNumber = undisruptedCopyNumber;
    }
}
