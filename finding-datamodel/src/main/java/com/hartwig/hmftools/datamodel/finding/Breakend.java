package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxGeneOrientation;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record Breakend(
    int id,
    int svId,
    @NotNull String gene,
    @NotNull String chromosome,
    @NotNull String chromosomeBand,
    @NotNull String transcript,
    boolean isCanonical,
    @NotNull LinxGeneOrientation geneOrientation,
    boolean disruptive,
    boolean reported,
    double undisruptedCopyNumber,
    @NotNull LinxBreakendType type,
    @NotNull TranscriptRegionType regionType,
    @NotNull TranscriptCodingType codingType,
    int nextSpliceExonRank,
    int orientation,
    int exonUp,
    int exonDown,
    double junctionCopyNumber)
{
}
