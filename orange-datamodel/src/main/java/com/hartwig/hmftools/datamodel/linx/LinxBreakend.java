package com.hartwig.hmftools.datamodel.linx;

import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LinxBreakend
{
    int id();

    int svId();

    @NotNull
    String gene();

    @NotNull
    String chromosome();

    @NotNull
    String chromosomeBand();

    @NotNull
    String transcript();

    boolean isCanonical();

    @NotNull
    String geneOrientation();

    boolean disruptive();

    ReportedStatus reportedStatus();

    double undisruptedCopyNumber();

    @NotNull
    LinxBreakendType type();

    @NotNull
    TranscriptRegionType regionType();

    @NotNull
    TranscriptCodingType codingType();

    int nextSpliceExonRank();

    int orientation();

    int exonUp();

    int exonDown();

    double junctionCopyNumber();
}
