package com.hartwig.hmftools.datamodel.linx;

import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.sv.LinxBreakendType;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(allParameters = true, passAnnotations = { NotNull.class, Nullable.class })
public interface LinxBreakend {
    int id();

    int svId();

    @NotNull
    String gene();

    @NotNull
    String transcriptId();

    boolean canonical();

    @NotNull
    String geneOrientation();

    boolean disruptive();

    boolean reportedDisruption();

    double undisruptedCopyNumber();

    @NotNull
    TranscriptRegionType regionType();

    @NotNull
    TranscriptCodingType codingType();

    int nextSpliceExonRank();

    // additional fields for patient report
    @NotNull
    LinxBreakendType type();

    @NotNull
    String chromosome();

    int orientation();

    int strand();

    @NotNull
    String chrBand();

    int exonUp();

    int exonDown();
    
    double junctionCopyNumber();
}
