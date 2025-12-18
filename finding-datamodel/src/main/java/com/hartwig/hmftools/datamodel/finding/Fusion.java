package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.LinxFusionType;
import com.hartwig.hmftools.datamodel.linx.LinxUnreportableReason;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface Fusion extends Driver
{
    @NotNull
    default String display()
    {
        return String.format("%s::%s", geneStart(), geneEnd());
    }

    @NotNull
    String geneStart();

    @NotNull
    String geneContextStart();

    @NotNull
    String geneTranscriptStart();

    @NotNull
    String geneEnd();

    @NotNull
    String geneContextEnd();

    @NotNull
    String geneTranscriptEnd();

    @NotNull LinxFusionType reportedType();

    @NotNull
    List<LinxUnreportableReason> unreportedReasons();

    @NotNull FusionPhasedType phased();

    int fusedExonUp();

    int fusedExonDown();

    int chainLinks();

    boolean chainTerminated();

    @NotNull
    String domainsKept();

    @NotNull
    String domainsLost();

    double junctionCopyNumber();
}
