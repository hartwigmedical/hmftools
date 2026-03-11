package com.hartwig.hmftools.datamodel.linx;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleAllelicDepth;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LinxFusion
{
    @NotNull
    default String display()
    {
        return String.format("%s::%s", geneUp(), geneDown());
    }

    @NotNull
    String geneUp();

    @NotNull
    String contextUp();

    @NotNull
    String transcriptUp();

    @NotNull
    String geneDown();

    @NotNull
    String contextDown();

    @NotNull
    String transcriptDown();

    @NotNull
    LinxFusionType reportedType();

    @NotNull
    FusionPhasedType phased();

    @NotNull
    DriverInterpretation driverInterpretation();

    int fusedExonUp();

    int fusedExonDown();

    int chainLinks();

    boolean chainTerminated();

    @NotNull
    String domainsKept();

    @NotNull
    String domainsLost();

    double junctionCopyNumber();

    @Nullable
    PurpleAllelicDepth rnaSupport();
}
