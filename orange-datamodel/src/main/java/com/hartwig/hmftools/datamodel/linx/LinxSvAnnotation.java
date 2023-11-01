package com.hartwig.hmftools.datamodel.linx;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LinxSvAnnotation
{
    @NotNull
    String vcfId();

    int svId();

    int clusterId();

    @NotNull
    String clusterReason();

    boolean fragileSiteStart();

    boolean fragileSiteEnd();

    boolean isFoldback();

    @NotNull
    String lineTypeStart();

    @NotNull
    String lineTypeEnd();

    double junctionCopyNumberMin();

    double junctionCopyNumberMax();

    @NotNull
    String geneStart();

    @NotNull
    String geneEnd();

    int localTopologyIdStart();

    int localTopologyIdEnd();

    @NotNull
    String localTopologyStart();

    @NotNull
    String localTopologyEnd();

    int localTICountStart();

    int localTICountEnd();
}
