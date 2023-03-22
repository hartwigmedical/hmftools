package com.hartwig.hmftools.datamodel.linx;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Gson.TypeAdapters
@Value.Immutable
public interface LinxSvAnnotation {

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
