package com.hartwig.hmftools.datamodel.linx;

import org.immutables.gson.Gson;
import org.immutables.value.Value;

@Gson.TypeAdapters
@Value.Immutable
public interface LinxSvAnnotation
{
    String vcfId();
    int svId();
    int clusterId();
    String clusterReason();
    boolean fragileSiteStart();
    boolean fragileSiteEnd();
    boolean isFoldback();
    String lineTypeStart();
    String lineTypeEnd();
    double junctionCopyNumberMin();
    double junctionCopyNumberMax();
    String geneStart();
    String geneEnd();
    int localTopologyIdStart();
    int localTopologyIdEnd();
    String localTopologyStart();
    String localTopologyEnd();
    int localTICountStart();
    int localTICountEnd();
}
