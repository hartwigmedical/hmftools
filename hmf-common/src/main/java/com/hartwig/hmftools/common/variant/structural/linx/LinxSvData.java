package com.hartwig.hmftools.common.variant.structural.linx;

import org.immutables.value.Value;

@Value.Immutable
public abstract class LinxSvData
{
    public abstract int svId();
    public abstract int clusterId();
    public abstract String clusterReason();
    public abstract boolean fragileSiteStart();
    public abstract boolean fragileSiteEnd();
    public abstract boolean isFoldback();
    public abstract String lineTypeStart();
    public abstract String lineTypeEnd();
    public abstract double junctionCopyNumberMin();
    public abstract double junctionCopyNumberMax();
    public abstract String geneStart();
    public abstract String geneEnd();
    public abstract double replicationTimingStart();
    public abstract double replicationTimingEnd();
    public abstract int localTopologyIdStart();
    public abstract int localTopologyIdEnd();
    public abstract String localTopologyStart();
    public abstract String localTopologyEnd();
    public abstract int localTICountStart();
    public abstract int localTICountEnd();

}
