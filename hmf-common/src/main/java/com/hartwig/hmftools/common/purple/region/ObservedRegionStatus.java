package com.hartwig.hmftools.common.purple.region;

public enum ObservedRegionStatus {
    GERMLINE_HOM_DELETION(true),
    GERMLINE_HET_DELETION(true),
    GERMLINE_AMPLIFICATION(true),
    SOMATIC(false),
    CENTROMERE(false),
    CLUSTER(false),
    UNKNOWN(false);

    private final boolean germlineEvent;


    ObservedRegionStatus(final boolean germlineEvent) {
        this.germlineEvent = germlineEvent;
    }

    public boolean isGermline() {
        return germlineEvent;
    }

}
