package com.hartwig.hmftools.esvee.common.saga;

public record SagaJunctionMatchInfo(
        SagaJunctionInfo junctionInfo,
        int alignOverlapLeft,
        int alignOverlapRight,
        boolean indelNearby
)
{
}
