package com.hartwig.hmftools.esvee.common.saga;

public record SagaLocationMatch(
        SagaVariant variant,
        SagaBreakend breakend,
        int distance
)
{
    public String variantId()
    {
        return variant.id();
    }
}
