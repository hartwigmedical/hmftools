package com.hartwig.hmftools.common.peach;

import static com.hartwig.hmftools.common.peach.PeachUtil.convertToHaplotypeString;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PeachGenotype
{
    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String allele();

    public abstract int alleleCount();

    @NotNull
    public abstract String function();

    @NotNull
    public abstract String linkedDrugs();

    @NotNull
    public abstract String urlPrescriptionInfo();

    @Deprecated
    @Value.Default
    public String panelVersion()
    {
        return "";
    };

    @Deprecated
    @Value.Default
    public String repoVersion()
    {
        return "";
    };

    @NotNull
    public String haplotype()
    {
        return convertToHaplotypeString(allele(), alleleCount());
    }
}