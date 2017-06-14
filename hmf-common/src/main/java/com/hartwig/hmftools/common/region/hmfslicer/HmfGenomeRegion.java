package com.hartwig.hmftools.common.region.hmfslicer;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class HmfGenomeRegion implements GenomeRegion {
    @Override
    @NotNull
    @Value.Parameter
    public abstract String chromosome();

    @Override
    @Value.Parameter
    public abstract long start();

    @Override
    @Value.Parameter
    public abstract long end();

    @Override
    @Value.Default
    public String annotation() {
        return "";
    }

    @NotNull
    @Value.Parameter
    public abstract String transcriptID();

    @Value.Parameter
    public abstract int transcriptVersion();

    @NotNull
    @Value.Parameter
    public abstract String gene();

    @NotNull
    @Value.Parameter
    public abstract String chromosomeBand();

    @NotNull
    @Value.Parameter
    public abstract String entrezId();

    @NotNull
    public String transcript() {
        return transcriptID() + "." + transcriptVersion();
    }
}
