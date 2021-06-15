package com.hartwig.hmftools.common.drivercatalog.dnds;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface DndsMutationalLoad
{

    @NotNull
    String sampleId();

    int snvBiallelic();

    int snvNonBiallelic();

    default int snvTotal()
    {
        return snvBiallelic() + snvNonBiallelic();
    }

    int indelBiallelic();

    int indelNonBiallelic();

    default int indelTotal()
    {
        return indelBiallelic() + indelNonBiallelic();
    }
}
