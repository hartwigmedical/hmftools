package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleData
{
    PurityContext purityContext();

    List<DriverCatalog> somaticDrivers();

    @Nullable
    List<DriverCatalog> germlineDrivers();

    List<PurpleVariantContext> somaticVariants();

    @Nullable
    List<PurpleVariantContext> germlineVariants();

    List<PurpleCopyNumber> somaticCopyNumbers();

    List<GeneCopyNumber> somaticGeneCopyNumbers();

    @Nullable
    List<GermlineDeletion> germlineDeletions();
}
