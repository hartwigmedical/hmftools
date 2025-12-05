package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.sv.StructuralVariant;

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

    List<PurpleVariantContext> allSomaticVariants();

    List<PurpleVariantContext> driverSomaticVariants();

    @Nullable
    List<PurpleVariantContext> allGermlineVariants();

    @Nullable
    List<PurpleVariantContext> driverGermlineVariants();

    List<StructuralVariant> allPassingSomaticStructuralVariants();

    @Nullable
    List<StructuralVariant> allPassingGermlineStructuralVariants();

    List<StructuralVariant> allInferredSomaticStructuralVariants();

    @Nullable
    List<StructuralVariant> allInferredGermlineStructuralVariants();

    List<PurpleCopyNumber> somaticCopyNumbers();

    List<GeneCopyNumber> somaticGeneCopyNumbers();

    @Nullable
    List<GermlineDeletion> allGermlineDeletions();

    @Nullable
    List<GermlineDeletion> driverGermlineDeletions();

    List<Segment> segments();
}
