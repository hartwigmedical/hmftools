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
    @NotNull
    PurityContext purityContext();

    @NotNull
    List<DriverCatalog> somaticDrivers();

    @Nullable
    List<DriverCatalog> germlineDrivers();

    @NotNull
    List<PurpleVariantContext> allSomaticVariants();

    @NotNull
    List<PurpleVariantContext> reportableSomaticVariants();

    @Nullable
    List<PurpleVariantContext> allGermlineVariants();

    @Nullable
    List<PurpleVariantContext> reportableGermlineVariants();

    @NotNull
    List<StructuralVariant> allPassingSomaticStructuralVariants();

    @Nullable
    List<StructuralVariant> allPassingGermlineStructuralVariants();

    @NotNull
    List<StructuralVariant> allInferredSomaticStructuralVariants();

    @Nullable
    List<StructuralVariant> allInferredGermlineStructuralVariants();

    @NotNull
    List<PurpleCopyNumber> allSomaticCopyNumbers();

    @NotNull
    List<GeneCopyNumber> allSomaticGeneCopyNumbers();

    @Nullable
    List<GermlineDeletion> allGermlineDeletions();

    @Nullable
    List<GermlineDeletion> reportableGermlineDeletions();

    @NotNull
    List<Segment> segments();
}
