package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Optional;

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
    default List<PurpleVariantContext> reportableSomaticVariants()
    {
        return allSomaticVariants().stream().filter(PurpleVariantContext::reported).toList();
    }

    @Nullable
    List<PurpleVariantContext> allGermlineVariants();

    @Nullable
    default List<PurpleVariantContext> reportableGermlineVariants()
    {
        return Optional.ofNullable(allGermlineVariants())
                .map(o -> o.stream().filter(PurpleVariantContext::reported).toList())
                .orElse(null);
    }

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
    default List<GermlineDeletion> reportableGermlineDeletions()
    {
        return Optional.ofNullable(allGermlineDeletions())
                .map(o -> o.stream().filter(germlineDeletion -> germlineDeletion.Reported).toList())
                .orElse(null);
    }

    @NotNull
    List<Segment> segments();
}
