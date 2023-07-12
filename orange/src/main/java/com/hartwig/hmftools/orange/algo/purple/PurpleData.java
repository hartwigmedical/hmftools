package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;

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
    List<PurpleSomaticVariant> allSomaticVariants();

    @NotNull
    List<PurpleSomaticVariant> reportableSomaticVariants();

    @Nullable
    List<PurpleSomaticVariant> allGermlineVariants();

    @Nullable
    List<PurpleSomaticVariant> reportableGermlineVariants();

    @NotNull
    List<StructuralVariant> allSomaticStructuralVariants();

    @Nullable
    List<StructuralVariant> allGermlineStructuralVariants();

    @NotNull
    List<PurpleCopyNumber> allSomaticCopyNumbers();

    @NotNull
    List<GeneCopyNumber> allSomaticGeneCopyNumbers();

    @Nullable
    List<GermlineDeletion> allGermlineDeletions();

    @Nullable
    List<GermlineDeletion> reportableGermlineDeletions();

}
