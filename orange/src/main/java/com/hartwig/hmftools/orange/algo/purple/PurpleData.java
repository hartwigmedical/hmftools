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

    List<PurpleVariantContext> otherSomaticVariants();

    @Nullable
    List<PurpleVariantContext> allGermlineVariants();

    @Nullable
    List<PurpleVariantContext> driverGermlineVariants();

    @Nullable
    List<PurpleVariantContext> otherGermlineVariants();

    List<StructuralVariant> allPassingSomaticStructuralVariants();

    @Nullable
    List<StructuralVariant> allPassingGermlineStructuralVariants();

    List<StructuralVariant> allInferredSomaticStructuralVariants();

    @Nullable
    List<StructuralVariant> allInferredGermlineStructuralVariants();

    // copy number of every segment in the genome
    List<PurpleCopyNumber> somaticCopyNumbers();

    // copy number of every gene
    List<GeneCopyNumber> somaticGeneCopyNumbers();

    // driver gene gain deletions
    List<DriverCatalog> driverSomaticGainDels();

    @Nullable
    List<GermlineDeletion> allGermlineDeletions();

    @Nullable
    List<GermlineDeletion> driverGermlineDeletions();

    @Nullable
    List<GermlineDeletion> otherGermlineDeletions();

    List<Segment> segments();
}
