package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;

@Value.Immutable
@Value.Style(passAnnotations = {NotNull.class, Nullable.class})
public abstract class PurpleInterpretedData {

    @NotNull
    public abstract PurpleFit fit();

    @NotNull
    public abstract PurpleCharacteristics characteristics();

    @NotNull
    public abstract List<DriverCatalog> somaticDrivers();

    @Nullable
    public abstract List<DriverCatalog> germlineDrivers();

    @NotNull
    public abstract List<PurpleVariant> allSomaticVariants();

    @NotNull
    public abstract List<PurpleVariant> reportableSomaticVariants();

    @NotNull
    public abstract List<PurpleVariant> additionalSuspectSomaticVariants();

    @Nullable
    public abstract List<PurpleVariant> allGermlineVariants();

    @Nullable
    public abstract List<PurpleVariant> reportableGermlineVariants();

    @Nullable
    public abstract List<PurpleVariant> additionalSuspectGermlineVariants();

    @NotNull
    public abstract List<PurpleCopyNumber> allSomaticCopyNumbers();

    @NotNull
    public abstract List<GeneCopyNumber> allSomaticGeneCopyNumbers();

    @NotNull
    public abstract List<GeneCopyNumber> suspectGeneCopyNumbersWithLOH();

    @NotNull
    public abstract List<PurpleGainLoss> allSomaticGainsLosses();

    @NotNull
    public abstract List<PurpleGainLoss> reportableSomaticGainsLosses();

    @NotNull
    public abstract List<PurpleGainLoss> nearReportableSomaticGains();

    @NotNull
    public abstract List<PurpleGainLoss> additionalSuspectSomaticGainsLosses();

    @Nullable
    public abstract List<GermlineDeletion> allGermlineDeletions();

    @Nullable
    public abstract List<PurpleGainLoss> allGermlineFullLosses();

    @Nullable
    public abstract List<PurpleGainLoss> reportableGermlineFullLosses();
}
