package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;
import com.hartwig.hmftools.patientreporter.variants.ImmutableVariantReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleAnalysis {

    @NotNull
    public abstract Gender gender();

    @NotNull
    public abstract FittedPurityStatus status();

    @NotNull
    public abstract FittedPurity fittedPurity();

    @NotNull
    public abstract FittedPurityScore fittedScorePurity();

    @NotNull
    public abstract List<PurpleCopyNumber> copyNumbers();

    @NotNull
    public abstract List<GeneCopyNumber> panelGeneCopyNumbers();

    public int genePanelSize() {
        return panelGeneCopyNumbers().size();
    }

    @NotNull
    public List<GeneCopyNumber> reportableGeneCopyNumbers() {
        return PurpleCopyNumberFilter.filterCopyNumbersForReport(fittedPurity().ploidy(), panelGeneCopyNumbers());
    }

    public double purityUncertainty() {
        return fittedScorePurity().maxPurity() - fittedScorePurity().minPurity();
    }

    @NotNull
    public String purityString() {
        return status() == FittedPurityStatus.NO_TUMOR
                ? "[below detection threshold]"
                : PatientReportFormat.formatPercent(fittedPurity().purity());
    }

    @NotNull
    public List<VariantReport> enrich(@NotNull final List<VariantReport> variantReports) {
        final List<VariantReport> result = Lists.newArrayList();
        final PurityAdjuster purityAdjuster = new PurityAdjuster(gender(), fittedPurity());
        final GenomeRegionSelector<PurpleCopyNumber> copyNumberSelector = GenomeRegionSelectorFactory.create(copyNumbers());

        for (final VariantReport variantReport : variantReports) {
            final Optional<PurpleCopyNumber> optionalCopyNumber = copyNumberSelector.select(variantReport.variant());
            if (optionalCopyNumber.isPresent()) {
                final PurpleCopyNumber copyNumber = optionalCopyNumber.get();
                double adjustedVAF = Math.min(1,
                        purityAdjuster.purityAdjustedVAF(copyNumber.chromosome(),
                                copyNumber.averageTumorCopyNumber(),
                                variantReport.alleleFrequency()));
                result.add(ImmutableVariantReport.builder()
                        .from(variantReport)
                        .baf(copyNumber.descriptiveBAF())
                        .impliedVAF(adjustedVAF)
                        .build());
            } else {
                result.add(variantReport);
            }
        }

        return result;
    }
}
