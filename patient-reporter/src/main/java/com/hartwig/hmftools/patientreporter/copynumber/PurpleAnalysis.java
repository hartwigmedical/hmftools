package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.patientreporter.genepanel.GeneModel;

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
    public abstract List<GeneCopyNumber> geneCopyNumbers();

    @NotNull
    public abstract List<GeneCopyNumber> panelGeneCopyNumbers();

    @NotNull
    public List<GeneCopyNumber> reportableGeneCopyNumbers(@NotNull GeneModel panelGeneModel) {
        return ReportableCopyNumbers.filterCopyNumbersForReport(fittedPurity().ploidy(), panelGeneCopyNumbers(), panelGeneModel);
    }
}
