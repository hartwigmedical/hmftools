package com.hartwig.hmftools.patientreporter.civic;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.patientreporter.report.data.Alteration;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;

public interface AlterationAnalyzer {
    List<Alteration> run(@NotNull final List<VariantReport> reportedVariants, @NotNull final List<GeneCopyNumber> copyNumbers,
            @NotNull final GeneModel geneModel, @NotNull final Set<String> tumorDoids);
}
