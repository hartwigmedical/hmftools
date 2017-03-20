package com.hartwig.hmftools.patientreporter.batch;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.patientreporter.algo.GenomeAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;

class RunStats {

    private final int allVariantCount;
    private final int passedVariantCount;
    private final int consensusPassedCount;
    private final int mutationalLoad;
    @NotNull
    private final Map<VariantConsequence, Integer> consequenceCounts;
    @NotNull
    private final List<VariantReport> variantFindings;
    @NotNull
    private final List<CopyNumberReport> copyNumberFindings;

    @NotNull
    static RunStats fromGenomeAnalysis(@NotNull GenomeAnalysis genomeAnalysis) {
        final VariantAnalysis variantAnalysis = genomeAnalysis.variantAnalysis();
        final CopyNumberAnalysis copyNumberAnalysis = genomeAnalysis.copyNumberAnalysis();

        return new RunStats(variantAnalysis.allVariants().size(), variantAnalysis.passedVariants().size(),
                variantAnalysis.consensusPassedVariants().size(), variantAnalysis.mutationalLoad(),
                ConsequenceCount.count(variantAnalysis.consensusPassedVariants()), variantAnalysis.findings(),
                copyNumberAnalysis.findings());
    }

    private RunStats(final int allVariantCount, final int passedVariantCount, final int consensusPassedCount,
            final int mutationalLoad, @NotNull final Map<VariantConsequence, Integer> consequenceCounts,
            @NotNull final List<VariantReport> variantFindings,
            @NotNull final List<CopyNumberReport> copyNumberFindings) {
        this.allVariantCount = allVariantCount;
        this.passedVariantCount = passedVariantCount;
        this.consensusPassedCount = consensusPassedCount;
        this.mutationalLoad = mutationalLoad;
        this.consequenceCounts = consequenceCounts;
        this.variantFindings = variantFindings;
        this.copyNumberFindings = copyNumberFindings;
    }

    int allVariantCount() {
        return allVariantCount;
    }

    int passedVariantCount() {
        return passedVariantCount;
    }

    int consensusPassedCount() {
        return consensusPassedCount;
    }

    int mutationalLoad() {
        return mutationalLoad;
    }

    @NotNull
    Map<VariantConsequence, Integer> consequenceCounts() {
        return consequenceCounts;
    }

    @NotNull
    List<VariantReport> variantFindings() {
        return variantFindings;
    }

    @NotNull
    List<CopyNumberReport> copyNumberFindings() {
        return copyNumberFindings;
    }
}
