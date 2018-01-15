package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.filter;
import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.passOnly;

import java.util.List;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantConsequence;

import org.jetbrains.annotations.NotNull;

public class VariantAnalyzer {

    @NotNull
    private final ConsequenceDeterminer determiner;

    @NotNull
    public static VariantAnalyzer fromSlicingRegions(@NotNull final GeneModel geneModel) {
        return new VariantAnalyzer(fromHmfSlicingRegion(geneModel));
    }

    @NotNull
    private static ConsequenceDeterminer fromHmfSlicingRegion(@NotNull final GeneModel geneModel) {
        return new ConsequenceDeterminer(geneModel);
    }

    private VariantAnalyzer(@NotNull final ConsequenceDeterminer determiner) {
        this.determiner = determiner;
    }

    @NotNull
    public VariantAnalysis run(@NotNull final List<SomaticVariant> variants) {
        final List<SomaticVariant> passedVariants = passOnly(variants);
        final List<SomaticVariant> missenseVariants = filter(passedVariants, isMissense());

        final ConsequenceOutput consequenceOutput = determiner.run(passedVariants);

        return ImmutableVariantAnalysis.of(variants, passedVariants, missenseVariants,
                consequenceOutput.consequentialVariants(),
                consequenceOutput.findings());
    }

    @NotNull
    private static Predicate<SomaticVariant> isMissense() {
        return variant -> variant.hasConsequence(VariantConsequence.MISSENSE_VARIANT);
    }
}
