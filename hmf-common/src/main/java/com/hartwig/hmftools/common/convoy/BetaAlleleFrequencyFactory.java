package com.hartwig.hmftools.common.convoy;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chromosome.Chromosomes;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class BetaAlleleFrequencyFactory {

    private static final Set<String> GENO_TYPE = Sets.newHashSet("0/1", "0|1");

    private final double minRefAlleleFrequency;
    private final double maxRefAlleleFrequency;
    private final long minCombinedDepth;
    private final long maxCombinedDepth;

    public BetaAlleleFrequencyFactory(double minRefAlleleFrequency, double maxRefAlleleFrequency, long minCombinedDepth, long maxCombinedDepth) {
        this.minRefAlleleFrequency = minRefAlleleFrequency;
        this.maxRefAlleleFrequency = maxRefAlleleFrequency;
        this.minCombinedDepth = minCombinedDepth;
        this.maxCombinedDepth = maxCombinedDepth;
    }

    public List<BetaAlleleFrequency> transform(List<GermlineVariant> variants) {
        return variants.stream().filter(this::filter).map(this::transform).collect(Collectors.toList());
    }

    private BetaAlleleFrequency transform(GermlineVariant variant) {
        int chrom = Chromosomes.asInt(variant.chromosome());
        double standardBAH = variant.tumorData().alleleFrequency();
        double modifiedBAF = 0.5 + Math.abs(standardBAH - 0.5);
        double chromPosition = chrom + variant.position() / (double) Chromosomes.length(variant.chromosome());
        return ImmutableBetaAlleleFrequency.of(chromPosition, standardBAH, modifiedBAF, variant.chromosome(), variant.position());
    }

    private boolean filter(GermlineVariant variant) {
        if (variant.tumorData() == null ||
                !GENO_TYPE.contains(variant.refData().genoType()) ||
                variant.type() != VariantType.SNP ||
                variant.refData().alleleFrequency() < minRefAlleleFrequency ||
                variant.refData().alleleFrequency() > maxRefAlleleFrequency ||
                variant.refData().combinedDepth() < minCombinedDepth ||
                variant.refData().combinedDepth() > maxCombinedDepth) {
            return false;
        }

        return true;
    }
}
