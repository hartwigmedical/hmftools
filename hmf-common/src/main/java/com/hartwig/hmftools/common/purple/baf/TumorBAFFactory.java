package com.hartwig.hmftools.common.purple.baf;

import java.util.List;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.GermlineSampleData;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public class TumorBAFFactory {

    private static final Set<String> HETEROZYGOUS_GENO_TYPES = Sets.newHashSet("0/1", "0|1");

    private final double minRefAlleleFrequency;
    private final double maxRefAlleleFrequency;
    private final long minCombinedDepth;
    private final long maxCombinedDepth;

    public TumorBAFFactory(final double minRefAlleleFrequency, final double maxRefAlleleFrequency, final long minCombinedDepth,
            final long maxCombinedDepth) {
        this.minRefAlleleFrequency = minRefAlleleFrequency;
        this.maxRefAlleleFrequency = maxRefAlleleFrequency;
        this.minCombinedDepth = minCombinedDepth;
        this.maxCombinedDepth = maxCombinedDepth;
    }

    public Multimap<String, TumorBAF> createBAF(@NotNull List<GermlineVariant> variants) {
        final Multimap<String, TumorBAF> result = ArrayListMultimap.create();
        for (final GermlineVariant variant : variants) {
            if (eligible(variant)) {
                result.put(variant.chromosome(), create(variant));
            }
        }

        return result;
    }

    private boolean eligible(@NotNull final GermlineVariant variant) {
        final GermlineSampleData tumorData = variant.tumorData();
        return tumorData != null
                && variant.type() == VariantType.SNP
                && isHetrozygous(variant) 
                && inRefAlleleFrequencyRange(variant)
                && inCombinedDepthRange(variant);
    }

    private boolean isHetrozygous(@NotNull final GermlineVariant variant) {
        return HETEROZYGOUS_GENO_TYPES.contains(variant.refData().genoType());
    }

    private boolean inRefAlleleFrequencyRange(@NotNull final GermlineVariant variant) {
        return variant.refData().alleleFrequency() > minRefAlleleFrequency && variant.refData().alleleFrequency() < maxRefAlleleFrequency;
    }

    private boolean inCombinedDepthRange(@NotNull final GermlineVariant variant) {
        return variant.refData().combinedDepth() > minCombinedDepth && variant.refData().combinedDepth() < maxCombinedDepth;
    }

    private TumorBAF create(final GermlineVariant variant) {
        final GermlineSampleData tumorData = variant.tumorData();
        assert (tumorData != null);

        double baf = tumorData.alleleFrequency();
        return ImmutableTumorBAF.builder().from(variant).baf(baf).build();
    }
}
