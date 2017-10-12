package com.hartwig.hmftools.strelka;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PotentialMNV {
    private static final Logger LOGGER = LogManager.getLogger(MNVDetector.class);

    public abstract String chromosome();

    public abstract int start();

    //MIVO: end position, non-inclusive
    public abstract int end();

    @Value.Derived
    public int lastPosition() {
        return end() - 1;
    }

    @Value.Derived
    public boolean containsOnlySNVs() {
        return variants().stream().allMatch(VariantContext::isSNP);
    }

    public abstract List<VariantContext> variants();

    public abstract List<Integer> gapPositions();

    static PotentialMNV addVariant(@NotNull final PotentialMNV potentialMnv, @NotNull final VariantContext variant) {
        if (potentialMnv.equals(PotentialMNV.empty())) {
            return fromVariant(variant);
        } else {
            final List<VariantContext> variants = Lists.newArrayList(potentialMnv.variants());
            variants.add(variant);
            final List<Integer> gaps = Lists.newArrayList(potentialMnv.gapPositions());
            if (potentialMnv.end() != variant.getStart()) {
                gaps.add(potentialMnv.end());
            }
            return ImmutablePotentialMNV.of(potentialMnv.chromosome(), potentialMnv.start(),
                    variant.getStart() + variant.getReference().getBaseString().length(), variants, gaps);
        }
    }

    static PotentialMNV fromVariant(@NotNull final VariantContext variant) {
        return ImmutablePotentialMNV.of(variant.getContig(), variant.getStart(),
                variant.getStart() + variant.getReference().getBaseString().length(), Lists.newArrayList(variant), Lists.newArrayList());
    }

    static PotentialMNV fromVariants(@NotNull final List<VariantContext> variants) {
        PotentialMNV state = PotentialMNV.empty();
        for (final VariantContext variant : variants) {
            state = addVariant(state, variant);
        }
        return state;
    }

    static PotentialMNV empty() {
        return ImmutablePotentialMNV.of("", -1, -1, Lists.newArrayList(), Lists.newArrayList());
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder(
                "chromosome " + chromosome() + ":\t" + (containsOnlySNVs() ? "SNV" : "INDEL") + "[" + start() + " - " + end() + "]: ");
        variants().forEach(variant -> {
            sb.append("[");
            sb.append(variant.getStart());
            sb.append(": ");
            sb.append(variant.getReference().getBaseString());
            sb.append(" -> ");
            sb.append(variant.getAlternateAlleles().stream().map(Allele::getBaseString).collect(Collectors.joining(",")));
            sb.append("] ");
        });
        sb.append("gaps: ");
        gapPositions().forEach(gap -> {
            sb.append(gap);
            sb.append(", ");
        });
        return sb.toString();
    }
}
