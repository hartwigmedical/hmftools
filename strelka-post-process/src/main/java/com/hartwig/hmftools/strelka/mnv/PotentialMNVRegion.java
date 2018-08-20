package com.hartwig.hmftools.strelka.mnv;

import static com.hartwig.hmftools.extensions.samtools.VariantContextUtils.splitMultiAlleleVariant;

import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PotentialMNVRegion {
    public abstract String chromosome();

    public abstract int start();

    //MIVO: end position, non-inclusive
    public abstract int end();

    public abstract List<VariantContext> variants();

    public abstract List<PotentialMNV> mnvs();

    @Value.Derived
    public List<PotentialMNV> potentialMnvs() {
        return mnvs().stream().filter(potentialMNV -> potentialMNV.variants().size() > 1).collect(Collectors.toList());
    }

    @Value.Derived
    Set<Integer> gapPositions() {
        return potentialMnvs().stream().flatMap(potentialMNV -> potentialMNV.gapPositions().stream()).collect(Collectors.toSet());
    }

    @Value.Lazy
    Set<VariantContext> variantsInPotentialMnvs() {
        return potentialMnvs().stream().flatMap(potentialMNV -> potentialMNV.variants().stream()).collect(Collectors.toSet());
    }

    @NotNull
    static PotentialMNVRegion addVariant(@NotNull final PotentialMNVRegion region, @NotNull final VariantContext variant,
            final int gapSize) {
        if (region.equals(PotentialMNVRegion.empty())) {
            return fromVariant(variant);
        } else {
            final List<PotentialMNV> mnvs = addVariantToPotentialMnvs(region.mnvs(), variant, gapSize);
            final List<VariantContext> variants = Lists.newArrayList(region.variants());
            variants.add(variant);
            final int mnvEnd = Math.max(region.end(), variant.getStart() + variant.getReference().getBaseString().length());
            return ImmutablePotentialMNVRegion.of(region.chromosome(), region.start(), mnvEnd, variants, mnvs);
        }
    }

    @NotNull
    static PotentialMNVRegion addVariants(@NotNull final PotentialMNVRegion region, @NotNull final List<VariantContext> variants,
            final int gapSize) {
        PotentialMNVRegion updatedRegion = region;
        for (final VariantContext variant : variants) {
            updatedRegion = addVariant(updatedRegion, variant, gapSize);
        }
        return updatedRegion;
    }

    @NotNull
    static PotentialMNVRegion fromVariant(@NotNull final VariantContext variant) {
        final List<PotentialMNV> mnvs = addVariantToPotentialMnvs(Lists.newArrayList(), variant, Integer.MAX_VALUE);
        return ImmutablePotentialMNVRegion.of(variant.getContig(), variant.getStart(),
                variant.getStart() + variant.getReference().getBaseString().length(), Lists.newArrayList(variant), mnvs);
    }

    @NotNull
    private static List<PotentialMNV> addVariantToPotentialMnvs(@NotNull final List<PotentialMNV> mnvs,
            @NotNull final VariantContext variant, final int gapSize) {
        if (variant.getAlternateAlleles().size() > 1) {
            return addVariantsToPotentialMnvs(mnvs, splitMultiAlleleVariant(variant), gapSize);
        } else {
            return addVariantsToPotentialMnvs(mnvs, Lists.newArrayList(variant), gapSize);
        }
    }

    @NotNull
    private static List<PotentialMNV> addVariantsToPotentialMnvs(@NotNull final List<PotentialMNV> mnvs,
            @NotNull final List<VariantContext> variants, final int gapSize) {
        final List<PotentialMNV> updatedMnvs = Lists.newArrayList(mnvs);
        variants.forEach(variant -> {
            updatedMnvs.add(PotentialMNV.fromVariant(variant));
            mnvs.stream()
                    .filter(potentialMNV -> potentialMNV.chromosome().equals(variant.getContig())
                            && variant.getStart() - potentialMNV.end() <= gapSize && variant.getStart() - potentialMNV.end() >= 0)
                    .forEach(potentialMNV -> updatedMnvs.add(PotentialMNV.addVariant(potentialMNV, variant)));
        });
        updatedMnvs.sort(Comparator.comparing(PotentialMNV::start).thenComparing(PotentialMNV::end));
        return updatedMnvs;
    }

    @NotNull
    public static PotentialMNVRegion empty() {
        return ImmutablePotentialMNVRegion.of("", -1, -1, Lists.newArrayList(), Lists.newArrayList());
    }
}
