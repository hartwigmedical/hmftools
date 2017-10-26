package com.hartwig.hmftools.strelka;

import static com.hartwig.hmftools.strelka.VariantContextUtils.splitMultiAlleleVariant;

import java.util.Comparator;
import java.util.List;
import java.util.Set;
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
public abstract class PotentialMNVRegion {
    private static final Logger LOGGER = LogManager.getLogger(PotentialMNVRegion.class);

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
    public Set<Integer> gapPositions() {
        return potentialMnvs().stream().flatMap(potentialMNV -> potentialMNV.gapPositions().stream()).collect(Collectors.toSet());
    }

    @NotNull
    static PotentialMNVRegion addVariant(@NotNull final PotentialMNVRegion region, @NotNull final VariantContext variant) {
        if (region.equals(PotentialMNVRegion.empty())) {
            return fromVariant(variant);
        } else {
            final List<PotentialMNV> mnvs = addVariantToPotentialMnvs(region.mnvs(), variant);
            final List<VariantContext> variants = Lists.newArrayList(region.variants());
            variants.add(variant);
            final int mnvEnd = Math.max(region.end(), variant.getStart() + variant.getReference().getBaseString().length());
            return ImmutablePotentialMNVRegion.of(region.chromosome(), region.start(), mnvEnd, variants, mnvs);
        }
    }

    @NotNull
    static PotentialMNVRegion addVariants(@NotNull final PotentialMNVRegion region, @NotNull final List<VariantContext> variants) {
        PotentialMNVRegion updatedRegion = region;
        for (final VariantContext variant : variants) {
            updatedRegion = addVariant(updatedRegion, variant);
        }
        return updatedRegion;
    }

    @NotNull
    static PotentialMNVRegion fromVariant(@NotNull final VariantContext variant) {
        final List<PotentialMNV> mnvs = addVariantToPotentialMnvs(Lists.newArrayList(), variant);
        return ImmutablePotentialMNVRegion.of(variant.getContig(), variant.getStart(),
                variant.getStart() + variant.getReference().getBaseString().length(), Lists.newArrayList(variant), mnvs);
    }

    @NotNull
    private static List<PotentialMNV> addVariantToPotentialMnvs(@NotNull final List<PotentialMNV> mnvs,
            @NotNull final VariantContext variant) {
        if (variant.getAlternateAlleles().size() > 1) {
            return addVariantsToPotentialMnvs(mnvs, splitMultiAlleleVariant(variant));
        } else {
            return addVariantsToPotentialMnvs(mnvs, Lists.newArrayList(variant));
        }
    }

    @NotNull
    private static List<PotentialMNV> addVariantsToPotentialMnvs(@NotNull final List<PotentialMNV> mnvs,
            @NotNull final List<VariantContext> variants) {
        final List<PotentialMNV> updatedMnvs = Lists.newArrayList(mnvs);
        variants.forEach(variant -> {
            updatedMnvs.add(PotentialMNV.fromVariant(variant));
            mnvs.stream()
                    .filter(potentialMNV -> potentialMNV.chromosome().equals(variant.getContig())
                            && variant.getStart() - potentialMNV.end() <= 1 && variant.getStart() - potentialMNV.end() >= 0)
                    .forEach(potentialMNV -> updatedMnvs.add(PotentialMNV.addVariant(potentialMNV, variant)));
        });
        updatedMnvs.sort(Comparator.comparing(PotentialMNV::start).thenComparing(PotentialMNV::end));
        return updatedMnvs;
    }

    @NotNull
    static PotentialMNVRegion empty() {
        return ImmutablePotentialMNVRegion.of("", -1, -1, Lists.newArrayList(), Lists.newArrayList());
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder(
                "chromosome " + chromosome() + ":\t" + (potentialMnvs().stream().allMatch(PotentialMNV::containsOnlySNVs) ? "SNV" : "INDEL")
                        + "[" + start() + " - " + end() + "]: ");
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
