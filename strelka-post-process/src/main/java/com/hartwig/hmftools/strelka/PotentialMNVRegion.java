package com.hartwig.hmftools.strelka;

import static com.hartwig.hmftools.strelka.VariantContextUtils.splitMultiAlleleVariant;

import java.util.Comparator;
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
public abstract class PotentialMNVRegion {
    private static final Logger LOGGER = LogManager.getLogger(PotentialMNVRegion.class);

    public abstract String chromosome();

    public abstract int start();

    //MIVO: end position, non-inclusive
    public abstract int end();

    public abstract List<PotentialMNV> mnvs();

    @Value.Derived
    public List<PotentialMNV> potentialMnvs() {
        return mnvs().stream().filter(potentialMNV -> potentialMNV.variants().size() > 1).collect(Collectors.toList());
    }

    static PotentialMNVRegion addVariant(@NotNull final PotentialMNVRegion region, @NotNull final VariantContext variant) {
        if (region.equals(PotentialMNVRegion.empty())) {
            return fromVariant(variant);
        } else if (variant.getAlternateAlleles().size() > 1) {
            return addVariants(region, splitMultiAlleleVariant(variant));
        } else {
            final List<PotentialMNV> variants = Lists.newArrayList(region.mnvs());
            variants.add(PotentialMNV.fromVariant(variant));
            for (final PotentialMNV potentialMNV : region.mnvs()) {
                if (potentialMNV.chromosome().equals(variant.getContig()) && variant.getStart() - potentialMNV.end() <= 1
                        && variant.getStart() - potentialMNV.end() >= 0) {
                    variants.add(PotentialMNV.addVariant(potentialMNV, variant));
                }
            }
            variants.sort(Comparator.comparing(PotentialMNV::start).thenComparing(PotentialMNV::end));
            final int mnvEnd = Math.max(region.end(), variant.getStart() + variant.getReference().getBaseString().length());
            return ImmutablePotentialMNVRegion.of(region.chromosome(), region.start(), mnvEnd, variants);
        }
    }

    static PotentialMNVRegion addVariants(@NotNull final PotentialMNVRegion region, @NotNull final List<VariantContext> variants) {
        PotentialMNVRegion updatedRegion = region;
        for (final VariantContext variant : variants) {
            updatedRegion = addVariant(updatedRegion, variant);
        }
        return updatedRegion;
    }

    @NotNull
    static PotentialMNVRegion fromVariant(@NotNull final VariantContext variant) {
        final List<PotentialMNV> mnvs = Lists.newArrayList(PotentialMNV.fromVariant(variant));
        return ImmutablePotentialMNVRegion.of(variant.getContig(), variant.getStart(),
                variant.getStart() + variant.getReference().getBaseString().length(), mnvs);
    }

        }
    }

    static PotentialMNVRegion empty() {
        return ImmutablePotentialMNVRegion.of("", -1, -1, Lists.newArrayList());
    }
}
