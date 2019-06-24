package com.hartwig.hmftools.common.drivercatalog;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class PromoterDrivers {

    private static final String GENE = "TERT";
    private static final String CANONICAL_EFFECT = "upstream gene variant";
    private static final Set<Long> PROMOTER_REGIONS = Sets.newHashSet(1295242L, 1295228L, 1295250L);

    static public <T extends SomaticVariant> List<DriverCatalog> drivers(@NotNull final List<T> variants) {
        final List<DriverCatalog> result = Lists.newArrayList();
        if (variants.stream().anyMatch(PromoterDrivers::tertPromoter)) {
            DriverCatalog driver = ImmutableDriverCatalog.builder()
                    .gene(GENE)
                    .missense(0)
                    .nonsense(0)
                    .inframe(0)
                    .frameshift(0)
                    .splice(0)
                    .dndsLikelihood(0)
                    .driverLikelihood(1)
                    .driver(DriverType.PROMOTER)
                    .category(DriverCategory.ONCO)
                    .build();
            result.add(driver);
        }

        return result;
    }

    private static boolean tertPromoter(@NotNull final SomaticVariant variant) {
        return variant.gene().equals(GENE) && PROMOTER_REGIONS.contains(variant.position()) && variant.canonicalEffect()
                .equals(CANONICAL_EFFECT);
    }

}
