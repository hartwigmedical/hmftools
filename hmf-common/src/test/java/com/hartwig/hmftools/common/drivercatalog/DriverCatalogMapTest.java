package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class DriverCatalogMapTest  {

    @Test
    public void canMapDriverCatalog() {
        String gene1 = "CDKN2A";
        String gene2 = "BRAF";
        double likelihood = 0.6;

        DriverCatalog driverGene1 = createCanonicalSomaticMutationEntryForGene(gene1, likelihood, "transcript1");
        DriverCatalog driverGene2 = createCanonicalSomaticMutationEntryForGene(gene1, likelihood, "transcript2");
        DriverCatalog driverGene3 = createCanonicalSomaticMutationEntryForGene(gene2, likelihood, "transcript3");
        List<DriverCatalog> mergedDriverCatalog = Lists.newArrayList(driverGene1, driverGene2, driverGene3);


        DriverCatalogKey driverCatalogKey1 = DriverCatalogKey.create(gene1, "transcript1");
        DriverCatalogKey driverCatalogKey2 = DriverCatalogKey.create(gene1, "transcript2");
        DriverCatalogKey driverCatalogKey3 = DriverCatalogKey.create(gene2, "transcript3");

        Map<DriverCatalogKey, DriverCatalog> driverMap = DriverCatalogMap.toDriverMap(mergedDriverCatalog);

        assertEquals(driverGene1, driverMap.get(driverCatalogKey1));
        assertEquals(driverGene2, driverMap.get(driverCatalogKey2));
        assertEquals(driverGene3, driverMap.get(driverCatalogKey3));
    }

    @NotNull
    public static DriverCatalog createCanonicalSomaticMutationEntryForGene(@NotNull String gene, double likelihood,
            @NotNull String transcript) {
        return create(gene, likelihood, DriverType.MUTATION, transcript);
    }

    @NotNull
    public static DriverCatalog createCanonicalGermlineMutationEntryForGene(@NotNull String gene, double likelihood,
            @NotNull String transcript) {
        return create(gene, likelihood, DriverType.GERMLINE_MUTATION, transcript);
    }

    private static DriverCatalog create(@NotNull String gene, double likelihood, @NotNull DriverType type, @NotNull String transcript) {
        return ImmutableDriverCatalog.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(gene)
                .transcript(transcript)
                .isCanonical(true)
                .driver(type)
                .category(DriverCategory.ONCO)
                .likelihoodMethod(LikelihoodMethod.DNDS)
                .driverLikelihood(likelihood)
                .missense(0)
                .nonsense(0)
                .splice(0)
                .inframe(0)
                .frameshift(0)
                .biallelic(false)
                .minCopyNumber(0)
                .maxCopyNumber(0)
                .build();
    }
}