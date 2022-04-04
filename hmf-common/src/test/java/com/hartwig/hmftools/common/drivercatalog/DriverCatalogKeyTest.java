package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.junit.Test;

public class DriverCatalogKeyTest {

    @Test
    public void canCreateUniqueKeySetMap() {
        String gene1 = "CDKN2A";
        String gene2 = "BRAF";
        double likelihood = 0.6;

        DriverCatalog driverGene1 = DriverCatalogMapTest.createCanonicalSomaticMutationEntryForGene(gene1, likelihood, "transcript1");
        DriverCatalog driverGene2 = DriverCatalogMapTest.createCanonicalSomaticMutationEntryForGene(gene1, likelihood, "transcript2");
        DriverCatalog driverGene3 = DriverCatalogMapTest.createCanonicalSomaticMutationEntryForGene(gene2, likelihood, "transcript3");
        List<DriverCatalog> mergedDriverCatalog = Lists.newArrayList(driverGene1, driverGene2, driverGene3);

        Map<DriverCatalogKey, DriverCatalog> driverMap = DriverCatalogMap.toDriverMap(mergedDriverCatalog);

        Set<DriverCatalogKey> keys = Sets.newHashSet(DriverCatalogKey.create(gene1, "transcript1"),
                DriverCatalogKey.create(gene1, "transcript2"),
                DriverCatalogKey.create(gene2, "transcript3"));

        assertEquals(keys, DriverCatalogKey.buildUniqueKeysSet(driverMap));
    }

    @Test
    public void canCreateUniqueKeySet() {
        String gene1 = "CDKN2A";
        String gene2 = "BRAF";
        double likelihood = 0.6;

        DriverCatalog driverGene1 = DriverCatalogMapTest.createCanonicalSomaticMutationEntryForGene(gene1, likelihood, "transcript1");
        DriverCatalog driverGene2 = DriverCatalogMapTest.createCanonicalSomaticMutationEntryForGene(gene1, likelihood, "transcript2");
        DriverCatalog driverGene3 = DriverCatalogMapTest.createCanonicalSomaticMutationEntryForGene(gene2, likelihood, "transcript3");
        List<DriverCatalog> mergedDriverCatalog = Lists.newArrayList(driverGene1, driverGene2, driverGene3);

        Set<DriverCatalogKey> keys = Sets.newHashSet(DriverCatalogKey.create(gene1, "transcript1"),
                DriverCatalogKey.create(gene1, "transcript2"),
                DriverCatalogKey.create(gene2, "transcript3"));
        assertEquals(keys, DriverCatalogKey.buildUniqueKeysSet(mergedDriverCatalog));
    }
}