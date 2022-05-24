package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import org.junit.Test;

public class DriverCatalogMapTest {

    @Test
    public void canMapDriverCatalog() {
        String gene1 = "CDKN2A";
        String gene2 = "BRAF";
        double likelihood = 0.6;

        DriverCatalog driverGene1 =
                DriverCatalogTestFactory.createCanonicalSomaticMutationEntryForGene(gene1, likelihood, "transcript1", DriverCategory.ONCO);
        DriverCatalog driverGene2 =
                DriverCatalogTestFactory.createCanonicalSomaticMutationEntryForGene(gene1, likelihood, "transcript2", DriverCategory.ONCO);
        DriverCatalog driverGene3 =
                DriverCatalogTestFactory.createCanonicalSomaticMutationEntryForGene(gene2, likelihood, "transcript3", DriverCategory.ONCO);
        List<DriverCatalog> mergedDriverCatalog = Lists.newArrayList(driverGene1, driverGene2, driverGene3);

        Map<DriverCatalogKey, DriverCatalog> driverMap = DriverCatalogMap.toDriverMap(mergedDriverCatalog);

        assertEquals(driverGene1, driverMap.get(DriverCatalogKey.create(gene1, "transcript1")));
        assertEquals(driverGene2, driverMap.get(DriverCatalogKey.create(gene1, "transcript2")));
        assertEquals(driverGene3, driverMap.get(DriverCatalogKey.create(gene2, "transcript3")));
    }
}