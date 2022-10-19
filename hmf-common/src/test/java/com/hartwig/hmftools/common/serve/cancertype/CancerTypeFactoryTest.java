package com.hartwig.hmftools.common.serve.cancertype;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CancerTypeFactoryTest {

    @Test
    public void canConvertSingle() {
        Set<CancerType> cancerTypes = Sets.newHashSet();
        cancerTypes.add(create("Hematologic cancer", "2531"));
        assertEquals("Hematologic cancer,2531", CancerTypeFactory.toString(cancerTypes));
    }

    @Test
    public void canConvertTwo() {
        Set<CancerType> cancerTypes = Sets.newHashSet();
        cancerTypes.add(create("Hematologic cancer", "2531"));
        cancerTypes.add(create("Skin Melanoma", "8923"));

        assertEquals("Hematologic cancer,2531;Skin Melanoma,8923", CancerTypeFactory.toString(cancerTypes));
    }

    @Test
    public void canConvertMultiple() {
        Set<CancerType> cancerTypes = Sets.newHashSet();
        cancerTypes.add(create("Hematologic cancer", "2531"));
        cancerTypes.add(create("Skin Melanoma", "8923"));
        cancerTypes.add(create("Bladder Cancer", "11054"));
        cancerTypes.add(create("Colorectal Cancer", "1520"));

        assertEquals("Hematologic cancer,2531;Colorectal Cancer,1520;Skin Melanoma,8923;Bladder Cancer,11054",
                CancerTypeFactory.toString(cancerTypes));
    }

    @Test
    public void canReadSingle() {
        String combinedNameAndDoid = "Hematologic cancer,2531";
        Set<CancerType> cancerTypes = CancerTypeFactory.fromString(combinedNameAndDoid);

        assertEquals(1, cancerTypes.size());
        CancerType tumorLocationBlacklisting = cancerTypes.iterator().next();
        assertEquals("Hematologic cancer", tumorLocationBlacklisting.name());
        assertEquals("2531", tumorLocationBlacklisting.doid());
    }

    @Test
    public void canReadTwo() {
        String combinedNamesAndDoids = "Hematologic cancer,2531;Skin Melanoma,8923";
        Set<CancerType> cancerTypes = CancerTypeFactory.fromString(combinedNamesAndDoids);

        assertEquals(2, cancerTypes.size());
        CancerType cancerType1 = findByName(cancerTypes, "Hematologic cancer");
        assertEquals("2531", cancerType1.doid());

        CancerType cancerType2 = findByName(cancerTypes, "Skin Melanoma");
        assertEquals("8923", cancerType2.doid());
    }

    @Test
    public void canReadMultiple() {
        String combinedNamesAndDoids = "Hematologic cancer,2531;Skin Melanoma,8923;Bladder Cancer,11054;Colorectal Cancer,1520";
        Set<CancerType> cancerTypes = CancerTypeFactory.fromString(combinedNamesAndDoids);

        assertEquals(4, cancerTypes.size());
        CancerType cancerType1 = findByName(cancerTypes, "Hematologic cancer");
        assertEquals("2531", cancerType1.doid());

        CancerType cancerType2 = findByName(cancerTypes, "Colorectal Cancer");
        assertEquals("1520", cancerType2.doid());

        CancerType cancerType3 = findByName(cancerTypes, "Skin Melanoma");
        assertEquals("8923", cancerType3.doid());

        CancerType cancerType4 = findByName(cancerTypes, "Bladder Cancer");
        assertEquals("11054", cancerType4.doid());
    }

    @NotNull
    private static CancerType findByName(@NotNull Iterable<CancerType> cancerTypes, @NotNull String nameToFind) {
        for (CancerType cancerType : cancerTypes) {
            if (cancerType.name().equals(nameToFind)) {
                return cancerType;
            }
        }

        throw new IllegalStateException("Could not find cancerType with name: " + nameToFind);

    }

    @Test
    public void canCreateDoidStrings() {
        Set<CancerType> cancerTypes = Sets.newHashSet();
        cancerTypes.add(create("Hematologic cancer", "2531"));
        cancerTypes.add(create("Skin Melanoma", "8923"));
        Set<String> doids = CancerTypeFactory.doidStrings(cancerTypes);
        assertEquals(Sets.newHashSet("2531", "8923"), doids);
    }

    @NotNull
    private static CancerType create(@NotNull String name, @NotNull String doid) {
        return ImmutableCancerType.builder().name(name).doid(doid).build();
    }
}