package com.hartwig.hmftools.serve.cancertype;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import static org.junit.Assert.assertEquals;

import com.beust.jcommander.internal.Sets;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class TumorLocationFactoryTest {

    @Test
    public void canTestTumorLocationBlacklistingSingle() {
        Set<CancerType> cancerTypeBlacklist = Sets.newHashSet();
        cancerTypeBlacklist.add(tumorLocationBlacklisting("Hematologic cancer", "2531"));
        assertEquals("Hematologic cancer,2531", CancerTypeFactory.extractCancerTypeBlacklist(cancerTypeBlacklist));
    }

    @Test
    public void canTestTumorLocationBlacklistingTwo() {
        Set<CancerType> cancerTypeBlacklist = Sets.newHashSet();
        cancerTypeBlacklist.add(tumorLocationBlacklisting("Hematologic cancer", "2531"));
        cancerTypeBlacklist.add(tumorLocationBlacklisting("Skin Melanoma", "8923"));

        assertEquals("Hematologic cancer,2531;Skin Melanoma,8923",
                CancerTypeFactory.extractCancerTypeBlacklist(cancerTypeBlacklist));

    }

    @Test
    public void canTestTumorLocationBlacklistingMultiple() {
        Set<CancerType> cancerTypeBlacklist = Sets.newHashSet();
        cancerTypeBlacklist.add(tumorLocationBlacklisting("Hematologic cancer", "2531"));
        cancerTypeBlacklist.add(tumorLocationBlacklisting("Skin Melanoma", "8923"));
        cancerTypeBlacklist.add(tumorLocationBlacklisting("Bladder Cancer", "11054"));
        cancerTypeBlacklist.add(tumorLocationBlacklisting("Colorectal Cancer", "1520"));

        assertEquals("Hematologic cancer,2531;Colorectal Cancer,1520;Skin Melanoma,8923;Bladder Cancer,11054",
                CancerTypeFactory.extractCancerTypeBlacklist(cancerTypeBlacklist));
    }

    @Test
    public void canTestReadTumorLocationBlacklistingSingle() {
        String combinedTumorLocations = "Hematologic cancer,2531";
        Set<CancerType> cancerTypeBlacklistSet = CancerTypeFactory.readCancerTypeBlacklistString(combinedTumorLocations);
        List<CancerType> cancerTypeBlacklistList = new ArrayList<>(cancerTypeBlacklistSet);

        assertEquals(1, cancerTypeBlacklistList.size());
        CancerType tumorLocationBlacklisting = cancerTypeBlacklistList.get(0);
        assertEquals("Hematologic cancer", tumorLocationBlacklisting.cancerType());
        assertEquals("2531", tumorLocationBlacklisting.doid());
    }

    @Test
    public void canTestReadTumorLocationBlacklistingTwo() {
        String combinedTumorLocations = "Hematologic cancer,2531;Skin Melanoma,8923";
        Set<CancerType> cancerTypeBlacklistSet = CancerTypeFactory.readCancerTypeBlacklistString(combinedTumorLocations);
        List<CancerType> cancerTypeBlacklistList = new ArrayList<>(cancerTypeBlacklistSet);

        assertEquals(2, cancerTypeBlacklistList.size());
        CancerType tumorLocationBlacklisting1 = cancerTypeBlacklistList.get(0);
        assertEquals("Hematologic cancer", tumorLocationBlacklisting1.cancerType());
        assertEquals("2531", tumorLocationBlacklisting1.doid());

        CancerType tumorLocationBlacklisting2 = cancerTypeBlacklistList.get(1);
        assertEquals("Skin Melanoma", tumorLocationBlacklisting2.cancerType());
        assertEquals("8923", tumorLocationBlacklisting2.doid());
    }

    @Test
    public void canTestReadTumorLocationBlacklistingMultiple() {
        String combinedTumorLocations = "Hematologic cancer,2531;Skin Melanoma,8923;Bladder Cancer,11054;Colorectal Cancer,1520";
        Set<CancerType> cancerTypeBlacklistSet = CancerTypeFactory.readCancerTypeBlacklistString(combinedTumorLocations);
        List<CancerType> cancerTypeBlacklistList = new ArrayList<>(cancerTypeBlacklistSet);

        assertEquals(4, cancerTypeBlacklistList.size());
        CancerType tumorLocationBlacklisting1 = cancerTypeBlacklistList.get(0);
        assertEquals("Hematologic cancer", tumorLocationBlacklisting1.cancerType());
        assertEquals("2531", tumorLocationBlacklisting1.doid());

        CancerType tumorLocationBlacklisting4 = cancerTypeBlacklistList.get(1);
        assertEquals("Colorectal Cancer", tumorLocationBlacklisting4.cancerType());
        assertEquals("1520", tumorLocationBlacklisting4.doid());

        CancerType tumorLocationBlacklisting2 = cancerTypeBlacklistList.get(2);
        assertEquals("Skin Melanoma", tumorLocationBlacklisting2.cancerType());
        assertEquals("8923", tumorLocationBlacklisting2.doid());

        CancerType tumorLocationBlacklisting3 = cancerTypeBlacklistList.get(3);
        assertEquals("Bladder Cancer", tumorLocationBlacklisting3.cancerType());
        assertEquals("11054", tumorLocationBlacklisting3.doid());
    }

    @NotNull
    public CancerType tumorLocationBlacklisting(@NotNull String tumorLocation, @NotNull String doid) {
        return ImmutableCancerType.builder().cancerType(tumorLocation).doid(doid).build();
    }
}