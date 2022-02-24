package com.hartwig.hmftools.serve.blacklisting;

import java.util.List;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class TumorLocationBlacklistTest {

    @Test
    public void canTestTumorLocationBlacklistingSingle() {
        List<TumorLocationBlacklisting> tumorLocationBlacklistings = Lists.newArrayList();
        tumorLocationBlacklistings.add(tumorLocationBlacklisting("Hematologic cancer", "2531"));
        assertEquals("Hematologic cancer,2531", TumorLocationBlacklist.extractTumorLocationBlacklisting(tumorLocationBlacklistings));
    }

    @Test
    public void canTestTumorLocationBlacklistingTwo() {
        List<TumorLocationBlacklisting> tumorLocationBlacklistings = Lists.newArrayList();
        tumorLocationBlacklistings.add(tumorLocationBlacklisting("Hematologic cancer", "2531"));
        tumorLocationBlacklistings.add(tumorLocationBlacklisting("Skin Melanoma", "8923"));

        assertEquals("Hematologic cancer,2531;Skin Melanoma,8923",
                TumorLocationBlacklist.extractTumorLocationBlacklisting(tumorLocationBlacklistings));

    }

    @Test
    public void canTestTumorLocationBlacklistingMultiple() {
        List<TumorLocationBlacklisting> tumorLocationBlacklistings = Lists.newArrayList();
        tumorLocationBlacklistings.add(tumorLocationBlacklisting("Hematologic cancer", "2531"));
        tumorLocationBlacklistings.add(tumorLocationBlacklisting("Skin Melanoma", "8923"));
        tumorLocationBlacklistings.add(tumorLocationBlacklisting("Bladder Cancer", "11054"));
        tumorLocationBlacklistings.add(tumorLocationBlacklisting("Colorectal Cancer", "1520"));

        assertEquals("Hematologic cancer,2531;Skin Melanoma,8923;Bladder Cancer,11054;Colorectal Cancer,1520",
                TumorLocationBlacklist.extractTumorLocationBlacklisting(tumorLocationBlacklistings));
    }

    @Test
    public void canTestReadTumorLocationBlacklistingSingle() {
        String combinedTumorLocations = "Hematologic cancer,2531";
        List<TumorLocationBlacklisting> tumorLocationBlacklistings = TumorLocationBlacklist.readTumorLocationBlacklistingString(combinedTumorLocations);

        assertEquals(1, tumorLocationBlacklistings.size());
        TumorLocationBlacklisting tumorLocationBlacklisting = tumorLocationBlacklistings.get(0);
        assertEquals("Hematologic cancer", tumorLocationBlacklisting.blacklistCancerType());
        assertEquals("2531", tumorLocationBlacklisting.blacklistedDoid());
    }

    @Test
    public void canTestReadTumorLocationBlacklistingTwo() {
        String combinedTumorLocations = "Hematologic cancer,2531;Skin Melanoma,8923";
        List<TumorLocationBlacklisting> tumorLocationBlacklistings = TumorLocationBlacklist.readTumorLocationBlacklistingString(combinedTumorLocations);

        assertEquals(2, tumorLocationBlacklistings.size());
        TumorLocationBlacklisting tumorLocationBlacklisting1 = tumorLocationBlacklistings.get(0);
        assertEquals("Hematologic cancer", tumorLocationBlacklisting1.blacklistCancerType());
        assertEquals("2531", tumorLocationBlacklisting1.blacklistedDoid());

        TumorLocationBlacklisting tumorLocationBlacklisting2 = tumorLocationBlacklistings.get(1);
        assertEquals("Skin Melanoma", tumorLocationBlacklisting2.blacklistCancerType());
        assertEquals("8923", tumorLocationBlacklisting2.blacklistedDoid());
    }

    @Test
    public void canTestReadTumorLocationBlacklistingMultiple() {
        String combinedTumorLocations = "Hematologic cancer,2531;Skin Melanoma,8923;Bladder Cancer,11054;Colorectal Cancer,1520";
        List<TumorLocationBlacklisting> tumorLocationBlacklistings = TumorLocationBlacklist.readTumorLocationBlacklistingString(combinedTumorLocations);

        assertEquals(4, tumorLocationBlacklistings.size());
        TumorLocationBlacklisting tumorLocationBlacklisting1 = tumorLocationBlacklistings.get(0);
        assertEquals("Hematologic cancer", tumorLocationBlacklisting1.blacklistCancerType());
        assertEquals("2531", tumorLocationBlacklisting1.blacklistedDoid());

        TumorLocationBlacklisting tumorLocationBlacklisting2 = tumorLocationBlacklistings.get(1);
        assertEquals("Skin Melanoma", tumorLocationBlacklisting2.blacklistCancerType());
        assertEquals("8923", tumorLocationBlacklisting2.blacklistedDoid());

        TumorLocationBlacklisting tumorLocationBlacklisting3 = tumorLocationBlacklistings.get(2);
        assertEquals("Bladder Cancer", tumorLocationBlacklisting3.blacklistCancerType());
        assertEquals("11054", tumorLocationBlacklisting3.blacklistedDoid());

        TumorLocationBlacklisting tumorLocationBlacklisting4 = tumorLocationBlacklistings.get(3);
        assertEquals("Colorectal Cancer", tumorLocationBlacklisting4.blacklistCancerType());
        assertEquals("1520", tumorLocationBlacklisting4.blacklistedDoid());
    }

    @NotNull
    public TumorLocationBlacklisting tumorLocationBlacklisting(@NotNull String tumorLocation, @NotNull String doid) {
        return ImmutableTumorLocationBlacklisting.builder().blacklistCancerType(tumorLocation).blacklistedDoid(doid).build();
    }
}