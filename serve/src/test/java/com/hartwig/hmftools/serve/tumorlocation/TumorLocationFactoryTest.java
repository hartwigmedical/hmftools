package com.hartwig.hmftools.serve.tumorlocation;

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
        Set<TumorLocation> tumorLocationBlacklistings = Sets.newHashSet();
        tumorLocationBlacklistings.add(tumorLocationBlacklisting("Hematologic cancer", "2531"));
        assertEquals("Hematologic cancer,2531", TumorLocationFactory.extractTumorLocationBlacklisting(tumorLocationBlacklistings));
    }

    @Test
    public void canTestTumorLocationBlacklistingTwo() {
        Set<TumorLocation> tumorLocationBlacklistings = Sets.newHashSet();
        tumorLocationBlacklistings.add(tumorLocationBlacklisting("Hematologic cancer", "2531"));
        tumorLocationBlacklistings.add(tumorLocationBlacklisting("Skin Melanoma", "8923"));

        assertEquals("Hematologic cancer,2531;Skin Melanoma,8923",
                TumorLocationFactory.extractTumorLocationBlacklisting(tumorLocationBlacklistings));

    }

    @Test
    public void canTestTumorLocationBlacklistingMultiple() {
        Set<TumorLocation> tumorLocationBlacklistings = Sets.newHashSet();
        tumorLocationBlacklistings.add(tumorLocationBlacklisting("Hematologic cancer", "2531"));
        tumorLocationBlacklistings.add(tumorLocationBlacklisting("Skin Melanoma", "8923"));
        tumorLocationBlacklistings.add(tumorLocationBlacklisting("Bladder Cancer", "11054"));
        tumorLocationBlacklistings.add(tumorLocationBlacklisting("Colorectal Cancer", "1520"));

        assertEquals("Hematologic cancer,2531;Colorectal Cancer,1520;Skin Melanoma,8923;Bladder Cancer,11054",
                TumorLocationFactory.extractTumorLocationBlacklisting(tumorLocationBlacklistings));
    }

    @Test
    public void canTestReadTumorLocationBlacklistingSingle() {
        String combinedTumorLocations = "Hematologic cancer,2531";
        Set<TumorLocation> tumorLocationBlacklistings = TumorLocationFactory.readTumorLocationBlacklistingString(combinedTumorLocations);
        List<TumorLocation> tumorLocationBlacklistingsList = new ArrayList<>(tumorLocationBlacklistings);

        assertEquals(1, tumorLocationBlacklistingsList.size());
        TumorLocation tumorLocationBlacklisting = tumorLocationBlacklistingsList.get(0);
        assertEquals("Hematologic cancer", tumorLocationBlacklisting.cancerType());
        assertEquals("2531", tumorLocationBlacklisting.doid());
    }

    @Test
    public void canTestReadTumorLocationBlacklistingTwo() {
        String combinedTumorLocations = "Hematologic cancer,2531;Skin Melanoma,8923";
        Set<TumorLocation> tumorLocationBlacklistings = TumorLocationFactory.readTumorLocationBlacklistingString(combinedTumorLocations);
        List<TumorLocation> tumorLocationBlacklistingsList = new ArrayList<>(tumorLocationBlacklistings);

        assertEquals(2, tumorLocationBlacklistingsList.size());
        TumorLocation tumorLocationBlacklisting1 = tumorLocationBlacklistingsList.get(0);
        assertEquals("Hematologic cancer", tumorLocationBlacklisting1.cancerType());
        assertEquals("2531", tumorLocationBlacklisting1.doid());

        TumorLocation tumorLocationBlacklisting2 = tumorLocationBlacklistingsList.get(1);
        assertEquals("Skin Melanoma", tumorLocationBlacklisting2.cancerType());
        assertEquals("8923", tumorLocationBlacklisting2.doid());
    }

    @Test
    public void canTestReadTumorLocationBlacklistingMultiple() {
        String combinedTumorLocations = "Hematologic cancer,2531;Skin Melanoma,8923;Bladder Cancer,11054;Colorectal Cancer,1520";
        Set<TumorLocation> tumorLocationBlacklistings = TumorLocationFactory.readTumorLocationBlacklistingString(combinedTumorLocations);
        List<TumorLocation> tumorLocationBlacklistingsList = new ArrayList<>(tumorLocationBlacklistings);

        assertEquals(4, tumorLocationBlacklistingsList.size());
        TumorLocation tumorLocationBlacklisting1 = tumorLocationBlacklistingsList.get(0);
        assertEquals("Hematologic cancer", tumorLocationBlacklisting1.cancerType());
        assertEquals("2531", tumorLocationBlacklisting1.doid());

        TumorLocation tumorLocationBlacklisting4 = tumorLocationBlacklistingsList.get(1);
        assertEquals("Colorectal Cancer", tumorLocationBlacklisting4.cancerType());
        assertEquals("1520", tumorLocationBlacklisting4.doid());

        TumorLocation tumorLocationBlacklisting2 = tumorLocationBlacklistingsList.get(2);
        assertEquals("Skin Melanoma", tumorLocationBlacklisting2.cancerType());
        assertEquals("8923", tumorLocationBlacklisting2.doid());

        TumorLocation tumorLocationBlacklisting3 = tumorLocationBlacklistingsList.get(3);
        assertEquals("Bladder Cancer", tumorLocationBlacklisting3.cancerType());
        assertEquals("11054", tumorLocationBlacklisting3.doid());
    }

    @NotNull
    public TumorLocation tumorLocationBlacklisting(@NotNull String tumorLocation, @NotNull String doid) {
        return ImmutableTumorLocation.builder().cancerType(tumorLocation).doid(doid).build();
    }
}