package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakendTestFactory;

import org.junit.Test;

public class VirusBlacklistModelTest {

    @Test
    public void canMatchIdToVirus() {
        VirusBlacklistModel virusBlacklistModel = new VirusBlacklistModel(Sets.newHashSet(1), Sets.newHashSet(2));

        VirusBreakend filteredGenus = VirusBreakendTestFactory.testBuilder().taxidGenus(1).build();
        assertTrue(virusBlacklistModel.isBlacklisted(filteredGenus));

        VirusBreakend filteredSpecies = VirusBreakendTestFactory.testBuilder().taxidSpecies(2).build();
        assertTrue(virusBlacklistModel.isBlacklisted(filteredSpecies));

        VirusBreakend pass = VirusBreakendTestFactory.testBuilder().taxidSpecies(1).taxidGenus(2).build();
        assertFalse(virusBlacklistModel.isBlacklisted(pass));
    }
}