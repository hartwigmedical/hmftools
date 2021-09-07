package com.hartwig.hmftools.virusinterpreter.algo;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusTestFactory;

import org.junit.Test;

public class VirusBlacklistModelTest {

    @Test
    public void canMatchIdToVirus() {
        VirusBlacklistModel virusBlacklistModel = new VirusBlacklistModel(Sets.newHashSet(1), Sets.newHashSet(2));

        VirusBreakend filteredGenus = VirusTestFactory.testVirusBreakendBuilder().taxidGenus(1).coverage(0).build();
        assertTrue(virusBlacklistModel.isBlacklisted(filteredGenus));

        VirusBreakend filteredSpecies = VirusTestFactory.testVirusBreakendBuilder().taxidSpecies(2).coverage(0).build();
        assertTrue(virusBlacklistModel.isBlacklisted(filteredSpecies));

        VirusBreakend pass = VirusTestFactory.testVirusBreakendBuilder().taxidSpecies(1).taxidGenus(2).coverage(0).build();
        assertFalse(virusBlacklistModel.isBlacklisted(pass));
    }
}