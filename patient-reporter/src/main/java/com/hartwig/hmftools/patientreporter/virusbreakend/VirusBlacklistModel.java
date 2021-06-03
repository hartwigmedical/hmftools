package com.hartwig.hmftools.patientreporter.virusbreakend;

import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.virus.VirusBreakend;

import org.jetbrains.annotations.NotNull;

public class VirusBlacklistModel {

    @NotNull
    private final Set<Integer> blacklistedGenera;
    @NotNull
    private final Set<Integer> blacklistedSpecies;

    public VirusBlacklistModel(@NotNull final Set<Integer> blacklistedGenera, @NotNull final Set<Integer> blacklistedSpecies) {
        this.blacklistedGenera = blacklistedGenera;
        this.blacklistedSpecies = blacklistedSpecies;
    }

    public boolean isBlacklisted(@NotNull VirusBreakend virusBreakend) {
        return blacklistedGenera.contains(virusBreakend.taxidGenus()) || blacklistedSpecies.contains(virusBreakend.taxidSpecies());
    }

    @VisibleForTesting
    int genusCount() {
        return blacklistedGenera.size();
    }

    @VisibleForTesting
    int speciesCount() {
        return blacklistedSpecies.size();
    }
}
