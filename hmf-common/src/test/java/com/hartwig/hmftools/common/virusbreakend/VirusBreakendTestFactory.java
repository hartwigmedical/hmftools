package com.hartwig.hmftools.common.virusbreakend;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class VirusBreakendTestFactory {

    private VirusBreakendTestFactory() {
    }

    @NotNull
    public static ImmutableVirusBreakend.Builder testBuilder() {
        return ImmutableVirusBreakend.builder()
                .taxidGenus(0)
                .nameGenus(Strings.EMPTY)
                .readsGenusTree(0)
                .taxidSpecies(0)
                .nameSpecies(Strings.EMPTY)
                .readsSpeciesTree(0)
                .taxidAssigned(0)
                .nameAssigned(Strings.EMPTY)
                .readsAssignedTree(0)
                .readsAssignedDirect(0)
                .reference(Strings.EMPTY)
                .referenceTaxid(0)
                .referenceKmerCount(0)
                .alternateKmerCount(0)
                .RName(Strings.EMPTY)
                .startPos(0)
                .endPos(0)
                .numReads(0)
                .covBases(0)
                .coverage(0)
                .meanDepth(0)
                .meanBaseQ(0)
                .meanMapQ(0)
                .integrations(0)
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES);
    }
}
