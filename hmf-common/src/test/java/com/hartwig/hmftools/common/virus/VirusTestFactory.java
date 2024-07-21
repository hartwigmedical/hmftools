package com.hartwig.hmftools.common.virus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class VirusTestFactory
{
    @NotNull
    public static VirusInterpreterData createEmptyData()
    {
        return ImmutableVirusInterpreterData.builder().build();
    }

    @NotNull
    public static VirusInterpreterData createMinimalData()
    {
        return createWithVirus(annotatedVirusBuilder().build());
    }

    @NotNull
    public static VirusInterpreterData createProperData()
    {
        return createWithVirus(annotatedVirusBuilder().interpretation(VirusType.MCV).expectedClonalCoverage(1.0).build());
    }

    @NotNull
    public static VirusInterpreterData createHHVInterpretationData()
    {
        return createWithVirus(annotatedVirusBuilder().interpretation(VirusType.HHV8).build());
    }

    @NotNull
    public static ImmutableVirusBreakend.Builder virusBreakendBuilder()
    {
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
                .meanDepth(0)
                .meanBaseQ(0)
                .meanMapQ(0)
                .integrations(0)
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES);
    }

    @NotNull
    public static ImmutableAnnotatedVirus.Builder annotatedVirusBuilder()
    {
        return ImmutableAnnotatedVirus.builder()
                .taxid(0)
                .name(Strings.EMPTY)
                .interpretation(null)
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                .integrations(0)
                .percentageCovered(1.0)
                .meanCoverage(1.0)
                .expectedClonalCoverage(1.0)
                .reported(false)
                .blacklisted(false)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH);
    }

    @NotNull
    private static VirusInterpreterData createWithVirus(@NotNull AnnotatedVirus virus)
    {
        return ImmutableVirusInterpreterData.builder().addReportableViruses(virus).addAllViruses(virus).build();
    }
}
