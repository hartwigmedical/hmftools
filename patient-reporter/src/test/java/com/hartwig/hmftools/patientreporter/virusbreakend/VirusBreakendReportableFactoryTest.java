package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.*;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.virusbreakend.ImmutableVirusBreakend;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakend;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakendQCStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VirusBreakendReportableFactoryTest {

    @NotNull
    private static List<VirusBreakend> virusBreakendData() {

        List<VirusBreakend> virusBreakends = Lists.newArrayList();

        virusBreakends.add(ImmutableVirusBreakend.builder()
                .taxidGenus(0)
                .nameGenus(Strings.EMPTY)
                .readsGenusTree(0)
                .taxidSpecies(1)
                .nameSpecies(Strings.EMPTY)
                .readsSpeciesTree(0)
                .taxidAssigned(0)
                .nameAssigned("Human papillomavirus type 16")
                .readsAssignedTree(0)
                .readsAssignedDirect(0)
                .reference(Strings.EMPTY)
                .referenceTaxid(1)
                .referenceKmerCount(0)
                .alternateKmerCount(0)
                .Rname(Strings.EMPTY)
                .startpos(0)
                .endpos(0)
                .numreads(0)
                .covbases(0)
                .coverage(0)
                .meandepth(0)
                .meanbaseq(0)
                .meanmapq(0)
                .integrations(2)
                .QCStatus(VirusBreakendQCStatus.ASSEMBLY_DOWNSAMPLED)
                .build());

        virusBreakends.add(ImmutableVirusBreakend.builder()
                .taxidGenus(0)
                .nameGenus(Strings.EMPTY)
                .readsGenusTree(0)
                .taxidSpecies(0)
                .nameSpecies(Strings.EMPTY)
                .readsSpeciesTree(0)
                .taxidAssigned(0)
                .nameAssigned("Human papillomavirus type 16")
                .readsAssignedTree(0)
                .readsAssignedDirect(0)
                .reference(Strings.EMPTY)
                .referenceTaxid(0)
                .referenceKmerCount(0)
                .alternateKmerCount(0)
                .Rname(Strings.EMPTY)
                .startpos(0)
                .endpos(0)
                .numreads(0)
                .covbases(0)
                .coverage(0)
                .meandepth(0)
                .meanbaseq(0)
                .meanmapq(0)
                .integrations(2)
                .QCStatus(VirusBreakendQCStatus.LOW_VIRAL_COVERAGE)
                .build());

        virusBreakends.add(ImmutableVirusBreakend.builder()
                .taxidGenus(0)
                .nameGenus(Strings.EMPTY)
                .readsGenusTree(0)
                .taxidSpecies(0)
                .nameSpecies(Strings.EMPTY)
                .readsSpeciesTree(0)
                .taxidAssigned(1)
                .nameAssigned("Human papillomavirus type 16")
                .readsAssignedTree(0)
                .readsAssignedDirect(0)
                .reference(Strings.EMPTY)
                .referenceTaxid(0)
                .referenceKmerCount(0)
                .alternateKmerCount(0)
                .Rname(Strings.EMPTY)
                .startpos(0)
                .endpos(0)
                .numreads(0)
                .covbases(0)
                .coverage(0)
                .meandepth(0)
                .meanbaseq(0)
                .meanmapq(0)
                .integrations(2)
                .QCStatus(VirusBreakendQCStatus.ASSEMBLY_DOWNSAMPLED)
                .build());

        virusBreakends.add(ImmutableVirusBreakend.builder()
                .taxidGenus(0)
                .nameGenus(Strings.EMPTY)
                .readsGenusTree(0)
                .taxidSpecies(0)
                .nameSpecies(Strings.EMPTY)
                .readsSpeciesTree(0)
                .taxidAssigned(0)
                .nameAssigned("Human papillomavirus type 16")
                .readsAssignedTree(0)
                .readsAssignedDirect(0)
                .reference(Strings.EMPTY)
                .referenceTaxid(1)
                .referenceKmerCount(0)
                .alternateKmerCount(0)
                .Rname(Strings.EMPTY)
                .startpos(0)
                .endpos(0)
                .numreads(0)
                .covbases(0)
                .coverage(0)
                .meandepth(0)
                .meanbaseq(0)
                .meanmapq(0)
                .integrations(0)
                .QCStatus(VirusBreakendQCStatus.ASSEMBLY_DOWNSAMPLED)
                .build());

        return virusBreakends;
    }

    @Test
    public void canInterpretVirusBreakendForReportingPos() {
        List<VirusBreakend> virusBreakends = virusBreakendData();

        Map<Integer, String> VirusDBIdMap = Maps.newHashMap();
        VirusDBIdMap.put(1, "Human papillomavirus type 16");
        VirusDbModel virusDbModel = new VirusDbModel(VirusDBIdMap);

        Map<Integer, String> VirusIdSummaryMap = Maps.newHashMap();
        VirusIdSummaryMap.put(1, "EBV");
        VirusIdSummaryMap.put(2, "HPV");
        VirusSummaryModel virusSummaryModel = new VirusSummaryModel(VirusIdSummaryMap);

        assertEquals(2,
                VirusBreakendReportableFactory.analyzeVirusBreakend(virusBreakends, virusDbModel, virusSummaryModel)
                        .reportableVirussen()
                        .size());

        ReportableVirusbreakend reportableVirusbreakend =
                VirusBreakendReportableFactory.analyzeVirusBreakend(virusBreakends, virusDbModel, virusSummaryModel)
                        .reportableVirussen()
                        .get(0);
        assertEquals("Human papillomavirus type 16", reportableVirusbreakend.virusName());
        assertEquals(2, reportableVirusbreakend.integrations());

        assertEquals("EBV positive, EBV negative",
                VirusBreakendReportableFactory.analyzeVirusBreakend(virusBreakends, virusDbModel, virusSummaryModel).virusNameSummary());

    }
}