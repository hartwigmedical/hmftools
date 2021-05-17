package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.assertEquals;

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

    @Test
    public void canInterpretVirusBreakendForReportingPos() {
        List<VirusBreakend> virusBreakends = virusBreakendData();

        Map<Integer, String> taxonomyMap = Maps.newHashMap();
        taxonomyMap.put(1, "Human papillomavirus type 16");
        TaxonomyDb taxonomyDb = new TaxonomyDb(taxonomyMap);

        Map<Integer, String> virusInterpretationMap = Maps.newHashMap();
        virusInterpretationMap.put(1, "EBV");
        virusInterpretationMap.put(2, "HPV");
        VirusInterpretationModel virusInterpretationModel = new VirusInterpretationModel(virusInterpretationMap);

        Map<Integer, String> virusBlacklistMap = Maps.newHashMap();
        virusBlacklistMap.put(1, "taxid_genus");
        virusBlacklistMap.put(2, "HPV");
        VirusBlacklistModel virusBlacklistModel = new VirusBlacklistModel(virusBlacklistMap);

        assertEquals(1,
                VirusBreakendReportableFactory.analyzeVirusBreakend("",
                        virusBreakends, taxonomyDb, virusInterpretationModel,
                        virusBlacklistModel).reportableViruses().size());

        ReportableVirusBreakend reportableVirusbreakend = VirusBreakendReportableFactory.analyzeVirusBreakend("",
                virusBreakends, taxonomyDb, virusInterpretationModel,
                virusBlacklistModel).reportableViruses().get(0);
        assertEquals("Human papillomavirus type 16", reportableVirusbreakend.virusName());
        assertEquals(2, reportableVirusbreakend.integrations());

        assertEquals("EBV positive, HPV negative",
                VirusBreakendReportableFactory.analyzeVirusBreakend("",
                        virusBreakends, taxonomyDb, virusInterpretationModel,
                        virusBlacklistModel).virusNameSummary());
    }

    @NotNull
    private static List<VirusBreakend> virusBreakendData() {
        List<VirusBreakend> virusBreakends = Lists.newArrayList();

        virusBreakends.add(testBuilder().nameAssigned("Human papillomavirus type 16")
                .taxidGenus(0)
                .integrations(2)
                .build());

        virusBreakends.add(testBuilder().nameAssigned("Human papillomavirus type 16")
                .taxidGenus(1)
                .integrations(2)
                .build());

        virusBreakends.add(testBuilder().nameAssigned("Human papillomavirus type 16")
                .taxidGenus(0)
                .qcStatus(VirusBreakendQCStatus.LOW_VIRAL_COVERAGE)
                .integrations(2)
                .build());

        virusBreakends.add(testBuilder().nameAssigned("Human papillomavirus type 16")
                .taxidGenus(0)
                .integrations(0)
                .build());

        return virusBreakends;
    }

    @NotNull
    private static ImmutableVirusBreakend.Builder testBuilder() {
        return ImmutableVirusBreakend.builder()
                .taxidGenus(1)
                .nameGenus(Strings.EMPTY)
                .readsGenusTree(0)
                .taxidSpecies(1)
                .nameSpecies(Strings.EMPTY)
                .readsSpeciesTree(0)
                .taxidAssigned(1)
                .nameAssigned(Strings.EMPTY)
                .readsAssignedTree(0)
                .readsAssignedDirect(0)
                .reference(Strings.EMPTY)
                .referenceTaxid(1)
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
                .integrations(2)
                .qcStatus(VirusBreakendQCStatus.PASS);
    }
}