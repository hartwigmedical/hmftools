package com.hartwig.hmftools.virusinterpreter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virus.VirusTestFactory;
import com.hartwig.hmftools.virusinterpreter.algo.ImmutableVirusReportingDb;
import com.hartwig.hmftools.virusinterpreter.algo.VirusReportingDb;
import com.hartwig.hmftools.virusinterpreter.algo.VirusReportingDbModel;
import com.hartwig.hmftools.virusinterpreter.coverages.CoveragesAnalysis;
import com.hartwig.hmftools.virusinterpreter.coverages.ImmutableCoveragesAnalysis;
import com.hartwig.hmftools.virusinterpreter.taxonomy.TaxonomyDb;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VirusInterpreterAlgoTest {

    @Test
    public void canAnalyzeVirusBreakends() {
        String name = "Human gammaherpesvirus 4";
        List<VirusBreakend> virusBreakends = createTestVirusBreakends();

        VirusInterpreterAlgo algo = createTestAlgo(name);

        List<AnnotatedVirus> annotatedViruses = algo.analyze(virusBreakends);
        assertEquals(7, annotatedViruses.size());
        assertEquals(4, annotatedViruses.stream().filter(x -> x.reported()).count());

        List<AnnotatedVirus> reportedVirus = Lists.newArrayList();
        for (AnnotatedVirus virus : annotatedViruses) {
            if (virus.reported()) {
                reportedVirus.add(virus);
            }
        }

        AnnotatedVirus reportedVirus1 = reportedVirus.get(0);
        assertEquals(name, reportedVirus1.name());
        assertEquals(2, reportedVirus1.integrations());
        assertEquals("EBV", reportedVirus1.interpretation());

        AnnotatedVirus reportedVirus2 = reportedVirus.get(1);
        assertEquals(name, reportedVirus2.name());
        assertEquals(0, reportedVirus2.integrations());
        assertEquals("EBV", reportedVirus2.interpretation());

        AnnotatedVirus reportedVirus3 = reportedVirus.get(2);
        assertEquals(name, reportedVirus3.name());
        assertEquals(0, reportedVirus3.integrations());
        assertEquals("EBV", reportedVirus3.interpretation());

        AnnotatedVirus reportedVirus4 = reportedVirus.get(3);
        assertEquals(name, reportedVirus4.name());
        assertEquals(0, reportedVirus4.integrations());
        assertEquals("EBV", reportedVirus4.interpretation());

        assertEquals(Integer.valueOf(90), algo.determineMinimalCoverageVirus(0, 1));
        assertNull(algo.determineMinimalCoverageVirus(1, 1));
    }

    private static VirusInterpreterAlgo createTestAlgo(@NotNull String name) {
        VirusReportingDb virusReporting1 = ImmutableVirusReportingDb.builder()
                .virusInterpretation("EBV")
                .integratedMinimalCoverage(null)
                .nonIntegratedMinimalCoverage(90)
                .build();

        VirusReportingDb virusReporting2 = ImmutableVirusReportingDb.builder()
                .virusInterpretation("EBV")
                .integratedMinimalCoverage(null)
                .nonIntegratedMinimalCoverage(null)
                .build();

        Map<Integer, String> taxonomyMap = Maps.newHashMap();
        taxonomyMap.put(1, name);
        TaxonomyDb taxonomyDb = new TaxonomyDb(taxonomyMap);

        Map<Integer, VirusReportingDb> virusReportingMap = Maps.newHashMap();
        virusReportingMap.put(1, virusReporting1);
        virusReportingMap.put(2, virusReporting2);

        VirusReportingDbModel virusReportingModel = new VirusReportingDbModel(virusReportingMap);

        CoveragesAnalysis coveragesAnalysis = ImmutableCoveragesAnalysis.builder().expectedClonalCoverage(34.5).build();

        return new VirusInterpreterAlgo(taxonomyDb, virusReportingModel, coveragesAnalysis);
    }

    @NotNull
    private static List<VirusBreakend> createTestVirusBreakends() {
        List<VirusBreakend> virusBreakends = Lists.newArrayList();

        // This one should be added --reported
        virusBreakends.add(VirusTestFactory.testVirusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(1)
                .integrations(2)
                .coverage(0)
                .build());

        // This virus not present in reporting model -- not reported
        virusBreakends.add(VirusTestFactory.testVirusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(1)
                .taxidSpecies(4)
                .integrations(2)
                .coverage(0)
                .build());
        //
        // This one has a failed QC -- not reported
        virusBreakends.add(VirusTestFactory.testVirusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(1)
                .qcStatus(VirusBreakendQCStatus.LOW_VIRAL_COVERAGE)
                .integrations(2)
                .coverage(0)
                .build());
        //
        // This one has no integrations -- reported
        virusBreakends.add(VirusTestFactory.testVirusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(1)
                .integrations(0)
                .meanDepth(35)
                .coverage(91)
                .build());

        // This one has no integrations -- not reported
        virusBreakends.add(VirusTestFactory.testVirusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(1)
                .integrations(0)
                .coverage(90)
                .build());

        // This one has no integrations -- reported
        virusBreakends.add(VirusTestFactory.testVirusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(2)
                .integrations(0)
                .coverage(91)
                .build());

        // This one has no integrations -- reported
        virusBreakends.add(VirusTestFactory.testVirusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(2)
                .integrations(0)
                .coverage(90)
                .build());

        return virusBreakends;
    }
}