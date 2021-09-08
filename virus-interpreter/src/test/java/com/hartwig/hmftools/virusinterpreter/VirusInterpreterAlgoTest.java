package com.hartwig.hmftools.virusinterpreter;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virus.VirusTestFactory;
import com.hartwig.hmftools.virusinterpreter.algo.ImmutableVirusReporting;
import com.hartwig.hmftools.virusinterpreter.algo.VirusReporting;
import com.hartwig.hmftools.virusinterpreter.algo.VirusReportingModel;
import com.hartwig.hmftools.virusinterpreter.taxonomy.TaxonomyDb;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VirusInterpreterAlgoTest {

    private static final String GENOMIC_DIR = Resources.getResource("genomic").getPath();

    private static final String PURPLE_QC_FILE = GENOMIC_DIR + File.separator + "sample.purple.qc";
    private static final String PURPLE_PURITY_TSV = GENOMIC_DIR + File.separator + "sample.purple.purity.tsv";
    private static final String TUMOR_SAMPLE_WGS_METRICS = GENOMIC_DIR + File.separator + "sample.wgsmetrics";

    @Test
    public void canAnalyzeVirusBreakends() throws IOException {
        List<VirusBreakend> virusBreakends = createTestVirusBreakends();

        VirusReporting virusReporting1 = ImmutableVirusReporting.builder()
                .reportOnSummary(true)
                .virusInterpretation("EBV")
                .integratedMinimalCoverage(null)
                .nonIntegratedMinimalCoverage(90)
                .build();

        VirusReporting virusReporting2 = ImmutableVirusReporting.builder()
                .reportOnSummary(true)
                .virusInterpretation("EBV")
                .integratedMinimalCoverage(null)
                .nonIntegratedMinimalCoverage(null)
                .build();

        String name = "Human papillomavirus type 16";
        Map<Integer, String> taxonomyMap = Maps.newHashMap();
        taxonomyMap.put(1, name);
        TaxonomyDb taxonomyDb = new TaxonomyDb(taxonomyMap);

        Map<Integer, VirusReporting> virusReportingMap = Maps.newHashMap();
        virusReportingMap.put(1, virusReporting1);
        virusReportingMap.put(2, virusReporting2);

        VirusReportingModel virusWhitelistModel = new VirusReportingModel(virusReportingMap);

        VirusInterpreterAlgo algo = new VirusInterpreterAlgo(taxonomyDb, virusWhitelistModel);
        List<AnnotatedVirus> annotatedViruses = algo.analyze(virusBreakends, PURPLE_PURITY_TSV, PURPLE_QC_FILE, TUMOR_SAMPLE_WGS_METRICS);
        assertEquals(5, annotatedViruses.size());
        assertEquals(3, annotatedViruses.stream().filter(x -> x.reported()).count());

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
        assertEquals(2, reportedVirus2.integrations());
        assertEquals("EBV", reportedVirus2.interpretation());
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

        // This one has a blacklisted genus taxid -- not reported
        virusBreakends.add(VirusTestFactory.testVirusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(1)
                .taxidSpecies(1)
                .integrations(2)
                .coverage(0)
                .build());

        // This one has a failed QC -- not reported
        virusBreakends.add(VirusTestFactory.testVirusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(1)
                .qcStatus(VirusBreakendQCStatus.LOW_VIRAL_COVERAGE)
                .integrations(2)
                .coverage(0)
                .build());

        // This one has no integrations -- reported
        virusBreakends.add(VirusTestFactory.testVirusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(1)
                .integrations(0)
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

        return virusBreakends;
    }
}