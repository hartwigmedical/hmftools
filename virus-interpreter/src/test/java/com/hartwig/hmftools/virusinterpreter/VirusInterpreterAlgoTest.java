package com.hartwig.hmftools.virusinterpreter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virus.VirusInterpretation;
import com.hartwig.hmftools.common.virus.VirusTestFactory;
import com.hartwig.hmftools.virusinterpreter.algo.ImmutableVirusWhitelist;
import com.hartwig.hmftools.virusinterpreter.algo.VirusBlacklistModel;
import com.hartwig.hmftools.virusinterpreter.algo.VirusWhitelist;
import com.hartwig.hmftools.virusinterpreter.algo.VirusWhitelistModel;
import com.hartwig.hmftools.virusinterpreter.taxonomy.TaxonomyDb;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VirusInterpreterAlgoTest {

    private static final String RUN_DIRECTORY = Resources.getResource("genomic").getPath();
    private static final String TUMOR_SAMPLE_WGS_METRICS = RUN_DIRECTORY + "/tumor_sample.wgsmetrics";
    private static final String PURPLE_PURITY_TSV = RUN_DIRECTORY + "/sample.purple.purity.tsv";
    private static final String PURPLE_QC_FILE = RUN_DIRECTORY + "/sample.purple.qc";

    @Test
    public void canAnalyzeVirusBreakends() throws IOException {
        List<VirusBreakend> virusBreakends = createTestVirusBreakends();

        VirusWhitelist virusWhitelist1 = ImmutableVirusWhitelist.builder()
                .taxidSpecies(1)
                .reportOnSummary(true)
                .virusInterpretation(VirusInterpretation.EBV)
                .nameSpecies("Human gammaherpesvirus 4")
                .integratedMinimalCoverage(null)
                .nonintegratedMinimalCoverage(90)
                .build();

        VirusWhitelist virusWhitelist2 = ImmutableVirusWhitelist.builder()
                .taxidSpecies(2)
                .reportOnSummary(true)
                .virusInterpretation(VirusInterpretation.EBV)
                .nameSpecies("Human gammaherpesvirus 4")
                .integratedMinimalCoverage(null)
                .nonintegratedMinimalCoverage(null)
                .build();

        String name = "Human papillomavirus type 16";
        Map<Integer, String> taxonomyMap = Maps.newHashMap();
        taxonomyMap.put(1, name);
        TaxonomyDb taxonomyDb = new TaxonomyDb(taxonomyMap);

        Map<Integer, VirusWhitelist> virusInterpretationMap = Maps.newHashMap();
        virusInterpretationMap.put(1, virusWhitelist1);
        virusInterpretationMap.put(2, virusWhitelist2);
        VirusWhitelistModel virusWhitelistModel = new VirusWhitelistModel(virusInterpretationMap);

        VirusBlacklistModel virusBlacklistModel = new VirusBlacklistModel(Sets.newHashSet(1), Sets.newHashSet());

        VirusInterpreterAlgo algo = new VirusInterpreterAlgo(taxonomyDb, virusWhitelistModel, virusBlacklistModel);
        List<AnnotatedVirus> annotatedViruses = algo.analyze(virusBreakends, PURPLE_PURITY_TSV, PURPLE_QC_FILE, TUMOR_SAMPLE_WGS_METRICS);
        assertEquals(5, annotatedViruses.size());
        assertEquals(2, annotatedViruses.stream().filter(x -> x.reported()).count());

        List<AnnotatedVirus> reportedVirus = Lists.newArrayList();
        for (AnnotatedVirus virus : annotatedViruses) {
            if (virus.reported()) {
                reportedVirus.add(virus);
            }
        }

        AnnotatedVirus reportedVirus1 = reportedVirus.get(0);
        assertEquals(name, reportedVirus1.name());
        assertEquals(2, reportedVirus1.integrations());
        assertEquals(VirusInterpretation.EBV, reportedVirus1.interpretation());

        AnnotatedVirus reportedVirus2 = reportedVirus.get(1);
        assertEquals(name, reportedVirus2.name());
        assertEquals(0, reportedVirus2.integrations());
        assertEquals(VirusInterpretation.EBV, reportedVirus2.interpretation());

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

        // This one has no integrations -- reported
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