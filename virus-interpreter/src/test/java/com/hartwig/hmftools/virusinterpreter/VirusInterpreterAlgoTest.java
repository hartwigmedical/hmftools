package com.hartwig.hmftools.virusinterpreter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.common.virus.VirusTestFactory;
import com.hartwig.hmftools.virusinterpreter.algo.ImmutableVirusReportingDb;
import com.hartwig.hmftools.virusinterpreter.algo.VirusReportingDb;
import com.hartwig.hmftools.virusinterpreter.algo.VirusReportingDbModel;
import com.hartwig.hmftools.virusinterpreter.coverages.CoveragesAnalysis;
import com.hartwig.hmftools.virusinterpreter.taxonomy.TaxonomyDb;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VirusInterpreterAlgoTest
{
    private static final String GENOMIC_DIR = Resources.getResource("genomic").getPath();

    private static final String PURPLE_QC_FILE = GENOMIC_DIR + File.separator + "sample.purple.qc";
    private static final String PURPLE_PURITY_TSV = GENOMIC_DIR + File.separator + "sample.purple.purity.tsv";

    @Test
    public void canTestIsPurpleQCPass()
    {
        assertTrue(VirusInterpreterAlgo.hasAcceptablePurpleQuality(Sets.newHashSet(PurpleQCStatus.WARN_DELETED_GENES)));
        assertTrue(VirusInterpreterAlgo.hasAcceptablePurpleQuality(Sets.newHashSet(PurpleQCStatus.PASS)));
        assertFalse(VirusInterpreterAlgo.hasAcceptablePurpleQuality(Sets.newHashSet(PurpleQCStatus.FAIL_CONTAMINATION)));
        assertFalse(VirusInterpreterAlgo.hasAcceptablePurpleQuality(Sets.newHashSet(PurpleQCStatus.FAIL_NO_TUMOR)));
    }

    @Test
    public void canAnalyzeVirusBreakends() throws IOException
    {
        String name = "Human gammaherpesvirus 4";
        List<VirusBreakend> virusBreakends = createTestVirusBreakends();

        VirusInterpreterAlgo algo = createTestAlgo(name);
        PurityContext purityContext = PurityContextFile.readWithQC(PURPLE_QC_FILE, PURPLE_PURITY_TSV);

        List<AnnotatedVirus> annotatedViruses = algo.analyze(virusBreakends, purityContext);
        assertEquals(7, annotatedViruses.size());
        assertEquals(4, annotatedViruses.stream().filter(x -> x.reported()).count());

        List<AnnotatedVirus> reportedVirus = Lists.newArrayList();
        for(AnnotatedVirus virus : annotatedViruses)
        {
            if(virus.reported())
            {
                reportedVirus.add(virus);
            }
        }

        AnnotatedVirus reportedVirus1 = reportedVirus.get(0);
        assertEquals(name, reportedVirus1.name());
        assertEquals(2, reportedVirus1.integrations());
        assertEquals(VirusType.EBV, reportedVirus1.interpretation());

        AnnotatedVirus reportedVirus2 = reportedVirus.get(1);
        assertEquals(name, reportedVirus2.name());
        assertEquals(0, reportedVirus2.integrations());
        assertEquals(VirusType.EBV, reportedVirus2.interpretation());

        AnnotatedVirus reportedVirus3 = reportedVirus.get(2);
        assertEquals(name, reportedVirus3.name());
        assertEquals(0, reportedVirus3.integrations());
        assertEquals(VirusType.EBV, reportedVirus3.interpretation());

        AnnotatedVirus reportedVirus4 = reportedVirus.get(3);
        assertEquals(name, reportedVirus4.name());
        assertEquals(0, reportedVirus4.integrations());
        assertEquals(VirusType.EBV, reportedVirus4.interpretation());

        assertEquals(Integer.valueOf(90), algo.determineMinimalCoverageVirus(0, 1));
        assertNull(algo.determineMinimalCoverageVirus(1, 1));
    }

    @Test
    public void canDetermineHighRiskVirus()
    {
        VirusInterpreterAlgo algo = createTestAlgo();

        assertEquals(algo.virusLikelihoodType(createTestVirusBreakendsForHighRiskVirus(1), true), VirusLikelihoodType.HIGH);
        assertEquals(algo.virusLikelihoodType(createTestVirusBreakendsForHighRiskVirus(2), true), VirusLikelihoodType.HIGH);
        assertEquals(algo.virusLikelihoodType(createTestVirusBreakendsForHighRiskVirus(3), false), VirusLikelihoodType.UNKNOWN);
    }

    @Test
    public void canReportTest()
    {
        VirusInterpreterAlgo algo = createTestAlgo();

        assertTrue(algo.report(createTestVirusBreakendsFail(1).get(0), 1.0, Sets.newHashSet(PurpleQCStatus.FAIL_NO_TUMOR)));
        assertTrue(algo.report(createTestVirusBreakendsFail(1).get(0), 1.0, Sets.newHashSet(PurpleQCStatus.FAIL_CONTAMINATION)));
        assertFalse(algo.report(createTestVirusBreakendsFail(0).get(0), 1.0, Sets.newHashSet(PurpleQCStatus.FAIL_NO_TUMOR)));
        assertFalse(algo.report(createTestVirusBreakendsFail(0).get(0), 1.0, Sets.newHashSet(PurpleQCStatus.FAIL_CONTAMINATION)));
    }

    @Test
    public void canDetermineBlacklistedVirus()
    {
        VirusInterpreterAlgo algo = createTestAlgo();
        assertTrue(algo.blacklist(11646, 0));
        assertFalse(algo.blacklist(10000, 0));
    }

    @NotNull
    private static VirusInterpreterAlgo createTestAlgo()
    {
        return createTestAlgo(Strings.EMPTY);
    }

    @NotNull
    private static VirusInterpreterAlgo createTestAlgo(@NotNull String name)
    {
        VirusReportingDb virusReporting1 = ImmutableVirusReportingDb.builder()
                .virusInterpretation(VirusType.EBV)
                .integratedMinimalCoverage(null)
                .nonIntegratedMinimalCoverage(90)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH)
                .build();

        VirusReportingDb virusReporting2 = ImmutableVirusReportingDb.builder()
                .virusInterpretation(VirusType.EBV)
                .integratedMinimalCoverage(null)
                .nonIntegratedMinimalCoverage(null)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH)
                .build();

        VirusReportingDb virusReporting3 = ImmutableVirusReportingDb.builder()
                .virusInterpretation(VirusType.EBV)
                .integratedMinimalCoverage(null)
                .nonIntegratedMinimalCoverage(null)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH)
                .build();

        Map<Integer, String> taxonomyMap = Maps.newHashMap();
        taxonomyMap.put(1, name);
        TaxonomyDb taxonomyDb = new TaxonomyDb(taxonomyMap);

        Map<Integer, VirusReportingDb> virusReportingMap = Maps.newHashMap();
        virusReportingMap.put(1, virusReporting1);
        virusReportingMap.put(2, virusReporting2);
        virusReportingMap.put(3, virusReporting3);

        VirusReportingDbModel virusReportingModel = new VirusReportingDbModel(virusReportingMap);

        CoveragesAnalysis coveragesAnalysis = new CoveragesAnalysis(34.5);

        return new VirusInterpreterAlgo(taxonomyDb, Lists.newArrayList(11646), virusReportingModel, coveragesAnalysis);
    }

    private static VirusBreakend createTestVirusBreakendsForHighRiskVirus(int taxidSpecies)
    {

        return VirusTestFactory.virusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(taxidSpecies)
                .integrations(2)
                .coverage(0)
                .build();
    }

    private static List<VirusBreakend> createTestVirusBreakends()
    {
        List<VirusBreakend> virusBreakends = Lists.newArrayList();

        // This one should be added --reported
        virusBreakends.add(VirusTestFactory.virusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(1)
                .integrations(2)
                .coverage(0)
                .build());

        // This virus not present in reporting model -- not reported
        virusBreakends.add(VirusTestFactory.virusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(1)
                .taxidSpecies(4)
                .integrations(2)
                .coverage(0)
                .build());
        //
        // This one has a failed QC -- not reported
        virusBreakends.add(VirusTestFactory.virusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(1)
                .qcStatus(VirusBreakendQCStatus.LOW_VIRAL_COVERAGE)
                .integrations(2)
                .coverage(0)
                .build());
        //
        // This one has no integrations -- reported
        virusBreakends.add(VirusTestFactory.virusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(1)
                .integrations(0)
                .meanDepth(35)
                .coverage(91)
                .build());

        // This one has no integrations -- not reported
        virusBreakends.add(VirusTestFactory.virusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(1)
                .integrations(0)
                .coverage(90)
                .build());

        // This one has no integrations -- reported
        virusBreakends.add(VirusTestFactory.virusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(2)
                .integrations(0)
                .coverage(91)
                .build());

        // This one has no integrations -- reported
        virusBreakends.add(VirusTestFactory.virusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(2)
                .integrations(0)
                .coverage(90)
                .build());

        return virusBreakends;
    }

    private static List<VirusBreakend> createTestVirusBreakendsFail(int integrations)
    {
        List<VirusBreakend> virusBreakends = Lists.newArrayList();

        // This one should be added --reported
        virusBreakends.add(VirusTestFactory.virusBreakendBuilder()
                .referenceTaxid(1)
                .taxidGenus(2)
                .taxidSpecies(1)
                .integrations(integrations)
                .coverage(0)
                .build());
        return virusBreakends;
    }
}