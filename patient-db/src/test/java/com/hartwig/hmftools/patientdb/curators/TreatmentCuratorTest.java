package com.hartwig.hmftools.patientdb.curators;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.patientdb.data.CuratedTreatment;

import org.junit.Test;

public class TreatmentCuratorTest {

    @Test
    public void canCreateFromProductionResource() throws IOException {
        assertNotNull(TreatmentCurator.fromProductionResource());
    }

    @Test
    public void matchesExactSingleWord() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final Optional<CuratedTreatment> curatedTreatment = curator.matchSingle("Zocor");
        assertTrue(curatedTreatment.isPresent());
        assertEquals("Zocor", curatedTreatment.get().name());
    }

    @Test
    public void removesSearchTermsFromUnused() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        assertEquals(22, curator.unusedSearchTerms().size());

        // KODU: Match with a canonical and make sure we can re-match after updating unused terms.
        assertEquals(1, curator.search("Zocor").size());
        assertEquals(21, curator.unusedSearchTerms().size());
        assertEquals(1, curator.search("Zocor").size());

        // KODU: Match with "other" name
        curator.search("amlodipine besylate");
        assertEquals(20, curator.unusedSearchTerms().size());

        // KODU: do an imperfect match but still remove the search term.
        curator.search("Avastine");
        assertEquals(19, curator.unusedSearchTerms().size());

        // KODU: Bogus example!!
        curator.search("This does not match at all!");
        assertEquals(19, curator.unusedSearchTerms().size());
    }

    @Test
    public void matchesIgnoringCaseSingleWord() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final Optional<CuratedTreatment> curatedTreatment = curator.matchSingle("zocor");
        assertTrue(curatedTreatment.isPresent());
        assertEquals("Zocor", curatedTreatment.get().name());
    }

    @Test
    public void matchesSingleWordWithTypo() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final Optional<CuratedTreatment> curatedTreatment = curator.matchSingle("lisinoprill");
        assertTrue(curatedTreatment.isPresent());
        assertEquals("Lisinopril", curatedTreatment.get().name());
    }

    @Test
    public void matchesExactMultiWord() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final Optional<CuratedTreatment> curatedTreatment = curator.matchSingle("amlodipine besylate");
        assertTrue(curatedTreatment.isPresent());
        assertEquals("Norvasc", curatedTreatment.get().name());
    }

    @Test
    public void matchesIgnoringCaseMultiWord() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final Optional<CuratedTreatment> curatedTreatment = curator.matchSingle("Amlodipine Besylate");
        assertTrue(curatedTreatment.isPresent());
        assertEquals("Norvasc", curatedTreatment.get().name());
    }

    @Test
    public void matchesMultiWordWithTypo() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final Optional<CuratedTreatment> curatedTreatment = curator.matchSingle("Amlodipin besylat");
        assertTrue(curatedTreatment.isPresent());
        assertEquals("Norvasc", curatedTreatment.get().name());
    }

    @Test
    public void matchesTermWithSpecialChars() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final Optional<CuratedTreatment> curatedTreatment = curator.matchSingle("z-pak");
        assertTrue(curatedTreatment.isPresent());
        assertEquals("Azithromycin", curatedTreatment.get().name());
    }

    @Test
    public void matchMultipleTreatments() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final List<CuratedTreatment> curatedTreatments = curator.matchMultiple("Prinivil,Zithromax/amlodipine besylate");
        assertEquals(3, curatedTreatments.size());
        final List<String> matches = curatedTreatments.stream().map(CuratedTreatment::name).collect(Collectors.toList());
        assertTrue(matches.contains("Lisinopril"));
        assertTrue(matches.contains("Norvasc"));
        assertTrue(matches.contains("Azithromycin"));
    }

    @Test
    public void matchMultipleTreatmentsWithTypos() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final List<CuratedTreatment> curatedTreatments = curator.matchMultiple("Prinivyl,Zithromaxx/amlodipin Besylate");
        assertEquals(3, curatedTreatments.size());
        final List<String> matches = curatedTreatments.stream().map(CuratedTreatment::name).collect(Collectors.toList());
        assertTrue(matches.contains("Lisinopril"));
        assertTrue(matches.contains("Norvasc"));
        assertTrue(matches.contains("Azithromycin"));
    }

    @Test
    public void matchMultipleSimilarTreatments() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final List<CuratedTreatment> curatedTreatments = curator.matchMultiple("amlodipine besylate, amlodipine acetate");
        assertEquals(2, curatedTreatments.size());
        final List<String> matches = curatedTreatments.stream().map(CuratedTreatment::name).collect(Collectors.toList());
        assertTrue(matches.contains("Norvasc"));
        assertTrue(matches.contains("Norvinopril"));
    }

    @Test
    public void doesNotMatchAmbiguousTerm() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final List<CuratedTreatment> acidCuratedTreatments = curator.search("acid");
        assertEquals(0, acidCuratedTreatments.size());
        final List<CuratedTreatment> amlodipineCuratedTreatments = curator.search("amlodipine");
        assertEquals(0, amlodipineCuratedTreatments.size());
    }

    @Test
    public void doesNotMatchAmbiguousMultiTerm() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final List<CuratedTreatment> curatedTreatments = curator.search("pain therapy");
        assertEquals(0, curatedTreatments.size());
    }

    @Test
    public void doesNotMatchNonExistentComposedTerm() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final List<CuratedTreatment> curatedTreatments = curator.search("amlodipine phosphate");
        assertEquals(0, curatedTreatments.size());
    }

    @Test
    public void matchesTermWithAlias() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final List<CuratedTreatment> curatedTreatments = curator.search("Zocor (simvastatin)");
        assertEquals(2, curatedTreatments.size());
        assertEquals("Zocor", curatedTreatments.get(0).name());
        assertEquals("Zocor", curatedTreatments.get(1).name());
    }

    @Test
    public void matchesTermWithNumbers() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final List<CuratedTreatment> curatedTreatments = curator.search("TNT 101");
        assertEquals(1, curatedTreatments.size());
        assertEquals("TNT-101", curatedTreatments.get(0).name());
    }

    @Test
    public void doesNotMatchSingleOccurrenceOfAmbiguousTerm() {
        final TreatmentCurator curator = TestCuratorFactory.treatmentCurator();
        final List<CuratedTreatment> curatedTreatments = curator.search("acetate");
        assertEquals(0, curatedTreatments.size());
    }
}
