package com.hartwig.hmftools.patientdb.clinical.curators;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedDrug;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class TreatmentCuratorTest {

    @Test
    public void canConvertStringFieldToSynonyms() {
        assertTrue(TreatmentCurator.toSynonyms(Strings.EMPTY).isEmpty());
        assertEquals("hi", TreatmentCurator.toSynonyms("hi").get(0));

        List<String> complex = TreatmentCurator.toSynonyms("hi, \"hi, there\", second, third, \"fourth, fifth\"");
        assertEquals("hi", complex.get(0));
        assertEquals("hi, there", complex.get(1));
        assertEquals("second", complex.get(2));
        assertEquals("third", complex.get(3));
        assertEquals("fourth, fifth", complex.get(4));
    }

    @Test
    public void matchesExactSingleWord() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        Optional<CuratedDrug> curatedTreatment = curator.matchSingle("Zocor");
        assertTrue(curatedTreatment.isPresent());
        assertEquals("Zocor", curatedTreatment.get().name());
    }

    @Test
    public void removesSearchTermsFromUnused() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        assertEquals(13, curator.unusedSearchTerms().size());

        // Match with a canonical
        assertEquals(1, curator.search("Zocor").size());
        assertEquals(13, curator.unusedSearchTerms().size());
        assertEquals(1, curator.search("Zocor").size());

        // Match with "other" name
        assertEquals(1, curator.search("amlodipine-besylate").size());
        assertEquals(12, curator.unusedSearchTerms().size());
        assertEquals(1, curator.search("amlodipine-besylate").size());

        // Do an imperfect match but still remove the search term.
        assertEquals(1, curator.search("Avastine").size());
        assertEquals(11, curator.unusedSearchTerms().size());
        assertEquals(1, curator.search("Avastine").size());

        // Bogus example!!
        assertEquals(0, curator.search("This does not match at all!").size());
        assertEquals(11, curator.unusedSearchTerms().size());
    }

    @Test
    public void matchesIgnoringCaseSingleWord() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        Optional<CuratedDrug> curatedTreatment = curator.matchSingle("zocor");
        assertTrue(curatedTreatment.isPresent());
        assertEquals("Zocor", curatedTreatment.get().name());
    }

    @Test
    public void matchesSingleWordWithTypo() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        Optional<CuratedDrug> curatedTreatment = curator.matchSingle("lisinoprill");
        assertTrue(curatedTreatment.isPresent());
        assertEquals("Lisinopril", curatedTreatment.get().name());
    }

    @Test
    public void matchesExactMultiWord() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        Optional<CuratedDrug> curatedTreatment = curator.matchSingle("amlodipine besylate");
        assertTrue(curatedTreatment.isPresent());
        assertEquals("Norvasc", curatedTreatment.get().name());
    }

    @Test
    public void matchesIgnoringCaseMultiWord() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        Optional<CuratedDrug> curatedTreatment = curator.matchSingle("Amlodipine Besylate");
        assertTrue(curatedTreatment.isPresent());
        assertEquals("Norvasc", curatedTreatment.get().name());
    }

    @Test
    public void matchesMultiWordWithTypo() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        Optional<CuratedDrug> curatedTreatment = curator.matchSingle("Amlodipin besylat");
        assertTrue(curatedTreatment.isPresent());
        assertEquals("Norvasc", curatedTreatment.get().name());
    }

    @Test
    public void matchesTermWithSpecialChars() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        Optional<CuratedDrug> curatedTreatment = curator.matchSingle("z-pak");
        assertTrue(curatedTreatment.isPresent());
        assertEquals("Azithromycin", curatedTreatment.get().name());
    }

    @Test
    public void matchMultipleTreatments() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        List<CuratedDrug> curatedDrugs = curator.matchMultiple("Prinivil,Zithromax/amlodipine besylate");
        assertEquals(3, curatedDrugs.size());
        List<String> matches = curatedDrugs.stream().map(CuratedDrug::name).collect(Collectors.toList());
        assertTrue(matches.contains("Lisinopril"));
        assertTrue(matches.contains("Norvasc"));
        assertTrue(matches.contains("Azithromycin"));
    }

    @Test
    public void matchMultipleTreatmentsWithTypos() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        List<CuratedDrug> curatedDrugs = curator.matchMultiple("Prinivyl,Zithromaxx/amlodipin Besylate");
        assertEquals(3, curatedDrugs.size());
        List<String> matches = curatedDrugs.stream().map(CuratedDrug::name).collect(Collectors.toList());
        assertTrue(matches.contains("Lisinopril"));
        assertTrue(matches.contains("Norvasc"));
        assertTrue(matches.contains("Azithromycin"));
    }

    @Test
    public void matchMultipleSimilarTreatments() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        List<CuratedDrug> curatedDrugs = curator.matchMultiple("amlodipine besylate, amlodipine acetate");
        assertEquals(2, curatedDrugs.size());
        List<String> matches = curatedDrugs.stream().map(CuratedDrug::name).collect(Collectors.toList());
        assertTrue(matches.contains("Norvasc"));
        assertTrue(matches.contains("Norvinopril"));
    }

    @Test
    public void doesNotMatchAmbiguousTerm() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        List<CuratedDrug> acidCuratedDrugs = curator.search("acid");
        assertEquals(0, acidCuratedDrugs.size());
        List<CuratedDrug> amlodipineCuratedDrugs = curator.search("amlodipine");
        assertEquals(0, amlodipineCuratedDrugs.size());
    }

    @Test
    public void doesNotMatchAmbiguousMultiTerm() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        List<CuratedDrug> curatedDrugs = curator.search("pain therapy");
        assertEquals(0, curatedDrugs.size());
    }

    @Test
    public void doesNotMatchNonExistentComposedTerm() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        List<CuratedDrug> curatedDrugs = curator.search("amlodipine phosphate");
        assertEquals(0, curatedDrugs.size());
    }

    @Test
    public void matchesTermWithAlias() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        List<CuratedDrug> curatedDrugs = curator.search("Zocor (simvastatin)");
        assertEquals(2, curatedDrugs.size());
        assertEquals("Zocor", curatedDrugs.get(0).name());
        assertEquals("Zocor", curatedDrugs.get(1).name());
    }

    @Test
    public void matchesTermWithNumbers() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        List<CuratedDrug> curatedDrugs = curator.search("TNT 101");
        assertEquals(1, curatedDrugs.size());
        assertEquals("TNT-101", curatedDrugs.get(0).name());
    }

    @Test
    public void doesNotMatchSingleOccurrenceOfAmbiguousTerm() {
        TreatmentCurator curator = CuratorTestFactory.treatmentCurator();
        List<CuratedDrug> curatedDrugs = curator.search("acetate");
        assertEquals(0, curatedDrugs.size());
    }
}
