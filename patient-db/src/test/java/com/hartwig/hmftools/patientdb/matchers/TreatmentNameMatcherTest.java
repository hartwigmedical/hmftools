package com.hartwig.hmftools.patientdb.matchers;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.lucene.queryparser.classic.ParseException;
import org.junit.Test;

public class TreatmentNameMatcherTest {
    private static final Logger LOGGER = LogManager.getLogger(TreatmentNameMatcher.class);
    private static final String TREATMENT_NAMES_CSV = Resources.getResource("treatment_name_mapping.csv").getPath();
    private static final TreatmentNameMatcher MATCHER;

    static {
        try {
            MATCHER = new TreatmentNameMatcher(TREATMENT_NAMES_CSV);
        } catch (IOException e) {
            LOGGER.error(e);
            throw new RuntimeException(e);
        }
    }

    @Test
    public void matchesExactSingleWord() throws IOException, ParseException {
        assertEquals("Zocor", MATCHER.match("Zocor"));
    }

    @Test
    public void matchesIgnoringCaseSingleWord() throws IOException, ParseException {
        assertEquals("Zocor", MATCHER.match("zocor"));
    }

    @Test
    public void matchesSingleWordWithTypo() throws IOException, ParseException {
        assertEquals("Zocor", MATCHER.match("zoco"));
    }

    @Test
    public void matchesExactMultiWord() throws IOException, ParseException {
        assertEquals("Norvasc", MATCHER.match("amlodipine besylate"));
    }

    @Test
    public void matchesIgnoringCaseMultiWord() throws IOException, ParseException {
        assertEquals("Norvasc", MATCHER.match("Amlodipine Besylate"));
    }

    @Test
    public void matchesMultiWordWithTypo() throws IOException, ParseException {
        assertEquals("Norvasc", MATCHER.match("Amlodipin besylat"));
    }

    @Test
    public void matchesTermWithSpecialChars() throws IOException, ParseException {
        assertEquals("Azithromycin", MATCHER.match("z-pak"));
    }

    @Test
    public void matchMultipleTreatments() throws IOException, ParseException {
        final List<String> matches = MATCHER.matchMultiple("Prinivil,Zithromax/amlodipine besylate");
        assertEquals(3, matches.size());
        assertTrue(matches.contains("Lisinopril"));
        assertTrue(matches.contains("Norvasc"));
        assertTrue(matches.contains("Azithromycin"));
    }

    @Test
    public void matchMultipleTreatmentsWithTypos() throws IOException, ParseException {
        final List<String> matches = MATCHER.matchMultiple("Prinivyl,Zithromaxx/amlodipin Besylate");
        assertEquals(3, matches.size());
        assertTrue(matches.contains("Lisinopril"));
        assertTrue(matches.contains("Norvasc"));
        assertTrue(matches.contains("Azithromycin"));
    }

    @Test
    public void matchMultipleSimilarTreatments() throws IOException, ParseException {
        final List<String> matches = MATCHER.matchMultiple("amlodipine besylate, amlodipine acetate");
        assertEquals(2, matches.size());
        assertTrue(matches.contains("Norvasc"));
        assertTrue(matches.contains("Norvinopril"));
    }

    @Test
    public void matchPartialName() throws IOException, ParseException {
        final List<String> matches = MATCHER.search("acetylsalicylic acid");
        assertEquals(1, matches.size());
        assertTrue(matches.contains("Aspirin"));
    }

    @Test
    public void doesNotMatchAmbiguousTerm() throws IOException, ParseException {
        final List<String> acidMatches = MATCHER.search("acid");
        assertEquals(0, acidMatches.size());
        final List<String> amlodipineMatches = MATCHER.search("amlodipine");
        assertEquals(0, amlodipineMatches.size());
    }

    @Test
    public void matchesUnambiguousShortTerm() throws IOException, ParseException {
        final List<String> matches = MATCHER.search("ASA");
        assertEquals(1, matches.size());
        assertTrue(matches.contains("Aspirin"));
    }
}
