package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.compar.common.CurationCache.isValidSourceLabel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.compar.MismatchFile;
import com.hartwig.hmftools.compar.common.CurationCache.Entry;
import com.hartwig.hmftools.compar.isofox.RnaFusionData;
import com.hartwig.hmftools.compar.isofox.TestRnaFusionDataBuilder;
import com.hartwig.hmftools.compar.purple.PurityData;
import com.hartwig.hmftools.compar.purple.TestPurityDataBuilder;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class CurationCacheTest
{
    @Rule
    public TemporaryFolder temp = new TemporaryFolder();

    private static final String SAMPLE_ID = "TEST";
    private static final String OTHER_SAMPLE_ID = "TEST2";
    private static final String FUSION = CategoryType.RNA_FUSION.toString();

    private static final String LABEL_A = "pipeline a production";
    private static final String LABEL_B = "pipeline b production";
    private static final String LABEL_C = "pipeline c production";

    private static final String REVIEWER = "reviewer a";
    private static final String CURATED_DATE = "2026-10-20";
    private static final String CURATION_HEADER = "CurationType\tVerdict\tCurationComment\tCuratedBy\tCuratedDate";

    private static final RnaFusionData FUSION_ITEM = TestRnaFusionDataBuilder.BUILDER.create();
    private static final RnaFusionData OTHER_FUSION_ITEM = TestRnaFusionDataBuilder.BUILDER.createWithAlternateDefaults();

    @Test
    public void sourceLabelsAcceptFreeText()
    {
        assertTrue(isValidSourceLabel("pipeline 5.34 production"));
        assertTrue(isValidSourceLabel("kind=pipeline;version=5.34;reference=GRCh38;config=production"));
        assertFalse(isValidSourceLabel(null));

        for(String label : List.of("", " ", "old\tnew", "old\nnew", "old\rnew"))
        {
            assertFalse(label, isValidSourceLabel(label));
        }
    }

    @Test
    public void entriesArePerSampleInFileOrderAndPairedKeysMerge() throws IOException
    {
        File cacheFile = writeCacheFile(
                "cache.tsv", List.of(
                        List.of(SAMPLE_ID, FUSION, "OLD_ONLY", FUSION_ITEM.key(), "", LABEL_A, LABEL_B, "EXPECTED", "known"),
                        List.of(SAMPLE_ID, FUSION, "FULL_MATCH", OTHER_FUSION_ITEM.key(), "", LABEL_A, LABEL_B, "NONE", ""),
                        List.of(OTHER_SAMPLE_ID, FUSION, "VALUE", FUSION_ITEM.key(), "SplitFrags(10/2)", LABEL_A, LABEL_B, "", ""),
                        List.of(SAMPLE_ID, FUSION, "NEW_ONLY", OTHER_FUSION_ITEM.key(), "", LABEL_B, LABEL_C, "NONE", ""),
                        List.of(SAMPLE_ID, FUSION, "INVALID_NEW", "", "", LABEL_B, LABEL_C, "NONE", "")));

        CurationCache cache = new CurationCache(cacheFile.getAbsolutePath());

        List<Entry> history = cache.findMismatches(
                SAMPLE_ID, new Mismatch(FUSION_ITEM, OTHER_FUSION_ITEM, MismatchType.VALUE, Collections.emptyList()));

        assertEquals(2, history.size());
        assertEquals(CurationType.EXPECTED, history.get(0).curation());
        assertEquals("known", history.get(0).comment());
        assertEquals(LABEL_B, history.get(1).oldSource());

        List<Entry> otherSampleHistory = cache.findMismatches(
                OTHER_SAMPLE_ID, new Mismatch(FUSION_ITEM, null, MismatchType.OLD_ONLY, Collections.emptyList()));

        assertEquals(1, otherSampleHistory.size());
        assertEquals(CurationType.NONE, otherSampleHistory.get(0).curation());
    }

    @Test
    public void cacheFileIsMatchedByFileIdentity() throws IOException
    {
        File cacheFile = writeCacheFile(
                "cache.tsv", List.of(
                        List.of(SAMPLE_ID, FUSION, "OLD_ONLY", FUSION_ITEM.key(), "", LABEL_A, LABEL_B, "NONE", "")));

        Path link = Files.createSymbolicLink(temp.getRoot().toPath().resolve("link.tsv"), cacheFile.toPath());

        CurationCache cache = new CurationCache(cacheFile.getAbsolutePath());

        assertTrue(cache.isCacheFile(cacheFile.getParent() + "/./" + cacheFile.getName()));
        assertTrue(cache.isCacheFile(link.toString()));
        assertFalse(cache.isCacheFile(temp.newFile("other.tsv").getAbsolutePath()));
        assertFalse(cache.isCacheFile(temp.getRoot() + "/missing.tsv"));
    }

    @Test
    public void differencesAreAcknowledgedChangedOrNewPerField() throws IOException
    {
        String fusionsDifference = "fusions(/coords(1-2 frags=1);coords(3-4 frags=1))";
        File cacheFile = writeCacheFile(
                "cache.tsv", List.of(
                        List.of(
                                SAMPLE_ID, FUSION, "VALUE", FUSION_ITEM.key(), "SplitFrags(10/2);Filter(A;B/PASS);" + fusionsDifference,
                                LABEL_A, LABEL_B, "NONE", "present in cache"),
                        List.of(SAMPLE_ID, FUSION, "VALUE", FUSION_ITEM.key(), "SplitFrags(20/2)", LABEL_B, LABEL_C, "ISSUE", "later")));
        CurationCache cache = new CurationCache(cacheFile.getAbsolutePath());
        Mismatch mismatch = new Mismatch(
                FUSION_ITEM, FUSION_ITEM, MismatchType.VALUE,
                List.of("SplitFrags(10/2)", "Filter(A;B/NEW)", "KnownType(A/B)", fusionsDifference));

        List<CurationCache.Match> matches = cache.match(SAMPLE_ID, mismatch);

        assertEquals(
                Lists.newArrayList(CurationStatus.ACKNOWLEDGED, CurationStatus.CHANGED_FROM_CACHE, null, CurationStatus.ACKNOWLEDGED),
                matches.stream().map(CurationCache.Match::status).toList());
        assertEquals(
                List.of(MismatchType.VALUE, MismatchType.VALUE, MismatchType.VALUE, MismatchType.VALUE),
                matches.stream().map(x -> x.mismatch().Type).toList());
        assertEquals(List.of("SplitFrags(10/2)"), matches.get(0).mismatch().DiffValues);
        assertEquals("SplitFrags(10/2)", matches.get(0).previous().differences());
        assertEquals(LABEL_A, matches.get(0).previous().oldSource());
        assertEquals("Filter(A;B/PASS)", matches.get(1).previous().differences());
        assertEquals(fusionsDifference, matches.get(3).previous().differences());
    }

    @Test
    public void sampleLevelCacheDoesNotAcknowledgeOtherFieldsOrSamples() throws IOException
    {
        File cacheFile = writeCacheFile(
                "purity.tsv", List.of(
                        List.of(SAMPLE_ID, "PURITY", "VALUE", "", "Purity(0.6/0.8)", LABEL_A, LABEL_B, "NONE", "")));
        CurationCache cache = new CurationCache(cacheFile.getAbsolutePath());
        PurityData purity = TestPurityDataBuilder.BUILDER.create();
        Mismatch mismatch = new Mismatch(purity, purity, MismatchType.VALUE, List.of("Purity(0.6/0.8)", "Ploidy(2/3)"));

        List<CurationCache.Match> matches = cache.match(SAMPLE_ID, mismatch);

        assertEquals(CurationStatus.ACKNOWLEDGED, matches.get(0).status());
        assertEquals(null, matches.get(1).status());
        assertEquals(null, cache.match(OTHER_SAMPLE_ID, mismatch).get(0).status());
    }

    @Test
    public void presenceDirectionAndInvalidOrFullMatchesArePreserved() throws IOException
    {
        File cacheFile = writeCacheFile(
                "cache.tsv", List.of(
                        List.of(SAMPLE_ID, FUSION, "OLD_ONLY", FUSION_ITEM.key(), "", LABEL_A, LABEL_B, "NONE", "")));
        CurationCache cache = new CurationCache(cacheFile.getAbsolutePath());

        Mismatch oldOnly = new Mismatch(FUSION_ITEM, null, MismatchType.OLD_ONLY, Collections.emptyList());
        Mismatch newOnly = new Mismatch(null, FUSION_ITEM, MismatchType.NEW_ONLY, Collections.emptyList());
        Mismatch otherItem = new Mismatch(null, OTHER_FUSION_ITEM, MismatchType.NEW_ONLY, Collections.emptyList());
        Mismatch fullMatch = new Mismatch(FUSION_ITEM, FUSION_ITEM, MismatchType.FULL_MATCH, Collections.emptyList());
        Mismatch invalid = new Mismatch(FUSION_ITEM, null, MismatchType.INVALID_NEW, Collections.emptyList());

        assertEquals(CurationStatus.ACKNOWLEDGED, cache.match(SAMPLE_ID, oldOnly).get(0).status());
        assertEquals(CurationStatus.CHANGED_FROM_CACHE, cache.match(SAMPLE_ID, newOnly).get(0).status());
        assertEquals(null, cache.match(SAMPLE_ID, otherItem).get(0).status());
        assertEquals(MismatchType.FULL_MATCH, cache.match(SAMPLE_ID, fullMatch).get(0).mismatch().Type);
        assertEquals(MismatchType.INVALID_NEW, cache.match(SAMPLE_ID, invalid).get(0).mismatch().Type);
    }

    @Test
    public void previousPrefersLatestCuratedEntry() throws IOException
    {
        File cacheFile = writeCacheFile(
                "cache.tsv", List.of(
                        List.of(SAMPLE_ID, FUSION, "OLD_ONLY", FUSION_ITEM.key(), "", LABEL_A, LABEL_B, "EXPECTED", "known"),
                        List.of(SAMPLE_ID, FUSION, "OLD_ONLY", FUSION_ITEM.key(), "", LABEL_B, LABEL_C, "NONE", "")));
        CurationCache cache = new CurationCache(cacheFile.getAbsolutePath());

        Mismatch oldOnly = new Mismatch(FUSION_ITEM, null, MismatchType.OLD_ONLY, Collections.emptyList());
        CurationCache.Match match = cache.match(SAMPLE_ID, oldOnly).get(0);

        assertEquals(MismatchType.OLD_ONLY, match.mismatch().Type);
        assertEquals(CurationStatus.ACKNOWLEDGED, match.status());
        assertEquals(LABEL_A, match.previous().oldSource());
        assertEquals("known", match.previous().comment());
    }

    @Test
    public void verdictCountsAsCuratedAndStampsAreLoaded() throws IOException
    {
        File cacheFile = writeCacheFile(
                "cache.tsv", List.of(
                        List.of(
                                SAMPLE_ID, FUSION, "OLD_ONLY", FUSION_ITEM.key(), "", LABEL_A, LABEL_B, "NONE", "regression", "OLD_CORRECT",
                                "2026-10-20T03:15:42Z"),
                        List.of(SAMPLE_ID, FUSION, "OLD_ONLY", FUSION_ITEM.key(), "", LABEL_B, LABEL_C, "NONE", "")));
        CurationCache cache = new CurationCache(cacheFile.getAbsolutePath());

        Mismatch oldOnly = new Mismatch(FUSION_ITEM, null, MismatchType.OLD_ONLY, Collections.emptyList());
        CurationCache.Entry previous = cache.match(SAMPLE_ID, oldOnly).get(0).previous();

        assertEquals(Verdict.OLD_CORRECT, previous.verdict());
        assertEquals(LABEL_A, previous.oldSource());
        assertEquals(REVIEWER, previous.curatedBy());
        assertEquals("2026-10-20T03:15:42Z", previous.curatedDate());
    }

    @Test
    public void bothCorrectVerdictIsLoaded() throws IOException
    {
        File cacheFile = writeCacheFile(
                "cache.tsv", List.of(
                        List.of(
                                SAMPLE_ID, FUSION, "OLD_ONLY", FUSION_ITEM.key(), "", LABEL_A, LABEL_B, "NONE", "preference", "BOTH_CORRECT",
                                "2026-10-23")));
        CurationCache cache = new CurationCache(cacheFile.getAbsolutePath());

        Mismatch oldOnly = new Mismatch(FUSION_ITEM, null, MismatchType.OLD_ONLY, Collections.emptyList());
        assertEquals(Verdict.BOTH_CORRECT, cache.match(SAMPLE_ID, oldOnly).get(0).previous().verdict());
    }

    private File writeCacheFile(final String filename, final List<List<String>> rows) throws IOException
    {
        List<String> lines = Lists.newArrayList(MismatchFile.header(true, false) + TSV_DELIM + CURATION_HEADER);

        for(List<String> row : rows)
        {
            lines.add(String.join(
                    TSV_DELIM, row.get(0), row.get(1), row.get(2), "", row.get(3), row.get(4), "", "", row.get(5),
                    row.get(6), row.get(7), row.size() > 9 ? row.get(9) : "", row.get(8), REVIEWER,
                    row.size() > 10 ? row.get(10) : CURATED_DATE));
        }

        File file = temp.newFile(filename);
        Files.write(file.toPath(), lines);
        return file;
    }
}
