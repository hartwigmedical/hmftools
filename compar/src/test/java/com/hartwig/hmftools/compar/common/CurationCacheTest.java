package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.compar.common.CurationCache.isValidSourceLabel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
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
    private static final String CURATION_HEADER = "CurationType\tSeverity\tVerdict\tCurationComment\tCuratedBy\tCuratedDate";

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
    public void latestEntryWinsAcrossPairedKeysPerSample() throws IOException
    {
        File cacheFile = writeCacheFile(
                "cache.tsv", List.of(
                        List.of(SAMPLE_ID, FUSION, "OLD_ONLY", FUSION_ITEM.key(), "", LABEL_A, LABEL_B, "LOW", "known"),
                        List.of(SAMPLE_ID, FUSION, "FULL_MATCH", OTHER_FUSION_ITEM.key(), "", LABEL_A, LABEL_B, "NONE", ""),
                        List.of(OTHER_SAMPLE_ID, FUSION, "VALUE", FUSION_ITEM.key(), "SplitFrags(10/2)", LABEL_A, LABEL_B, "NONE", ""),
                        List.of(SAMPLE_ID, FUSION, "NEW_ONLY", OTHER_FUSION_ITEM.key(), "", LABEL_B, LABEL_C, "NONE", ""),
                        List.of(SAMPLE_ID, FUSION, "INVALID_NEW", "", "", LABEL_B, LABEL_C, "NONE", "")));

        CurationCache cache = new CurationCache(cacheFile.getAbsolutePath());

        Entry latest = cache.match(
                SAMPLE_ID, new Mismatch(FUSION_ITEM, OTHER_FUSION_ITEM, MismatchType.VALUE, Collections.emptyList())).previous();

        assertEquals(MismatchType.NEW_ONLY, latest.type());
        assertEquals(LABEL_B, latest.oldSource());

        Entry otherSampleLatest = cache.match(
                OTHER_SAMPLE_ID, new Mismatch(FUSION_ITEM, null, MismatchType.OLD_ONLY, Collections.emptyList())).previous();

        assertEquals(MismatchType.VALUE, otherSampleLatest.type());
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
    public void rowIsAcknowledgedOnlyWhenAllDifferencesMatch() throws IOException
    {
        String differences = "SplitFrags(10/2);Filter(A;B/PASS);fusions(/coords(1-2 frags=1);coords(3-4 frags=1))";
        File cacheFile = writeCacheFile(
                "cache.tsv", List.of(
                        List.of(SAMPLE_ID, FUSION, "VALUE", FUSION_ITEM.key(), differences, LABEL_A, LABEL_B, "NONE", "present in cache")));
        CurationCache cache = new CurationCache(cacheFile.getAbsolutePath());

        Mismatch sameDifferences = new Mismatch(
                FUSION_ITEM, FUSION_ITEM, MismatchType.VALUE,
                List.of("SplitFrags(10/2)", "Filter(A;B/PASS)", "fusions(/coords(1-2 frags=1);coords(3-4 frags=1))"));
        assertEquals(CurationStatus.ACKNOWLEDGED, cache.match(SAMPLE_ID, sameDifferences).status());

        Mismatch fewerDifferences = new Mismatch(FUSION_ITEM, FUSION_ITEM, MismatchType.VALUE, List.of("SplitFrags(10/2)"));
        assertEquals(CurationStatus.CHANGED_FROM_CACHE, cache.match(SAMPLE_ID, fewerDifferences).status());
    }

    @Test
    public void extraDifferenceOrOtherSampleIsNotAcknowledged() throws IOException
    {
        File cacheFile = writeCacheFile(
                "purity.tsv", List.of(
                        List.of(SAMPLE_ID, "PURITY", "VALUE", "", "Purity(0.6/0.8)", LABEL_A, LABEL_B, "NONE", "")));
        CurationCache cache = new CurationCache(cacheFile.getAbsolutePath());
        PurityData purity = TestPurityDataBuilder.BUILDER.create();
        Mismatch mismatch = new Mismatch(purity, purity, MismatchType.VALUE, List.of("Purity(0.6/0.8)", "Ploidy(2/3)"));

        assertEquals(CurationStatus.CHANGED_FROM_CACHE, cache.match(SAMPLE_ID, mismatch).status());
        assertNull(cache.match(OTHER_SAMPLE_ID, mismatch).status());
    }

    @Test
    public void presenceDirectionCountsAndEveryTypeIsMatched() throws IOException
    {
        File cacheFile = writeCacheFile(
                "cache.tsv", List.of(
                        List.of(SAMPLE_ID, FUSION, "OLD_ONLY", FUSION_ITEM.key(), "", LABEL_A, LABEL_B, "NONE", ""),
                        List.of(SAMPLE_ID, FUSION, "INVALID_NEW", "", "", LABEL_A, LABEL_B, "NONE", "")));
        CurationCache cache = new CurationCache(cacheFile.getAbsolutePath());

        Mismatch oldOnly = new Mismatch(FUSION_ITEM, null, MismatchType.OLD_ONLY, Collections.emptyList());
        Mismatch newOnly = new Mismatch(null, FUSION_ITEM, MismatchType.NEW_ONLY, Collections.emptyList());
        Mismatch otherItem = new Mismatch(null, OTHER_FUSION_ITEM, MismatchType.NEW_ONLY, Collections.emptyList());
        Mismatch fullMatch = new Mismatch(FUSION_ITEM, FUSION_ITEM, MismatchType.FULL_MATCH, Collections.emptyList());
        Mismatch invalid = new Mismatch(
                new InvalidDataItem(CategoryType.RNA_FUSION), null, MismatchType.INVALID_NEW, Collections.emptyList());

        assertEquals(CurationStatus.ACKNOWLEDGED, cache.match(SAMPLE_ID, oldOnly).status());
        assertEquals(CurationStatus.CHANGED_FROM_CACHE, cache.match(SAMPLE_ID, newOnly).status());
        assertNull(cache.match(SAMPLE_ID, otherItem).status());
        assertEquals(CurationStatus.CHANGED_FROM_CACHE, cache.match(SAMPLE_ID, fullMatch).status());
        assertEquals(MismatchType.OLD_ONLY, cache.match(SAMPLE_ID, fullMatch).previous().type());
        assertEquals(CurationStatus.ACKNOWLEDGED, cache.match(SAMPLE_ID, invalid).status());
    }

    @Test
    public void statusAndPreviousComeFromLatestEntry() throws IOException
    {
        File cacheFile = writeCacheFile(
                "cache.tsv", List.of(
                        List.of(SAMPLE_ID, FUSION, "VALUE", FUSION_ITEM.key(), "SplitFrags(10/2)", LABEL_A, LABEL_B, "LOW", "known"),
                        List.of(SAMPLE_ID, FUSION, "VALUE", FUSION_ITEM.key(), "SplitFrags(10/5)", LABEL_B, LABEL_C, "NONE", "")));
        CurationCache cache = new CurationCache(cacheFile.getAbsolutePath());

        CurationCache.Match earlierDifference = cache.match(
                SAMPLE_ID, new Mismatch(FUSION_ITEM, FUSION_ITEM, MismatchType.VALUE, List.of("SplitFrags(10/2)")));
        assertEquals(CurationStatus.CHANGED_FROM_CACHE, earlierDifference.status());
        assertEquals(LABEL_B, earlierDifference.previous().oldSource());

        CurationCache.Match latestDifference = cache.match(
                SAMPLE_ID, new Mismatch(FUSION_ITEM, FUSION_ITEM, MismatchType.VALUE, List.of("SplitFrags(10/5)")));
        assertEquals(CurationStatus.ACKNOWLEDGED, latestDifference.status());
        assertEquals(LABEL_B, latestDifference.previous().oldSource());
    }

    @Test
    public void curationColumnsAreLoaded() throws IOException
    {
        File cacheFile = writeCacheFile(
                "cache.tsv", List.of(
                        List.of(
                                SAMPLE_ID, FUSION, "OLD_ONLY", FUSION_ITEM.key(), "", LABEL_A, LABEL_B, "NONE", "regression", "OLD_CORRECT",
                                "2026-10-20T03:15:42Z", "BATCH_INFERRED")));
        CurationCache cache = new CurationCache(cacheFile.getAbsolutePath());

        Mismatch oldOnly = new Mismatch(FUSION_ITEM, null, MismatchType.OLD_ONLY, Collections.emptyList());
        Entry previous = cache.match(SAMPLE_ID, oldOnly).previous();

        assertEquals(CurationType.BATCH_INFERRED, previous.curationType());
        assertEquals(Verdict.OLD_CORRECT, previous.verdict());
        assertEquals(LABEL_A, previous.oldSource());
        assertEquals(REVIEWER, previous.curatedBy());
        assertEquals("2026-10-20T03:15:42Z", previous.curatedDate());
    }

    private File writeCacheFile(final String filename, final List<List<String>> rows) throws IOException
    {
        List<String> lines = Lists.newArrayList(MismatchFile.header(true, false) + TSV_DELIM + CURATION_HEADER);

        for(List<String> row : rows)
        {
            lines.add(String.join(
                    TSV_DELIM, row.get(0), row.get(1), row.get(2), "", row.get(3), row.get(4), "", "", row.get(5),
                    row.get(6), row.size() > 11 ? row.get(11) : "MANUAL", row.get(7), row.size() > 9 ? row.get(9) : "", row.get(8),
                    REVIEWER, row.size() > 10 ? row.get(10) : CURATED_DATE));
        }

        File file = temp.newFile(filename);
        Files.write(file.toPath(), lines);
        return file;
    }
}
