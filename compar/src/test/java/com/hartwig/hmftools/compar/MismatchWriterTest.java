package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;
import com.hartwig.hmftools.compar.common.WriteType;
import com.hartwig.hmftools.compar.isofox.RnaFusionComparer;
import com.hartwig.hmftools.compar.isofox.RnaFusionData;
import com.hartwig.hmftools.compar.isofox.TestRnaFusionDataBuilder;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class MismatchWriterTest
{
    @Rule
    public TemporaryFolder temp = new TemporaryFolder();

    private static final String SAMPLE_ID = "TEST";

    private static final String LABEL_A = "pipeline a production";
    private static final String LABEL_B = "kind=pipeline;version=unfinished";
    private static final String LABEL_C = "kind=pipeline;version=c;reference=GRCh38;config=none";

    private static final String REVIEWER = "reviewer a";
    private static final String CURATED_DATE = "2026-10-20";
    private static final String CURATION_HEADER = "CurationType\tVerdict\tCurationComment\tCuratedBy\tCuratedDate";

    private static final RnaFusionData FUSION_ITEM = TestRnaFusionDataBuilder.BUILDER.create();
    private static final RnaFusionData OTHER_FUSION_ITEM = TestRnaFusionDataBuilder.BUILDER.createWithAlternateDefaults();

    private static final Mismatch FUSION_OLD_ONLY = new Mismatch(FUSION_ITEM, null, MismatchType.OLD_ONLY, Collections.emptyList());
    private static final Mismatch OTHER_FUSION_OLD_ONLY = new Mismatch(
            OTHER_FUSION_ITEM, null, MismatchType.OLD_ONLY, Collections.emptyList());

    private static final List<String> PREV_HEADER = List.of(
            "PrevComparison", "PrevMismatchType", "PrevDifferences", "PrevCurationType", "PrevVerdict", "PrevCurationComment",
            "PrevCuratedBy", "PrevCuratedDate");

    private static final List<String> NO_HISTORY = Collections.nCopies(PREV_HEADER.size(), "");

    @Test
    public void defaultOutputHasNoCurationColumns() throws IOException
    {
        List<String> lines = writeOutput(createConfig(temp.newFolder(), null, null, null), List.of(FUSION_OLD_ONLY));

        assertEquals("SampleId\tCategory\tMismatchType\tKey\tDifferences\tOldValues\tNewValues", lines.get(0));
        assertEquals(2, lines.size());
    }

    @Test
    public void sourceLabelsWriteArchiveColumns() throws IOException
    {
        List<String> lines = writeOutput(createConfig(temp.newFolder(), LABEL_A, LABEL_B, null), List.of(FUSION_OLD_ONLY));

        List<String> header = Arrays.asList(lines.get(0).split(TSV_DELIM, -1));
        assertEquals("SampleId", header.get(0));
        assertEquals("CurationStatus", header.get(3));
        assertEquals(List.of("OldSource", "NewSource"), header.subList(8, header.size()));

        List<String> row = Arrays.asList(lines.get(1).split(TSV_DELIM, -1));
        assertEquals(header.size(), row.size());
        assertEquals(SAMPLE_ID, row.get(0));
        assertEquals("", row.get(3));
        assertEquals(List.of(LABEL_A, LABEL_B), row.subList(8, row.size()));
    }

    @Test
    public void historyPrefersLatestCuratedRowForTheChangedField() throws IOException
    {
        File cacheFile = writeCacheFile("curation.tsv", List.of(
                List.of("VALUE", FUSION_ITEM.key(), "SplitFrags(10/3)", LABEL_A, LABEL_B, "EXPECTED", "known", "NEW_CORRECT"),
                List.of("VALUE", FUSION_ITEM.key(), "SplitFrags(10/2)", LABEL_B, LABEL_C, "NONE", ""),
                List.of("VALUE", FUSION_ITEM.key(), "Filter(PASS/FRAGMENT_COUNT)", LABEL_B, LABEL_C, "NONE", "")));

        List<String> lines = writeOutput(
                createConfig(temp.newFolder(), null, null, cacheFile.getAbsolutePath()),
                List.of(new Mismatch(FUSION_ITEM, FUSION_ITEM, MismatchType.VALUE, List.of("SplitFrags(100/10)")), OTHER_FUSION_OLD_ONLY));

        List<String> header = Arrays.asList(lines.get(0).split(TSV_DELIM, -1));
        assertEquals(PREV_HEADER, header.subList(header.size() - PREV_HEADER.size(), header.size()));

        assertEquals(
                List.of(LABEL_A + ">" + LABEL_B, "VALUE", "SplitFrags(10/3)", "EXPECTED", "NEW_CORRECT", "known", REVIEWER, CURATED_DATE),
                history(lines.get(1)));

        assertEquals(NO_HISTORY, history(lines.get(2)));
    }

    @Test
    public void reviewedOutputIsReusedAsCache() throws IOException
    {
        List<String> firstRunLines = writeOutput(createConfig(temp.newFolder(), LABEL_A, LABEL_B, null), List.of(FUSION_OLD_ONLY));

        File cacheFile = temp.newFile("001_reviewed.tsv");
        String reviewedRow = review(firstRunLines.get(0), firstRunLines.get(1), "EXPECTED", "known");
        Files.write(cacheFile.toPath(), List.of(firstRunLines.get(0) + TSV_DELIM + CURATION_HEADER, reviewedRow));

        List<String> lines = writeOutput(
                createConfig(temp.newFolder(), LABEL_B, LABEL_C, cacheFile.getAbsolutePath()), List.of(FUSION_OLD_ONLY));

        assertEquals(List.of(LABEL_A + ">" + LABEL_B, "OLD_ONLY", "", "EXPECTED", "NONE", "known", REVIEWER, CURATED_DATE),
                history(lines.get(1)));
    }

    @Test
    public void rowsFromCacheRunsAppendToFirstArchive() throws IOException
    {
        Mismatch firstMismatch = new Mismatch(FUSION_ITEM, FUSION_ITEM, MismatchType.VALUE, List.of("SplitFrags(10/2)"));
        Mismatch secondMismatch = new Mismatch(FUSION_ITEM, FUSION_ITEM, MismatchType.VALUE, List.of("SplitFrags(10/3)"));

        List<String> archived = writeOutput(createConfig(temp.newFolder(), LABEL_A, LABEL_B, null), List.of(firstMismatch));
        File cacheFile = temp.newFile("pipeline_curation.tsv");
        String reviewedRow = review(archived.get(0), archived.get(1), "NONE", "");
        Files.write(cacheFile.toPath(), List.of(archived.get(0) + TSV_DELIM + CURATION_HEADER, reviewedRow));

        List<String> changed = writeOutput(
                createConfig(temp.newFolder(), LABEL_B, LABEL_C, cacheFile.getAbsolutePath()), List.of(secondMismatch));

        Files.write(cacheFile.toPath(), List.of(review(changed.get(0), changed.get(1), "NONE", "")), StandardOpenOption.APPEND);

        List<String> lines = writeOutput(
                createConfig(temp.newFolder(), LABEL_C, LABEL_A, cacheFile.getAbsolutePath()), List.of(secondMismatch));

        List<String> header = Arrays.asList(lines.get(0).split(TSV_DELIM, -1));
        List<String> changedRow = Arrays.asList(changed.get(1).split(TSV_DELIM, -1));
        List<String> repeatedRow = Arrays.asList(lines.get(1).split(TSV_DELIM, -1));
        assertEquals("VALUE", changedRow.get(header.indexOf("MismatchType")));
        assertEquals("CHANGED_FROM_CACHE", changedRow.get(header.indexOf("CurationStatus")));
        assertEquals("VALUE", repeatedRow.get(header.indexOf("MismatchType")));
        assertEquals("ACKNOWLEDGED", repeatedRow.get(header.indexOf("CurationStatus")));
        assertEquals(List.of(LABEL_B + ">" + LABEL_C, "VALUE", "SplitFrags(10/3)", "NONE", "NONE", "", REVIEWER, CURATED_DATE),
                history(lines.get(1)));
    }

    @Test
    public void mixedDifferencesWriteSeparateRowsInBothOutputFormats() throws IOException
    {
        File cacheFile = writeCacheFile("mixed.tsv", List.of(
                List.of("VALUE", FUSION_ITEM.key(), "SplitFrags(10/2);Filter(PASS/A)", LABEL_A, LABEL_B, "NONE", "reviewed")));
        ComparConfig config = createConfig(temp.newFolder(), LABEL_B, LABEL_C, cacheFile.getAbsolutePath());
        config.WriteTypes.add(WriteType.TYPE_SPECIFIC);
        Mismatch mismatch = new Mismatch(FUSION_ITEM, FUSION_ITEM, MismatchType.VALUE,
                List.of("SplitFrags(10/2)", "Filter(PASS/B)", "KnownType(A/B)"));

        List<String> lines = writeOutput(config, List.of(mismatch));
        List<String> detailed = Files.readAllLines(Paths.get(config.formOutputFilePrefix() + ".rna_fusion.tsv"));

        List<String> header = Arrays.asList(lines.get(0).split(TSV_DELIM, -1));
        List<String> detailedHeader = Arrays.asList(detailed.get(0).split(TSV_DELIM, -1));
        assertEquals(4, lines.size());
        assertEquals(4, detailed.size());
        assertEquals(List.of("VALUE", "VALUE", "VALUE"),
                lines.subList(1, 4).stream().map(x -> x.split(TSV_DELIM, -1)[header.indexOf("MismatchType")]).toList());
        assertEquals(List.of("ACKNOWLEDGED", "CHANGED_FROM_CACHE", ""),
                lines.subList(1, 4).stream().map(x -> x.split(TSV_DELIM, -1)[header.indexOf("CurationStatus")]).toList());
        assertEquals(List.of("ACKNOWLEDGED", "CHANGED_FROM_CACHE", ""),
                detailed.subList(1, 4).stream().map(x -> x.split(TSV_DELIM, -1)[detailedHeader.indexOf("CurationStatus")]).toList());
        assertEquals(List.of("SplitFrags(10/2)", "Filter(PASS/B)", "KnownType(A/B)"),
                lines.subList(1, 4).stream().map(x -> x.split(TSV_DELIM, -1)[header.indexOf("Differences")]).toList());
        assertEquals(List.of(LABEL_A + ">" + LABEL_B, "VALUE", "SplitFrags(10/2)", "NONE", "NONE", "reviewed", REVIEWER, CURATED_DATE),
                history(lines.get(1)));
        assertEquals(List.of(LABEL_A + ">" + LABEL_B, "VALUE", "Filter(PASS/A)", "NONE", "NONE", "reviewed", REVIEWER, CURATED_DATE),
                history(lines.get(2)));
        assertEquals(NO_HISTORY, history(lines.get(3)));
    }

    @Test
    public void withoutCacheDifferencesStayOnOneRow() throws IOException
    {
        Mismatch mismatch = new Mismatch(FUSION_ITEM, FUSION_ITEM, MismatchType.VALUE,
                List.of("SplitFrags(10/2)", "Filter(PASS/B)"));
        List<String> lines = writeOutput(createConfig(temp.newFolder(), null, null, null), List.of(mismatch));

        List<String> header = Arrays.asList(lines.get(0).split(TSV_DELIM, -1));
        List<String> row = Arrays.asList(lines.get(1).split(TSV_DELIM, -1));
        assertEquals(2, lines.size());
        assertEquals(SAMPLE_ID, row.get(header.indexOf("SampleId")));
        assertEquals("VALUE", row.get(header.indexOf("MismatchType")));
        assertEquals("SplitFrags(10/2);Filter(PASS/B)", row.get(header.indexOf("Differences")));
        assertFalse(header.contains("CurationStatus"));
    }

    @Test
    public void outputMatchingCacheFileIsRejected() throws IOException
    {
        String filePrefix = "compar_cohort";

        for(String outputFile : List.of(filePrefix + ".tsv", filePrefix + ".rna_fusion.tsv", FieldConfigFile.generateFileName(filePrefix)))
        {
            File cacheFile = writeCacheFile(outputFile, List.of(
                    List.of("OLD_ONLY", FUSION_ITEM.key(), "", LABEL_A, LABEL_B, "EXPECTED", "reviewed")));

            List<String> reviewedLines = Files.readAllLines(cacheFile.toPath());

            ComparConfig config = createConfig(temp.getRoot(), null, null, temp.getRoot() + "/./" + cacheFile.getName());
            config.WriteTypes.add(WriteType.TYPE_SPECIFIC);
            config.WriteTypes.add(WriteType.FIELD_CONFIG);
            MismatchWriter writer = new MismatchWriter(config);

            assertFalse(outputFile, writer.initialiseOutputFiles(new FieldCheckCache()));
            writer.close();

            assertEquals(reviewedLines, Files.readAllLines(cacheFile.toPath()));
        }
    }

    private static List<String> history(final String line)
    {
        List<String> row = Arrays.asList(line.split(TSV_DELIM, -1));
        return row.subList(row.size() - PREV_HEADER.size(), row.size());
    }

    private static String review(final String header, final String line, final String curationType, final String comment)
    {
        int newSourceIndex = Arrays.asList(header.split(TSV_DELIM, -1)).indexOf("NewSource");
        List<String> values = Lists.newArrayList(Arrays.asList(line.split(TSV_DELIM, -1)).subList(0, newSourceIndex + 1));
        values.addAll(List.of(curationType, "NONE", comment, REVIEWER, CURATED_DATE));
        return String.join(TSV_DELIM, values);
    }

    private static ComparConfig createConfig(
            final File outputDir, final String oldLabel, final String newLabel, final String cacheFilePath)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        ComparConfig.addConfig(configBuilder);
        configBuilder.setValue("sample", SAMPLE_ID);
        configBuilder.setValue("sample_dir_old", outputDir.getAbsolutePath());
        configBuilder.setValue("sample_dir_new", outputDir.getAbsolutePath());
        configBuilder.setValue("output_dir", outputDir.getAbsolutePath());
        configBuilder.setValue(ComparConfig.CATEGORIES, CategoryType.RNA_FUSION.toString());
        configBuilder.setValue(ComparConfig.WRITE_TYPES, "GENERIC");

        if(oldLabel != null)
        {
            configBuilder.setValue(ComparConfig.SOURCE_LABEL + "_old", oldLabel);
            configBuilder.setValue(ComparConfig.SOURCE_LABEL + "_new", newLabel);
        }

        if(cacheFilePath != null)
            configBuilder.setValue(ComparConfig.CURATION_CACHE, cacheFilePath);

        return new ComparConfig(configBuilder);
    }

    private List<String> writeOutput(final ComparConfig config, final List<Mismatch> mismatches) throws IOException
    {
        MismatchWriter writer = new MismatchWriter(config);

        try
        {
            writer.initialiseOutputFiles(new FieldCheckCache());
            writer.writeSampleMismatches(SAMPLE_ID, new RnaFusionComparer(config, Collections.emptyMap()), Lists.newArrayList(mismatches));
        }
        finally
        {
            writer.close();
        }

        return Files.readAllLines(Paths.get(config.formOutputFilePrefix() + ".tsv"));
    }

    private File writeCacheFile(final String filename, final List<List<String>> rows) throws IOException
    {
        List<String> lines = Lists.newArrayList(MismatchFile.header(true, false) + TSV_DELIM + CURATION_HEADER);

        for(List<String> row : rows)
        {
            lines.add(String.join(TSV_DELIM,
                    SAMPLE_ID, CategoryType.RNA_FUSION.toString(), row.get(0), "", row.get(1), row.get(2), "", "", row.get(3),
                    row.get(4), row.get(5), row.size() > 7 ? row.get(7) : "", row.get(6), REVIEWER, CURATED_DATE));
        }

        File file = temp.newFile(filename);
        Files.write(file.toPath(), lines);
        return file;
    }
}
