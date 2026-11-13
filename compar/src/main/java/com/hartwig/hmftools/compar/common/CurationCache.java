package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.MismatchFile.FLD_CURATED_BY;
import static com.hartwig.hmftools.compar.MismatchFile.FLD_CURATED_DATE;
import static com.hartwig.hmftools.compar.MismatchFile.FLD_CURATION_COMMENT;
import static com.hartwig.hmftools.compar.MismatchFile.FLD_CURATION_TYPE;
import static com.hartwig.hmftools.compar.MismatchFile.FLD_DIFFERENCES;
import static com.hartwig.hmftools.compar.MismatchFile.FLD_MISMATCH_TYPE;
import static com.hartwig.hmftools.compar.MismatchFile.FLD_NEW_SOURCE;
import static com.hartwig.hmftools.compar.MismatchFile.FLD_OLD_SOURCE;
import static com.hartwig.hmftools.compar.MismatchFile.FLD_SEVERITY;
import static com.hartwig.hmftools.compar.MismatchFile.FLD_VERDICT;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_CATEGORY;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_ITEM_KEY;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_SAMPLE_ID;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.Instant;
import java.time.LocalDate;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class CurationCache
{
    private final Map<String, Map<CategoryType, Map<String, Entry>>> mLatestEntries;
    private final Path mFilePath;
    private int mLoadCount;

    private static final List<String> REQUIRED_FIELDS = List.of(
            FLD_SAMPLE_ID, FLD_CATEGORY, FLD_MISMATCH_TYPE, FLD_ITEM_KEY, FLD_DIFFERENCES, FLD_OLD_SOURCE, FLD_NEW_SOURCE,
            FLD_CURATION_TYPE, FLD_SEVERITY, FLD_VERDICT, FLD_CURATION_COMMENT, FLD_CURATED_BY, FLD_CURATED_DATE);

    public CurationCache(final String filename)
    {
        mLatestEntries = Maps.newHashMap();
        mFilePath = Paths.get(filename);
        mLoadCount = 0;
        loadFile(filename);
    }

    public static boolean isValidSourceLabel(final String label)
    {
        return label != null && !label.isBlank() && label.chars().noneMatch(x -> x == '\t' || x == '\r' || x == '\n');
    }

    public boolean isCacheFile(final String filename)
    {
        Path path = Paths.get(filename);

        if(!Files.exists(path))
            return false;

        try
        {
            return Files.isSameFile(mFilePath, path);
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to compare file({}) with curation cache file: {}", filename, e.toString());
            return true;
        }
    }

    public record Entry(
            MismatchType type, String differences, CurationType curationType, Severity severity, Verdict verdict, String comment,
            String curatedBy, String curatedDate, String oldSource, String newSource, int loadIndex)
    {
    }

    public record Match(Entry previous, CurationStatus status)
    {
    }

    public Match match(final String sampleId, final Mismatch mismatch)
    {
        Map<String, Entry> keyEntries = mLatestEntries.getOrDefault(sampleId, Collections.emptyMap())
                .getOrDefault(mismatch.nonNullItem().category(), Collections.emptyMap());

        Entry oldEntry = mismatch.OldItem != null ? keyEntries.get(mismatch.OldItem.key()) : null;
        Entry newEntry = mismatch.NewItem != null ? keyEntries.get(mismatch.NewItem.key()) : null;

        Entry latest = oldEntry;

        if(newEntry != null && (oldEntry == null || newEntry.loadIndex() > oldEntry.loadIndex()))
            latest = newEntry;

        if(latest == null)
            return new Match(null, null);

        String differences = String.join(ITEM_DELIM, mismatch.DiffValues);
        boolean acknowledged = latest.type() == mismatch.Type && latest.differences().equals(differences);

        return new Match(latest, acknowledged ? CurationStatus.ACKNOWLEDGED : CurationStatus.CHANGED_FROM_CACHE);
    }

    private void loadFile(final String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));

            if(lines.isEmpty())
                throw new IOException("curation cache file is empty");

            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);

            for(String field : REQUIRED_FIELDS)
            {
                if(!fieldsIndexMap.containsKey(field))
                {
                    CMP_LOGGER.error("curation cache file({}) missing column({})", filename, field);
                    System.exit(1);
                }
            }

            int sampleIndex = fieldsIndexMap.get(FLD_SAMPLE_ID);
            int categoryIndex = fieldsIndexMap.get(FLD_CATEGORY);
            int mismatchTypeIndex = fieldsIndexMap.get(FLD_MISMATCH_TYPE);
            int keyIndex = fieldsIndexMap.get(FLD_ITEM_KEY);
            int differencesIndex = fieldsIndexMap.get(FLD_DIFFERENCES);
            int oldSourceIndex = fieldsIndexMap.get(FLD_OLD_SOURCE);
            int newSourceIndex = fieldsIndexMap.get(FLD_NEW_SOURCE);
            int curationTypeIndex = fieldsIndexMap.get(FLD_CURATION_TYPE);
            int severityIndex = fieldsIndexMap.get(FLD_SEVERITY);
            int verdictIndex = fieldsIndexMap.get(FLD_VERDICT);
            int curationCommentIndex = fieldsIndexMap.get(FLD_CURATION_COMMENT);
            int curatedByIndex = fieldsIndexMap.get(FLD_CURATED_BY);
            int curatedDateIndex = fieldsIndexMap.get(FLD_CURATED_DATE);

            for(int lineIndex = 1; lineIndex < lines.size(); ++lineIndex)
            {
                String[] values = lines.get(lineIndex).split(TSV_DELIM, -1);

                try
                {
                    MismatchType mismatchType = MismatchType.valueOf(values[mismatchTypeIndex]);
                    String sampleId = values[sampleIndex];
                    CategoryType category = CategoryType.valueOf(values[categoryIndex]);
                    String oldSource = values[oldSourceIndex];
                    String newSource = values[newSourceIndex];

                    if(sampleId.isEmpty())
                        throw new IllegalArgumentException("missing SampleId");

                    if(!isValidSourceLabel(oldSource) || !isValidSourceLabel(newSource))
                        throw new IllegalArgumentException("source labels must be nonblank single-line values without tabs");

                    CurationType curationType = CurationType.valueOf(values[curationTypeIndex]);
                    Severity severity = Severity.valueOf(values[severityIndex]);
                    String verdictStr = values[verdictIndex];
                    Verdict verdict = verdictStr.isEmpty() ? Verdict.NONE : Verdict.valueOf(verdictStr);
                    String curatedBy = values[curatedByIndex];
                    String curatedDate = values[curatedDateIndex];

                    if(curatedBy.isBlank() || curatedDate.isBlank())
                        throw new IllegalArgumentException("missing CuratedBy or CuratedDate");

                    if(curatedDate.contains("T"))
                        Instant.parse(curatedDate);
                    else
                        LocalDate.parse(curatedDate);

                    Entry entry = new Entry(
                            mismatchType, values[differencesIndex], curationType, severity, verdict, values[curationCommentIndex],
                            curatedBy, curatedDate, oldSource, newSource, mLoadCount++);

                    mLatestEntries.computeIfAbsent(sampleId, x -> Maps.newHashMap())
                            .computeIfAbsent(category, x -> Maps.newHashMap())
                            .put(values[keyIndex], entry);
                }
                catch(Exception e)
                {
                    CMP_LOGGER.error("invalid curation cache file({}) row({}): {}", filename, lineIndex + 1, e.getMessage());
                    System.exit(1);
                }
            }

            CMP_LOGGER.info("loaded {} curation cache entries from file({})", mLoadCount, filename);
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to load curation cache file({}): {}", filename, e.toString());
            System.exit(1);
        }
    }
}
