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
import static com.hartwig.hmftools.compar.MismatchFile.FLD_VERDICT;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_CATEGORY;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_ITEM_KEY;
import static com.hartwig.hmftools.compar.common.ComparConstants.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.compar.common.MismatchType.FULL_MATCH;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_BOTH;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_ERROR;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_NEW;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_OLD;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.Instant;
import java.time.LocalDate;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class CurationCache
{
    private final Map<String, Map<CategoryType, Map<String, List<Entry>>>> mSampleMismatches;
    private final Path mFilePath;
    private int mLoadCount;

    private static final Set<MismatchType> INVALID_TYPES = Set.of(INVALID_OLD, INVALID_NEW, INVALID_BOTH, INVALID_ERROR);

    private static final List<String> REQUIRED_FIELDS = List.of(
            FLD_SAMPLE_ID, FLD_CATEGORY, FLD_MISMATCH_TYPE, FLD_ITEM_KEY, FLD_DIFFERENCES, FLD_OLD_SOURCE, FLD_NEW_SOURCE,
            FLD_CURATION_TYPE, FLD_VERDICT, FLD_CURATION_COMMENT, FLD_CURATED_BY, FLD_CURATED_DATE);

    public CurationCache(final String filename)
    {
        mSampleMismatches = Maps.newHashMap();
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
            MismatchType type, String differences, CurationType curation, Verdict verdict, String comment, String curatedBy,
            String curatedDate, String oldSource, String newSource, int loadIndex)
    {
    }

    public record Match(Mismatch mismatch, Entry previous, CurationStatus status)
    {
    }

    public List<Match> match(final String sampleId, final Mismatch mismatch)
    {
        if(mismatch.Type == FULL_MATCH || INVALID_TYPES.contains(mismatch.Type))
            return List.of(new Match(mismatch, null, null));

        List<Entry> entries = findMismatches(sampleId, mismatch);
        List<Match> matches = Lists.newArrayList();
        List<String> differences = mismatch.DiffValues.isEmpty() ? List.of("") : mismatch.DiffValues;

        for(String difference : differences)
        {
            String field = differenceField(difference);
            List<Entry> fieldHistory = entries.stream()
                    .filter(x -> differenceField(x.differences()).equals(field)).toList();
            List<Entry> exactHistory = fieldHistory.stream()
                    .filter(x -> x.type() == mismatch.Type && x.differences().equals(difference)).toList();

            CurationStatus status = !exactHistory.isEmpty() ? CurationStatus.ACKNOWLEDGED :
                    !fieldHistory.isEmpty() ? CurationStatus.CHANGED_FROM_CACHE : null;
            Mismatch result = new Mismatch(
                    mismatch.OldItem, mismatch.NewItem, mismatch.Type,
                    difference.isEmpty() ? Collections.emptyList() : List.of(difference));
            List<Entry> matchingHistory = !exactHistory.isEmpty() ? exactHistory : fieldHistory;
            List<Entry> curatedHistory = matchingHistory.stream()
                    .filter(x -> x.curation() != CurationType.NONE || x.verdict() != Verdict.NONE).toList();
            List<Entry> previousHistory = !curatedHistory.isEmpty() ? curatedHistory : matchingHistory;
            Entry previous = previousHistory.isEmpty() ? null : previousHistory.get(previousHistory.size() - 1);
            matches.add(new Match(result, previous, status));
        }

        return matches;
    }

    private static String differenceField(final String difference)
    {
        int valueStart = difference.indexOf('(');
        return valueStart >= 0 ? difference.substring(0, valueStart) : difference;
    }

    private static List<String> splitDifferences(final String differences)
    {
        List<String> items = Lists.newArrayList();
        int depth = 0;
        int itemStart = 0;

        for(int index = 0; index < differences.length(); ++index)
        {
            char character = differences.charAt(index);

            if(character == '(')
                ++depth;
            else if(character == ')')
                --depth;
            else if(depth == 0 && differences.startsWith(ITEM_DELIM, index))
            {
                items.add(differences.substring(itemStart, index));
                itemStart = index + ITEM_DELIM.length();
            }
        }

        items.add(differences.substring(itemStart));
        return items;
    }

    public List<Entry> findMismatches(final String sampleId, final Mismatch mismatch)
    {
        Map<String, List<Entry>> keyMismatches = mSampleMismatches.getOrDefault(sampleId, Collections.emptyMap())
                .getOrDefault(mismatch.nonNullItem().category(), Collections.emptyMap());

        Set<String> keys = Sets.newLinkedHashSet();

        if(mismatch.OldItem != null)
            keys.add(mismatch.OldItem.key());

        if(mismatch.NewItem != null)
            keys.add(mismatch.NewItem.key());

        List<Entry> entries = Lists.newArrayList();
        keys.forEach(x -> entries.addAll(keyMismatches.getOrDefault(x, Collections.emptyList())));
        entries.sort(Comparator.comparingInt(Entry::loadIndex));

        return entries;
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

                    if(mismatchType == FULL_MATCH || INVALID_TYPES.contains(mismatchType))
                        continue;

                    String sampleId = values[sampleIndex];
                    CategoryType category = CategoryType.valueOf(values[categoryIndex]);
                    String oldSource = values[oldSourceIndex];
                    String newSource = values[newSourceIndex];

                    if(sampleId.isEmpty())
                        throw new IllegalArgumentException("missing SampleId");

                    if(!isValidSourceLabel(oldSource) || !isValidSourceLabel(newSource))
                        throw new IllegalArgumentException("source labels must be nonblank single-line values without tabs");

                    String curationTypeStr = values[curationTypeIndex];
                    CurationType curationType = curationTypeStr.isEmpty() ? CurationType.NONE : CurationType.valueOf(curationTypeStr);
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

                    List<String> differences = splitDifferences(values[differencesIndex]);

                    for(String difference : differences)
                    {
                        Entry entry = new Entry(
                                mismatchType, difference, curationType, verdict, values[curationCommentIndex], curatedBy, curatedDate,
                                oldSource, newSource, mLoadCount++);

                        mSampleMismatches.computeIfAbsent(sampleId, x -> Maps.newHashMap())
                                .computeIfAbsent(category, x -> Maps.newHashMap())
                                .computeIfAbsent(values[keyIndex], x -> Lists.newArrayList())
                                .add(entry);
                    }
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
