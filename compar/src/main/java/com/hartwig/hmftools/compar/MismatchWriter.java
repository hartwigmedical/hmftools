package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CurationCache;
import com.hartwig.hmftools.compar.common.FieldDisplayInfo;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.WriteType;
import com.hartwig.hmftools.compar.common.field.FieldValue;

public class MismatchWriter
{
    private final ComparConfig mConfig;
    private BufferedWriter mCombinedWriter;
    private final Map<CategoryType,BufferedWriter> mCategoryWriters;

    private final String mOldSourceLabel;
    private final String mNewSourceLabel;
    private final CurationCache mCurationCache;

    public MismatchWriter(final ComparConfig config)
    {
        mConfig = config;
        mCombinedWriter = null;
        mCategoryWriters = Maps.newHashMap();

        mOldSourceLabel = config.getSourceData(SourceType.OLD).Label;
        mNewSourceLabel = config.getSourceData(SourceType.NEW).Label;

        mCurationCache = config.CurationCacheFile != null ? new CurationCache(config.CurationCacheFile) : null;
    }

    public boolean initialiseOutputFiles(final FieldCheckCache fieldCheckCache)
    {
        String filePrefix = mConfig.formOutputFilePrefix();
        String combinedFile = filePrefix + TSV_EXTENSION;
        Map<String,ItemComparer> detailedFileComparers = Maps.newLinkedHashMap();

        if(mConfig.WriteTypes.contains(WriteType.TYPE_SPECIFIC))
        {
            for(ItemComparer comparer : ComparerUtils.buildComparers(mConfig, fieldCheckCache))
            {
                detailedFileComparers.put(filePrefix + "." + comparer.category().toString().toLowerCase() + TSV_EXTENSION, comparer);
            }
        }

        if(mCurationCache != null)
        {
            List<String> outputFiles = Lists.newArrayList(detailedFileComparers.keySet());
            outputFiles.add(combinedFile);

            if(mConfig.WriteTypes.contains(WriteType.FIELD_CONFIG))
                outputFiles.add(FieldConfigFile.generateFileName(filePrefix));

            for(String outputFile : outputFiles)
            {
                if(mCurationCache.isCacheFile(outputFile))
                {
                    CMP_LOGGER.error("output file({}) is also the curation cache file", outputFile);
                    return false;
                }
            }
        }

        try
        {
            for(Map.Entry<String,ItemComparer> entry : detailedFileComparers.entrySet())
            {
                String detailedFile = entry.getKey();
                ItemComparer comparer = entry.getValue();

                CMP_LOGGER.debug("writing output results: {}", detailedFile);

                BufferedWriter writer = createBufferedWriter(detailedFile, false);

                writer.write(MismatchFile.commonHeader(false, mCurationCache != null));

                List<String> compareFields = comparer.displayFieldNames();

                for(String field : compareFields)
                {
                    writer.write(String.format("\tOld%s\tNew%s", field, field));
                }

                writer.newLine();
                mCategoryWriters.put(comparer.category(), writer);
            }

            if(mConfig.WriteTypes.contains(WriteType.GENERIC))
            {
                CMP_LOGGER.debug("writing output results: {}", combinedFile);

                mCombinedWriter = createBufferedWriter(combinedFile, false);

                boolean writeSourceLabels = mOldSourceLabel != null;

                mCombinedWriter.write(MismatchFile.header(writeSourceLabels, mCurationCache != null));
                mCombinedWriter.newLine();
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to initialise Compar output files: {}", e.toString());
            return false;
        }

        return true;
    }

    public void close()
    {
        closeBufferedWriter(mCombinedWriter);
        mCategoryWriters.values().forEach(x -> closeBufferedWriter(x));
    }

    public synchronized void writeSampleMismatches(final String sampleId, final ItemComparer comparer, final List<Mismatch> mismatches)
    {
        if(mismatches.isEmpty())
            return;

        checkRemoveIgnoredGenes(mismatches);

        if(mCategoryWriters.isEmpty() && mCombinedWriter == null)
            return;

        try
        {
            CategoryType category = comparer.category();

            BufferedWriter categoryWriter = mCategoryWriters.get(category);

            for(Mismatch mismatch : mismatches)
            {
                List<FieldDisplayInfo> fieldDisplayValues = extractDisplayData(comparer, mismatch);

                List<CurationCache.Match> matches = mCurationCache != null ?
                        mCurationCache.match(sampleId, mismatch) : List.of(new CurationCache.Match(mismatch, null, null));

                for(CurationCache.Match match : matches)
                {
                    writeMismatch(sampleId, match, fieldDisplayValues, categoryWriter);
                }
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to write sample data: {}", e.toString());
        }
    }

    private List<FieldDisplayInfo> extractDisplayData(final ItemComparer comparer, final Mismatch mismatch)
    {
        List<FieldDisplayInfo> fieldDisplayValues = Lists.newArrayList();

        List<String> fieldNames = comparer.displayFieldNames();

        Map<String,FieldValue> oldFieldValues = mismatch.OldItem != null ? mismatch.OldItem.fieldValues() : Collections.emptyMap();
        Map<String,FieldValue> newFieldValues = mismatch.NewItem != null ? mismatch.NewItem.fieldValues() : Collections.emptyMap();

        for(String fieldName : fieldNames)
        {
            FieldValue oldValue = oldFieldValues.get(fieldName);
            FieldValue newValue = newFieldValues.get(fieldName);

            fieldDisplayValues.add(new FieldDisplayInfo(
                    fieldName,
                    oldValue != null ? oldValue.displayValue() : null,
                    newValue != null ? newValue.displayValue() : null));
        }

        return fieldDisplayValues;
    }

    private void writeMismatch(
            final String sampleId, final CurationCache.Match match, final List<FieldDisplayInfo> fieldDisplayValues,
            final BufferedWriter categoryWriter) throws IOException
    {
        Mismatch mismatch = match.mismatch();
        String statusStr = match.status() != null ? match.status().toString() : "";

        if(mCombinedWriter != null)
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(sampleId);

            String curationStatus = mOldSourceLabel != null || mCurationCache != null ? statusStr : null;
            sj.add(MismatchFile.toTsv(mismatch, false, fieldDisplayValues, curationStatus));

            if(mOldSourceLabel != null)
                sj.add(mOldSourceLabel).add(mNewSourceLabel);

            if(mCurationCache != null)
                sj.add(MismatchFile.historyTsv(match.previous()));

            mCombinedWriter.write(sj.toString());
            mCombinedWriter.newLine();
        }

        if(categoryWriter != null)
        {
            categoryWriter.write(String.format("%s\t", sampleId));

            String curationStatus = mCurationCache != null ? statusStr : null;
            categoryWriter.write(MismatchFile.toTsv(mismatch, true, fieldDisplayValues, curationStatus));
            categoryWriter.newLine();
        }
    }

    private void checkRemoveIgnoredGenes(final List<Mismatch> mismatches)
    {
        if(mConfig.IgnoreGenes.isEmpty())
            return;

        int index = 0;

        while(index < mismatches.size())
        {
            Mismatch mismatch = mismatches.get(index);

            ComparableItem item = mismatch.nonNullItem();

            if(!item.geneName().isEmpty() && mConfig.IgnoreGenes.contains(item.geneName()))
            {
                mismatches.remove(index);
            }
            else
            {
                ++index;
            }
        }
    }
}
