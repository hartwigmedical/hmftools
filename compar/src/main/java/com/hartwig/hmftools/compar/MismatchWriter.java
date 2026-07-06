package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.compar.MismatchFile.loadSampleCurations;
import static com.hartwig.hmftools.compar.common.CommonUtils.buildComparers;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CurationType.NONE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import javax.annotation.Nullable;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.KnownMismatch;
import com.hartwig.hmftools.compar.common.CurationInfo;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.WriteType;
import com.hartwig.hmftools.compar.common.field.Field;

public class MismatchWriter
{
    private final ComparConfig mConfig;
    private BufferedWriter mCombinedWriter;
    private final Map<CategoryType,BufferedWriter> mCategoryWriters;

    private final Map<String,List<KnownMismatch>> mSampleKnownMismatches;
    private final boolean mWriteCurations;
    private final boolean mWriteCurationComment;

    public MismatchWriter(final ComparConfig config)
    {
        mConfig = config;
        mCombinedWriter = null;
        mCategoryWriters = Maps.newHashMap();

        String sampleId = mConfig.SampleIds.size() == 1 ? mConfig.SampleIds.get(0) : null;

        mSampleKnownMismatches = config.KnownMismatchFile != null ?
                loadSampleCurations(config.KnownMismatchFile, sampleId) : Collections.emptyMap();

        mWriteCurations = !mSampleKnownMismatches.isEmpty();
        mWriteCurationComment = mSampleKnownMismatches.values().stream().anyMatch(x -> x.stream().anyMatch(y -> !y.Curation.Comment.isEmpty()));
    }

    public boolean initialiseOutputFiles()
    {
        String filePrefix = mConfig.OutputDir;

        if(mConfig.singleSample())
            filePrefix += mConfig.SampleIds.get(0) + ".cmp";
        else
            filePrefix += "compar_cohort";

        if(mConfig.OutputId != null)
            filePrefix += "." + mConfig.OutputId;

        try
        {
            if(mConfig.WriteTypes.contains(WriteType.TYPE_SPECIFIC))
            {
                List<ItemComparer> comparers = buildComparers(mConfig);

                for(ItemComparer comparer : comparers)
                {
                    String detailedFile = filePrefix + "." + comparer.category().toString().toLowerCase() + TSV_EXTENSION;

                    CMP_LOGGER.debug("writing output results: {}", detailedFile);

                    BufferedWriter writer = createBufferedWriter(detailedFile, false);

                    writer.write(MismatchFile.commonHeader(mConfig.multiSample(), false));

                    List<String> compareFields = comparer.displayFieldNames();

                    for(String field : compareFields)
                    {
                        writer.write(String.format("\tOld%s\tNew%s", field, field));
                    }

                    writer.newLine();
                    mCategoryWriters.put(comparer.category(), writer);
                }
            }

            if(mConfig.WriteTypes.contains(WriteType.GENERIC))
            {
                String outputFile = filePrefix + TSV_EXTENSION;

                CMP_LOGGER.debug("writing output results: {}", outputFile);

                mCombinedWriter = createBufferedWriter(outputFile, false);

                mCombinedWriter.write(MismatchFile.header(mConfig.multiSample(), mWriteCurations, mWriteCurationComment));
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
        // checkRemoveExpectedMismatches(sampleId, mismatches);

        if(mCategoryWriters.isEmpty() && mCombinedWriter == null)
            return;

        try
        {
            List<Field> displayFields = determineDisplayFields(comparer);
            CategoryType category = comparer.category();

            BufferedWriter categoryWriter = mCategoryWriters.get(category);

            List<KnownMismatch> knownMismatches = mSampleKnownMismatches.getOrDefault(sampleId, Collections.emptyList());

            for(Mismatch mismatch : mismatches)
            {
                // check or any expected mismatches / curations
                Map<String,CurationInfo> matchCurations = KnownMismatch.matchCurations(mismatch, knownMismatches);

                if(!matchCurations.isEmpty() && !mismatch.DiffValues.isEmpty() && matchCurations.size() != mismatch.DiffValues.size())
                {
                    // need to write separate lines to distinguish between the matches
                    for(String diff : mismatch.DiffValues)
                    {
                        Mismatch singleMismatch = new Mismatch(mismatch.OldItem, mismatch.NewItem, mismatch.Type, List.of(diff));
                        writeMismatch(sampleId, singleMismatch, displayFields, categoryWriter, matchCurations.get(diff));
                    }
                }
                else
                {
                    CurationInfo curationInfo = !matchCurations.isEmpty() ? matchCurations.entrySet().iterator().next().getValue() : null;
                    writeMismatch(sampleId, mismatch, displayFields, categoryWriter, curationInfo);
                }
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to write sample data: {}", e.toString());
        }
    }

    private void writeMismatch(
            final String sampleId, final Mismatch mismatch, final List<Field> displayFields, final BufferedWriter categoryWriter,
            @Nullable final CurationInfo curationInfo) throws IOException
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        if(sampleId != null && mConfig.multiSample())
        {
            sj.add(sampleId);

            if(categoryWriter != null)
            {
                categoryWriter.write(String.format("%s\t", sampleId));
            }
        }

        if(mCombinedWriter != null)
        {
            sj.add(MismatchFile.toTsv(mismatch, false, displayFields));

            if(mWriteCurations)
            {
                if(curationInfo != null)
                {
                    sj.add(curationInfo.Type.toString());

                    if(mWriteCurationComment)
                        sj.add(curationInfo.Comment);
                }
                else
                {
                    sj.add(NONE.toString());

                    if(mWriteCurationComment)
                        sj.add("");
                }
            }

            mCombinedWriter.write(sj.toString());
            mCombinedWriter.newLine();
        }

        if(categoryWriter != null)
        {
            categoryWriter.write(MismatchFile.toTsv(mismatch, true, displayFields));
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

    private static List<Field> determineDisplayFields(final ItemComparer comparer)
    {
        Map<String, Field> fieldNameToField = comparer.fields().stream().collect(Collectors.toMap(Field::name, f -> f));

        return comparer.displayFieldNames().stream().map(fieldNameToField::get).collect(Collectors.toList());
    }
}
