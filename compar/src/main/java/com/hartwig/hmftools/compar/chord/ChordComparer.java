package com.hartwig.hmftools.compar.chord;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.common.CategoryType.CHORD;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.fileExists;
import static com.hartwig.hmftools.compar.common.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class ChordComparer extends ItemComparer
{
    private static final String OLD_CHORD_FILE_EXTENSION = "_chord_prediction.txt";

    protected enum Fields
    {
        BRCA1,
        BRCA2,
        Score,
        Type,
        Status;
    }

    public ChordComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.BRCA1.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.BRCA1.toString(), 0.1, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.BRCA2.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.BRCA2.toString(), 0.1, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.Score.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.Score.toString(), 0.1, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.Type.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Type.toString()), null));

        mFields.add(new FieldInfo(
                Fields.Status.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Status.toString()), null));
    }

    @Override
    public CategoryType category() { return CHORD; }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            ChordData chordData = ChordDataFile.read(determineChordFilePath(sampleId, fileSources));
            comparableItems.add(new ChordComparData(chordData, mFields));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Chord data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }

    private static String determineChordFilePath(final String sampleId, final PipelineSourcePaths fileSources)
    {
        final String currentFilePath = ChordDataFile.generateFilename(fileSources.Chord, sampleId);
        final String oldFilePath = checkAddDirSeparator(fileSources.Chord) + sampleId + OLD_CHORD_FILE_EXTENSION;

        if(!fileExists(currentFilePath) && fileExists(oldFilePath))
        {
            return oldFilePath;
        }
        else
        {
            return currentFilePath;
        }
    }
}
