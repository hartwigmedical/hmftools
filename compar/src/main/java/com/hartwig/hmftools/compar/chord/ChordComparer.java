package com.hartwig.hmftools.compar.chord;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.common.CategoryType.CHORD;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.fileExists;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.DoubleField;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.StringField;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class ChordComparer implements ItemComparer
{
    private static final String OLD_CHORD_FILE_EXTENSION = "_chord_prediction.txt";

    protected static final String FLD_BRCA1 = "BRCA1";
    protected static final String FLD_BRCA2 = "BRCA2";
    protected static final String FLD_STATUS = "Status";
    protected static final String FLD_TYPE = "Type";
    protected static final String FLD_SCORE = "Score";

    private final ComparConfig mConfig;

    public ChordComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category()
    {
        return CHORD;
    }

    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return List.of(
                new DoubleField(FLD_BRCA1, i -> ((ChordComparData) i).Chord.BRCA1Value(), true,
                        0.1, null, "%.2f"),
                new DoubleField(FLD_BRCA2, i -> ((ChordComparData) i).Chord.BRCA2Value(), true,
                        0.1, null, "%.2f"),
                new DoubleField(FLD_SCORE, i -> ((ChordComparData) i).Chord.hrdValue(), true,
                        0.1, null, "%.2f"),
                new StringField(FLD_TYPE, i -> ((ChordComparData) i).Chord.hrdType(), true),
                new StringField(FLD_STATUS, i -> ((ChordComparData) i).Chord.hrStatus().toString(), true)
        );
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches, final FieldConfig fieldConfig)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches, fieldConfig);
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Lists.newArrayList(FLD_BRCA1, FLD_BRCA2, FLD_SCORE, FLD_STATUS, FLD_TYPE);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        ChordData chordData = dbAccess.readChord(sampleId);
        final List<ComparableItem> comparableItems = Lists.newArrayList();
        comparableItems.add(new ChordComparData(chordData));
        return comparableItems;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            ChordData chordData = ChordDataFile.read(determineChordFilePath(sampleId, fileSources));
            comparableItems.add(new ChordComparData(chordData));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Chord data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }

    private static String determineChordFilePath(final String sampleId, final FileSources fileSources)
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
