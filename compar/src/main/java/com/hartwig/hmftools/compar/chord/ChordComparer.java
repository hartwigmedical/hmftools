package com.hartwig.hmftools.compar.chord;

import static com.hartwig.hmftools.compar.common.Category.CHORD;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.chord.ChordComparData.FLD_BRCA1;
import static com.hartwig.hmftools.compar.chord.ChordComparData.FLD_BRCA2;
import static com.hartwig.hmftools.compar.chord.ChordComparData.FLD_SCORE;
import static com.hartwig.hmftools.compar.chord.ChordComparData.FLD_STATUS;
import static com.hartwig.hmftools.compar.chord.ChordComparData.FLD_TYPE;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class ChordComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public ChordComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return CHORD; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_BRCA1, 0.1, 0);
        thresholds.addFieldThreshold(FLD_BRCA2, 0.1, 0);
        thresholds.addFieldThreshold(FLD_SCORE, 0.1, 0);
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_BRCA1, FLD_BRCA2, FLD_SCORE, FLD_STATUS, FLD_TYPE);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        ChordData chordData = dbAccess.readChord(sampleId);
        final List<ComparableItem> comparableItems = Lists.newArrayList();
        comparableItems.add(new ChordComparData(chordData));
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            ChordData chordData = ChordDataFile.read(ChordDataFile.generateFilename(fileSources.Chord, sampleId));
            comparableItems.add(new ChordComparData(chordData));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Chord data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }
}
