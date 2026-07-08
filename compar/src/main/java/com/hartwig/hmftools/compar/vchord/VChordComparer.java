package com.hartwig.hmftools.compar.vchord;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.V_CHORD;
import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_BREAST;
import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_OTHER;
import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_OVARIAN;
import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_PANCREATIC;
import static com.hartwig.hmftools.compar.vchord.VChordData.FLD_PROSTATE;

import java.io.UncheckedIOException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.vchord.VChordPrediction;
import com.hartwig.hmftools.common.vchord.VChordPredictionFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.DoubleField;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class VChordComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public VChordComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category() { return V_CHORD; }

    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return List.of(
                new DoubleField(FLD_BREAST, i -> ((VChordData) i).VChord().breastCancerHrdScore(), true, 0.1, null, "%.2f"),
                new DoubleField(FLD_OVARIAN, i -> ((VChordData) i).VChord().ovarianCancerHrdScore(), true, 0.1, null, "%.2f"),
                new DoubleField(FLD_PANCREATIC, i -> ((VChordData) i).VChord().pancreaticCancerScore(), true, 0.1, null, "%.2f"),
                new DoubleField(FLD_PROSTATE, i -> ((VChordData) i).VChord().prostateCancerScore(), true, 0.1, null, "%.2f"),
                new DoubleField(FLD_OTHER, i -> ((VChordData) i).VChord().otherCancerScore(), true, 0.1, null, "%.2f")
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
        return Lists.newArrayList(FLD_BREAST, FLD_OVARIAN, FLD_PANCREATIC, FLD_PROSTATE, FLD_OTHER);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        // currently unsupported
        return Collections.emptyList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();
        try
        {
            VChordPrediction vChordData = VChordPredictionFile.read(VChordPredictionFile.generateFilename(fileSources.VChord, sampleId));
            comparableItems.add(new VChordData(vChordData));
        }
        catch(UncheckedIOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load vChord data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }
}
