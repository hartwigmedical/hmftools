package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.TruthsetValue;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class PurityComparer extends ItemComparer
{
    protected enum Fields
    {
        Purity,
        Ploidy,
        Contamination,
        TmbPerMb,
        MsIndelsPerMb,
        Tml,
        CopyNumberSegments,
        UnsupportedCopyNumberSegments,
        SvTmb,
        QcStatus,
        Gender,
        GermlineAberrations,
        FitMethod,
        MsStatus,
        TmbStatus,
        TmlStatus,
        TincLevel;
    }

    public PurityComparer(final ComparConfig config, final Map<String,FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.Purity.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.Purity.toString(), 0.04, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.Ploidy.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.Ploidy.toString(), 0.1, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.Contamination.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.Contamination.toString(), 0.005, null),
                "%.4f"));

        mFields.add(new FieldInfo(
                Fields.TmbPerMb.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.TmbPerMb.toString(), 0.1, 0.05),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.MsIndelsPerMb.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.MsIndelsPerMb.toString(), 0.1, 0.05),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.Tml.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.Tml.toString(), 1.0, 0.05),
                null));

        mFields.add(new FieldInfo(
                Fields.CopyNumberSegments.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.CopyNumberSegments.toString(), 5.0, 0.2),
                null));

        mFields.add(new FieldInfo(
                Fields.UnsupportedCopyNumberSegments.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.UnsupportedCopyNumberSegments.toString(), 5.0, 0.2),
                null));

        mFields.add(new FieldInfo(
                Fields.SvTmb.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.SvTmb.toString(), 5.0, 0.05),
                null));

        mFields.add(new FieldInfo(
                Fields.QcStatus.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.QcStatus.toString()), null));

        mFields.add(new FieldInfo(
                Fields.Gender.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.Gender.toString()),null));

        mFields.add(new FieldInfo(
                Fields.GermlineAberrations.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.GermlineAberrations.toString()), null));

        mFields.add(new FieldInfo(
                Fields.FitMethod.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.FitMethod.toString()),
                null));

        mFields.add(new FieldInfo(
                Fields.MsStatus.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.MsStatus.toString()),
                null));

        mFields.add(new FieldInfo(
                Fields.TmbStatus.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.TmbStatus.toString()),
                null));

        mFields.add(new FieldInfo(
                Fields.TmlStatus.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.TmlStatus.toString()),
                null));

        mFields.add(new FieldInfo(
                Fields.TincLevel.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.TincLevel.toString(), 0.1, null),
                "%.3f"));
    }

    @Override
    public CategoryType category() { return PURITY; }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(
            final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            PurityContext purityContext = PurityContextFile.read(fileSources.Purple, sampleId);
            comparableItems.add(new PurityData(purityContext, mFields));

        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Purple purity data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }

    public List<ComparableItem> loadFromTruthset(final Map<String,List<TruthsetValue>> valuesByKey)
    {
        List<ComparableItem> comparableItems = Lists.newArrayList();

        if(valuesByKey.size() == 1)
            comparableItems.add(new PurityData(valuesByKey.get(PURITY.toString()), mFields));

        return comparableItems;
    }
}
