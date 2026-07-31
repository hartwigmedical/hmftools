package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.FieldCheckCache.getOrMakeFieldCheck;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.FieldCheckCache;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.TruthsetValue;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class PurityComparer implements ItemComparer
{
    public static final String FLD_PURITY = PurityFields.Purity.toString();
    public static final String FLD_PLOIDY = PurityFields.Ploidy.toString();

    private static List<FieldInfo> mFields;

    protected enum PurityFields
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

    private final ComparConfig mConfig;

    public PurityComparer(final ComparConfig config, final Map<String,FieldCheck> fieldCheckMap)
    {
        mConfig = config;

        mFields = Lists.newArrayList();

        mFields.add(new FieldInfo(
                PurityFields.Purity.toString(),
                getOrMakeFieldCheck(fieldCheckMap, PurityFields.Purity.toString(), 0.04, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                PurityFields.Ploidy.toString(),
                getOrMakeFieldCheck(fieldCheckMap, PurityFields.Ploidy.toString(), 0.1, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                PurityFields.Contamination.toString(),
                getOrMakeFieldCheck(fieldCheckMap, PurityFields.Contamination.toString(), 0.005, null),
                "%.4f"));

        mFields.add(new FieldInfo(
                PurityFields.TmbPerMb.toString(),
                getOrMakeFieldCheck(fieldCheckMap, PurityFields.TmbPerMb.toString(), 0.1, 0.05),
                "%.2f"));

        mFields.add(new FieldInfo(
                PurityFields.MsIndelsPerMb.toString(),
                getOrMakeFieldCheck(fieldCheckMap, PurityFields.MsIndelsPerMb.toString(), 0.1, 0.05),
                "%.2f"));

        mFields.add(new FieldInfo(
                PurityFields.Tml.toString(),
                getOrMakeFieldCheck(fieldCheckMap, PurityFields.Tml.toString(), 1.0, 0.05),
                null));

        mFields.add(new FieldInfo(
                PurityFields.CopyNumberSegments.toString(),
                getOrMakeFieldCheck(fieldCheckMap, PurityFields.CopyNumberSegments.toString(), 5.0, 0.2),
                null));

        mFields.add(new FieldInfo(
                PurityFields.UnsupportedCopyNumberSegments.toString(),
                getOrMakeFieldCheck(fieldCheckMap, PurityFields.UnsupportedCopyNumberSegments.toString(), 5.0, 0.2),
                null));

        mFields.add(new FieldInfo(
                PurityFields.SvTmb.toString(),
                getOrMakeFieldCheck(fieldCheckMap, PurityFields.SvTmb.toString(), 5.0, 0.05),
                null));

        mFields.add(new FieldInfo(
                PurityFields.QcStatus.toString(),
                getOrMakeFieldCheck(fieldCheckMap, PurityFields.QcStatus.toString()), null));

        mFields.add(new FieldInfo(
                PurityFields.Gender.toString(),
                getOrMakeFieldCheck(fieldCheckMap, PurityFields.Gender.toString()),null));

        mFields.add(new FieldInfo(
                PurityFields.GermlineAberrations.toString(),
                getOrMakeFieldCheck(fieldCheckMap, PurityFields.GermlineAberrations.toString()), null));

        mFields.add(new FieldInfo(
                PurityFields.FitMethod.toString(), getOrMakeFieldCheck(fieldCheckMap, PurityFields.FitMethod.toString()),
                null));

        mFields.add(new FieldInfo(
                PurityFields.MsStatus.toString(), getOrMakeFieldCheck(fieldCheckMap, PurityFields.MsStatus.toString()),
                null));

        mFields.add(new FieldInfo(
                PurityFields.TmbStatus.toString(), getOrMakeFieldCheck(fieldCheckMap, PurityFields.TmbStatus.toString()),
                null));

        mFields.add(new FieldInfo(
                PurityFields.TmlStatus.toString(), getOrMakeFieldCheck(fieldCheckMap, PurityFields.TmlStatus.toString()),
                null));

        mFields.add(new FieldInfo(
                PurityFields.TincLevel.toString(),
                getOrMakeFieldCheck(fieldCheckMap, PurityFields.TincLevel.toString(), 0.1, null),
                "%.3f"));
    }

    @Override
    public CategoryType category() { return PURITY; }

    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return Collections.emptyList();
    }

    // TODO: rename to fields and make part of interface
    public List<FieldInfo> fieldsList() { return mFields; }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches, final FieldCheckCache fieldConfig)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches, fieldConfig);
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(PurityFields.values()).map(x -> x.toString()).collect(Collectors.toList());
        /*
        return Lists.newArrayList(
                FLD_PURITY, FLD_PLOIDY, FLD_CONTAMINATION, FLD_TMB, FLD_TML, FLD_MS_INDELS, FLD_SV_TMB, FLD_CN_SEGS ,FLD_UNS_CN_SEGS,
                FLD_QC_STATUS, FLD_GENDER, FLD_GERM_ABS, FLD_FIT_METHOD, FLD_MS_STATUS, FLD_TMB_STATUS, FLD_TML_STATUS, FLD_TINC_LEVEL);
        */
    }

    @Override
    public List<ComparableItem> loadFromFile(
            final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources,
            final Map<String,FieldCheck> fieldCheckMap)
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

    public List<ComparableItem> loadFromTruthset(final List<TruthsetValue> truthsetValues, final Map<String, FieldCheck> fieldCheckMap)
    {
        List<ComparableItem> comparableItems = Lists.newArrayList();
        comparableItems.add(new PurityData(truthsetValues, mFields));
        return comparableItems;
    }

    static List<String> germlineAberrations(final PurityData purityData)
    {
        return purityData.Purity.qc().germlineAberrations().stream()
                .map(a -> a.toString())
                .sorted()
                .toList();
    }

    static List<String> qcStatus(final PurityData purityData)
    {
        return purityData.Purity.qc().status().stream()
                .map(s -> s.toString())
                .sorted()
                .toList();
    }
}
