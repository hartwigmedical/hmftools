package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

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

public class PurityComparer implements ItemComparer
{
    public static final String FLD_PURITY = PurityFields.Purity.toString();
    public static final String FLD_PLOIDY = PurityFields.Ploidy.toString();

    /*
    protected static final String FLD_CONTAMINATION = "Contamination";
    protected static final String FLD_TMB = "TmbPerMb";
    protected static final String FLD_MS_INDELS = "MsIndelsPerMb";
    protected static final String FLD_TML = "Tml";
    protected static final String FLD_CN_SEGS = "CopyNumberSegments";
    protected static final String FLD_UNS_CN_SEGS = "UnsupportedCopyNumberSegments";
    protected static final String FLD_SV_TMB = "SvTmb";
    protected static final String FLD_QC_STATUS = "QcStatus";
    protected static final String FLD_GENDER = "Gender";
    protected static final String FLD_GERM_ABS = "GermlineAberrations";
    protected static final String FLD_FIT_METHOD = "FitMethod";
    protected static final String FLD_MS_STATUS = "MsStatus";
    protected static final String FLD_TMB_STATUS = "TmbStatus";
    protected static final String FLD_TML_STATUS = "TmlStatus";
    protected static final String FLD_TINC_LEVEL = "TincLevel";
    */

    protected enum PurityFields
    {
        Purity(true, 0.04, null, "%.2f"),
        Ploidy(true, 0.1, null, "%.2f"),
        Contamination(true, 0.005, null, "%.4f"),
        TmbPerMb(true, 0.1, 0.05, "%.2f"),
        MsIndelsPerMb(true, 0.1, 0.05, "%.2f"),
        Tml(true, 1.0, 0.05, null),
        CopyNumberSegments(true, 5.0, 0.2, null),
        UnsupportedCopyNumberSegments(true, 5.0, 0.2, null),
        SvTmb(true, 5.0, 0.05, null),
        QcStatus(true, null, null, null),
        Gender(true, null, null, null),
        GermlineAberrations(true, null, null, null),
        FitMethod(true, null, null, null),
        MsStatus(true, null, null, null),
        TmbStatus(true, null, null, null),
        TmlStatus(true, null, null, null),
        TincLevel(true, 0.1, null, "%.3f");

        public final boolean IsCompared;
        public final Double DefaultAbsoluteThreshold;
        public final Double DefaultPercentageThreshold;
        public final String FormatString;

        PurityFields(
                final boolean isCompared, final Double defaultAbsoluteThreshold, final Double defaultPercentageThreshold,
                final String formatString)
        {
            IsCompared = isCompared;
            DefaultAbsoluteThreshold = defaultAbsoluteThreshold;
            DefaultPercentageThreshold = defaultPercentageThreshold;
            FormatString = formatString;
        }
    }

    private final ComparConfig mConfig;

    public PurityComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category() { return PURITY; }

    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        /*
        return List.of(
                new DoubleField(FLD_PURITY, i -> ((PurityData) i).Purity.bestFit().purity(),
                        true, 0.04, null, "%.2f"),
                new DoubleField(FLD_PLOIDY, i -> ((PurityData) i).Purity.bestFit().ploidy(),
                        true, 0.1, null, "%.2f"),
                new DoubleField(FLD_CONTAMINATION, i -> ((PurityData) i).Purity.qc().contamination(),
                        true, 0.005, null, "%.4f"),
                new DoubleField(FLD_TMB, i -> ((PurityData) i).Purity.tumorMutationalBurdenPerMb(),
                        true, 0.1, 0.05, "%.2f"),
                new DoubleField(FLD_MS_INDELS, i -> ((PurityData) i).Purity.microsatelliteIndelsPerMb(),
                        true, 0.1, 0.05, "%.2f"),
                new IntField(FLD_TML, i -> ((PurityData) i).Purity.tumorMutationalLoad(),
                        true, 1., 0.05),
                new IntField(FLD_CN_SEGS, i -> ((PurityData) i).Purity.qc().copyNumberSegments(),
                        true, 5., 0.2),
                new IntField(FLD_UNS_CN_SEGS, i -> ((PurityData) i).Purity.qc().unsupportedCopyNumberSegments(),
                        true, 5., 0.2),
                new IntField(FLD_SV_TMB, i -> ((PurityData) i).Purity.svTumorMutationalBurden(),
                        true, 5., 0.05),
                new StringListField(FLD_QC_STATUS, i -> qcStatus(((PurityData) i)), true),
                new StringField(FLD_GENDER, i -> ((PurityData) i).Purity.gender().toString(), true),
                new StringListField(FLD_GERM_ABS, i -> germlineAberrations(((PurityData) i)), true),
                new StringField(FLD_FIT_METHOD, i -> ((PurityData) i).Purity.method().toString(), true),
                new StringField(FLD_MS_STATUS, i -> ((PurityData) i).Purity.microsatelliteStatus().toString(), true),
                new StringField(FLD_TMB_STATUS, i -> ((PurityData) i).Purity.tumorMutationalBurdenStatus().toString(),
                        true),
                new StringField(FLD_TML_STATUS, i -> ((PurityData) i).Purity.tumorMutationalLoadStatus().toString(),
                        true),
                new DoubleField(FLD_TINC_LEVEL, i -> ((PurityData) i).Purity.qc().tincLevel(),
                        true, 0.1, null, "%.3f")
        );
        */

        return Collections.emptyList();
    }

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
            comparableItems.add(new PurityData(purityContext, fieldCheckMap));

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
        comparableItems.add(new PurityData(truthsetValues, fieldCheckMap));
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
