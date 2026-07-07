package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_CN_SEGS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_CONTAMINATION;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_FIT_METHOD;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_GENDER;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_GERM_ABS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_MS_INDELS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_MS_STATUS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_PLOIDY;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_PURITY;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_QC_STATUS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_SV_TMB;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_TINC_LEVEL;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_TMB;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_TMB_STATUS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_TML;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_TML_STATUS;
import static com.hartwig.hmftools.compar.purple.PurityData.FLD_UNS_CN_SEGS;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.DoubleField;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.IntField;
import com.hartwig.hmftools.compar.common.field.StringField;
import com.hartwig.hmftools.compar.common.field.StringListField;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class PurityComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public PurityComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category() { return PURITY; }

    @Override
    public List<Field> fields(MatchLevel matchLevel)
    {
        return List.of(
                new DoubleField(FLD_PURITY, i -> ((PurityData) i).Purity.bestFit().purity(), true, 0.04, null, "%.2f"),
                new DoubleField(FLD_PLOIDY, i -> ((PurityData) i).Purity.bestFit().ploidy(), true, 0.1, null, "%.2f"),
                new DoubleField(FLD_CONTAMINATION, i -> ((PurityData) i).Purity.qc().contamination(), true, 0.005, null, "%.4f"),
                new DoubleField(FLD_TMB, i -> ((PurityData) i).Purity.tumorMutationalBurdenPerMb(), true, 0.1, 0.05, "%.2f"),
                new DoubleField(FLD_MS_INDELS, i -> ((PurityData) i).Purity.microsatelliteIndelsPerMb(), true, 0.1, 0.05, "%.4f"),
                new IntField(FLD_TML, i -> ((PurityData) i).Purity.tumorMutationalLoad(), true, 1., 0.05),
                new IntField(FLD_CN_SEGS, i -> ((PurityData) i).Purity.qc().copyNumberSegments(), true, 5., 0.2),
                new IntField(FLD_UNS_CN_SEGS, i -> ((PurityData) i).Purity.qc().unsupportedCopyNumberSegments(), true, 5., 0.2),
                new IntField(FLD_SV_TMB, i -> ((PurityData) i).Purity.svTumorMutationalBurden(), true, 5., 0.05),
                new StringListField(FLD_QC_STATUS, i -> qcStatus(((PurityData) i)), true),
                new StringField(FLD_GENDER, i -> ((PurityData) i).Purity.gender().toString(), true),
                new StringListField(FLD_GERM_ABS, i -> germlineAberrations(((PurityData) i)), true),
                new StringField(FLD_FIT_METHOD, i -> ((PurityData) i).Purity.method().toString(), true),
                new StringField(FLD_MS_STATUS, i -> ((PurityData) i).Purity.microsatelliteStatus().toString(), true),
                new StringField(FLD_TMB_STATUS, i -> ((PurityData) i).Purity.tumorMutationalBurdenStatus().toString(), true),
                new StringField(FLD_TML_STATUS, i -> ((PurityData) i).Purity.tumorMutationalLoadStatus().toString(), true),
                new DoubleField(FLD_TINC_LEVEL, i -> ((PurityData) i).Purity.qc().tincLevel(), true, 0.1, null, "%.2f")
        );
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Lists.newArrayList(
                FLD_PURITY, FLD_PLOIDY, FLD_CONTAMINATION, FLD_TMB, FLD_TML, FLD_MS_INDELS, FLD_SV_TMB, FLD_CN_SEGS ,FLD_UNS_CN_SEGS,
                FLD_QC_STATUS, FLD_GENDER, FLD_GERM_ABS, FLD_FIT_METHOD, FLD_MS_STATUS, FLD_TMB_STATUS, FLD_TML_STATUS, FLD_TINC_LEVEL);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        PurityContext purityContext = dbAccess.readPurityContext(sampleId);

        List<ComparableItem> items = Lists.newArrayList();
        items.add(new PurityData(purityContext));
        return items;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            PurityContext purityContext = PurityContextFile.read(fileSources.Purple, sampleId);
            comparableItems.add(new PurityData(purityContext));

        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Purple purity data: {}", sampleId, e.toString());
            return null;
        }

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
