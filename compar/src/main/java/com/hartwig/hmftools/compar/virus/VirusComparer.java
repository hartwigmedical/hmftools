package com.hartwig.hmftools.compar.virus;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.VIRUS;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.virus.AnnotatedVirusFile;
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
import com.hartwig.hmftools.compar.common.field.BooleanField;
import com.hartwig.hmftools.compar.common.field.DoubleField;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.IntField;
import com.hartwig.hmftools.compar.common.field.StringField;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class VirusComparer implements ItemComparer
{
    protected static final String FLD_INTEGRATIONS = "Integrations";
    protected static final String FLD_MEAN_COVERAGE = "MeanCoverage";
    protected static final String FLD_DRIVER_LIKELIHOOD = "DriverLikelihood";

    private final ComparConfig mConfig;

    public VirusComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category()
    {
        return VIRUS;
    }

    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return List.of(
                new BooleanField(FLD_REPORTED, i -> ((VirusData) i).Virus.reported(), true),
                new IntField(FLD_INTEGRATIONS, i -> ((VirusData) i).Virus.integrations(),
                        true, null, 0.20),
                new DoubleField(FLD_MEAN_COVERAGE, i -> ((VirusData) i).Virus.meanCoverage(),
                        true, null, 0.15, "%.2f"),
                new StringField(FLD_DRIVER_LIKELIHOOD, i -> String.valueOf(((VirusData) i).Virus.virusDriverLikelihoodType()),
                        true)
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
        return Lists.newArrayList(FLD_REPORTED, FLD_INTEGRATIONS, FLD_MEAN_COVERAGE, FLD_DRIVER_LIKELIHOOD);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        // currently unsupported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();
        try
        {
            AnnotatedVirusFile.read(AnnotatedVirusFile.generateFileName(fileSources.Virus, sampleId))
                    .forEach(v -> comparableItems.add(new VirusData(v)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Virus interpreter data: {}", sampleId, e.toString());
            return null;
        }
        return comparableItems;
    }

}
