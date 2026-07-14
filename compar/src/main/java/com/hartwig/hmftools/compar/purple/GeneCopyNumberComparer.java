package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.CategoryType.GENE_COPY_NUMBER;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
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
import com.hartwig.hmftools.compar.common.field.DisplayOnlyField;
import com.hartwig.hmftools.compar.common.field.DoubleField;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class GeneCopyNumberComparer implements ItemComparer
{
    protected static final String FLD_MIN_COPY_NUMBER = "MinCopyNumber";
    protected static final String FLD_MAX_COPY_NUMBER = "MaxCopyNumber";
    protected static final String FLD_MIN_REGION_START = "MinRegionStart";
    protected static final String FLD_MIN_REGION_END = "MinRegionEnd";

    private final ComparConfig mConfig;

    public GeneCopyNumberComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category() { return GENE_COPY_NUMBER; }

    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return List.of(
                new DoubleField(FLD_MIN_COPY_NUMBER, i -> ((GeneCopyNumberData) i).CopyNumber.minCopyNumber(),
                        true, 0.5, 0.15, "%.2f"),
                new DoubleField(FLD_MAX_COPY_NUMBER, i -> ((GeneCopyNumberData) i).CopyNumber.maxCopyNumber(),
                        true, 0.5, 0.15, "%.2f"),
                new DisplayOnlyField(FLD_MIN_REGION_START, i -> String.valueOf(((GeneCopyNumberData) i).CopyNumber.MinRegionStart), i -> true),
                new DisplayOnlyField(FLD_MIN_REGION_END, i -> String.valueOf(((GeneCopyNumberData) i).CopyNumber.MinRegionEnd), i -> true)
        );
    }

    @Override
    public boolean hasReportable() { return false; }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches, final FieldConfig fieldConfig)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches, fieldConfig);
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Lists.newArrayList(FLD_MIN_COPY_NUMBER, FLD_MAX_COPY_NUMBER, FLD_MIN_REGION_START, FLD_MIN_REGION_END);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        List<ComparableItem> items = Lists.newArrayList();

        final Set<String> driverGenes = mConfig.DriverGenes;

        if(driverGenes.isEmpty())
            return items;

        final List<GeneCopyNumber> copyNumbers = dbAccess.readGeneCopynumbers(sampleId, new ArrayList<>(driverGenes));

        copyNumbers.forEach(x -> items.add(new GeneCopyNumberData(x)));
        return items;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        final List<ComparableItem> items = Lists.newArrayList();

        // load only the genes in the driver catalog if the driver category is being run

        final Set<String> driverGenes = mConfig.DriverGenes;

        if(driverGenes.isEmpty())
            return items;

        try
        {
            List<GeneCopyNumber> copyNumbers = GeneCopyNumberFile.read(GeneCopyNumberFile.generateFilename(
                    fileSources.Purple, sampleId));

            copyNumbers.stream().filter(x -> driverGenes.contains(x.geneName())).forEach(x -> items.add(new GeneCopyNumberData(x)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to read gene copy number data: {}", sampleId, e.toString());
            return null;
        }

        return items;
    }
}
