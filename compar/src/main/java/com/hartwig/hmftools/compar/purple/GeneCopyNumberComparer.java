package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.Category.GENE_COPY_NUMBER;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.purple.GeneCopyNumberData.FLD_MAX_COPY_NUMBER;
import static com.hartwig.hmftools.compar.purple.GeneCopyNumberData.FLD_MIN_COPY_NUMBER;
import static com.hartwig.hmftools.compar.purple.GeneCopyNumberData.FLD_MIN_REGION_END;
import static com.hartwig.hmftools.compar.purple.GeneCopyNumberData.FLD_MIN_REGION_START;

import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class GeneCopyNumberComparer implements ItemComparer
{
    private final ComparConfig mConfig;
    private final Set<String> mDriverGenes;

    public GeneCopyNumberComparer(final ComparConfig config)
    {
        mConfig = config;
        mDriverGenes = Sets.newHashSet();
    }

    public void addDriverGenes(final Set<String> driverGenes) { mDriverGenes.addAll(driverGenes); }

    @Override
    public Category category() { return GENE_COPY_NUMBER; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_MIN_COPY_NUMBER, 0.5, 0.15);
        thresholds.addFieldThreshold(FLD_MAX_COPY_NUMBER, 0.5, 0.15);
    }

    @Override
    public boolean hasReportable() { return false; }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_MIN_COPY_NUMBER, FLD_MAX_COPY_NUMBER, FLD_MIN_REGION_START, FLD_MIN_REGION_END);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        List<ComparableItem> items = Lists.newArrayList();

        final Set<String> driverGenes = !mDriverGenes.isEmpty() ? mDriverGenes : mConfig.DriverGenes;

        if(driverGenes.isEmpty())
            return items;

        final List<GeneCopyNumber> copyNumbers = dbAccess.readGeneCopynumbers(
                sampleId, driverGenes.stream().collect(Collectors.toList()));

        copyNumbers.forEach(x -> items.add(new GeneCopyNumberData(x)));
        return items;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> items = Lists.newArrayList();

        // load only the genes in the driver catalog if the driver category is being run

        final Set<String> driverGenes = !mDriverGenes.isEmpty() ? mDriverGenes : mConfig.DriverGenes;

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
