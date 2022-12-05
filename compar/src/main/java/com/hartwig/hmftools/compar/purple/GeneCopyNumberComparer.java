package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.Category.GENE_COPY_NUMBER;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.purple.CopyNumberData.FLD_COPY_NUMBER;
import static com.hartwig.hmftools.compar.purple.CopyNumberData.FLD_MAJOR_ALLELE_CN;
import static com.hartwig.hmftools.compar.purple.CopyNumberData.FLD_METHOD;
import static com.hartwig.hmftools.compar.purple.GeneCopyNumberData.FLD_MAX_COPY_NUMBER;
import static com.hartwig.hmftools.compar.purple.GeneCopyNumberData.FLD_MIN_COPY_NUMBER;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class GeneCopyNumberComparer implements ItemComparer
{
    private final ComparConfig mConfig;
    private final List<String> mDriverGeneNames;

    public GeneCopyNumberComparer(final ComparConfig config)
    {
        mConfig = config;
        mDriverGeneNames = mConfig.DriverGenes.stream().collect(Collectors.toList());
    }

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
        return Lists.newArrayList(FLD_COPY_NUMBER, FLD_MAJOR_ALLELE_CN, FLD_METHOD);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        List<ComparableItem> items = Lists.newArrayList();

        if(mDriverGeneNames.isEmpty())
            return items;

        final List<GeneCopyNumber> copyNumbers = dbAccess.readGeneCopynumbers(sampleId, mDriverGeneNames);
        copyNumbers.forEach(x -> items.add(new GeneCopyNumberData(x)));
        return items;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> items = Lists.newArrayList();

        if(mDriverGeneNames.isEmpty())
            return items;

        try
        {
            List<GeneCopyNumber> copyNumbers = GeneCopyNumberFile.read(GeneCopyNumberFile.generateFilenameForReading(
                    fileSources.Purple, sampleId));

            copyNumbers.stream().filter(x -> mConfig.DriverGenes.contains(x.geneName())).forEach(x -> items.add(new GeneCopyNumberData(x)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to read gene copy number data: {}", sampleId, e.toString());
            return null;
        }

        return items;
    }
}
