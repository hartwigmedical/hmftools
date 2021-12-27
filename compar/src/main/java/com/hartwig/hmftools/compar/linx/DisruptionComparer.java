package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.DISRUPTION;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class DisruptionComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public DisruptionComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return DISRUPTION; }

    @Override
    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        return dbAccess.readBreakends(sampleId).stream()
                .filter(x -> mConfig.DriverGenes.isEmpty() || mConfig.DriverGenes.contains(x.gene()))
                .map(x -> new DisruptionData(x)).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            List<LinxBreakend> breakends = LinxBreakend.read(LinxBreakend.generateFilename(fileSources.Linx, sampleId));

            breakends.stream()
                    .filter(x -> mConfig.DriverGenes.isEmpty() || mConfig.DriverGenes.contains(x.gene()))
                    .forEach(x -> comparableItems.add(new DisruptionData(x)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.info("sample({}) failed to load Linx breakend disruption data: {}", sampleId, e.toString());
        }

        return comparableItems;
    }}
