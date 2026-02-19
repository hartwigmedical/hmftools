package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_ALT_SJ_TYPE;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_FRAG_COUNT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_START;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public record NovelSpliceJunctionComparer(ComparConfig mConfig) implements ItemComparer
{
    @Override
    public CategoryType category()
    {
        return CategoryType.NOVEL_SPLICE_JUNCTION;
    }

    @Override
    public boolean hasReportable()
    {
        return false;
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_FRAG_COUNT, 5, 0.05);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return List.of(FLD_ALT_SJ_TYPE, FLD_FRAG_COUNT, FLD_REGION_START, FLD_REGION_END);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        // Not currently supported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        String filename = NovelSpliceJunctionFile.generateFilename(fileSources.Isofox, sampleId);
        List<NovelSpliceJunction> junctions = NovelSpliceJunctionFile.read(filename);
        if(junctions == null){
            CMP_LOGGER.warn("sample({}) failed to load Isofox Gene data", sampleId);
            return null;
        }

        List<ComparableItem> comparableItems = Lists.newArrayList();

        for(NovelSpliceJunction junction : junctions)
        {
            BasePosition comparisonPositionStart = CommonUtils.determineComparisonGenomePosition(
                    junction.chromosome(), junction.junctionStart(), fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

            BasePosition comparisonPositionEnd = CommonUtils.determineComparisonGenomePosition(
                    junction.chromosome(), junction.junctionEnd(), fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

            comparableItems.add(new NovelSpliceJunctionData(junction, comparisonPositionStart, comparisonPositionEnd));
        }

        return comparableItems;
    }
}
