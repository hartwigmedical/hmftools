package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.rna.RnaFusionFile.FLD_DISCORD_FRAGS;
import static com.hartwig.hmftools.common.rna.RnaFusionFile.FLD_KNOWN_TYPE;
import static com.hartwig.hmftools.common.rna.RnaFusionFile.FLD_REALIGN_FRAGS;
import static com.hartwig.hmftools.common.rna.RnaFusionFile.FLD_SPLIT_FRAGS;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.isofox.RnaFusionData.FLD_JUNC_TYPE_DOWN;
import static com.hartwig.hmftools.compar.isofox.RnaFusionData.FLD_JUNC_TYPE_UP;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.rna.RnaFusionFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public record RnaFusionComparer(ComparConfig mConfig) implements ItemComparer
{
    @Override
    public CategoryType category()
    {
        return CategoryType.RNA_FUSION;
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
        thresholds.addFieldThreshold(FLD_SPLIT_FRAGS, 5, 0.05);
        thresholds.addFieldThreshold(FLD_REALIGN_FRAGS, 5, 0.05);
        thresholds.addFieldThreshold(FLD_DISCORD_FRAGS, 5, 0.05);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return List.of(FLD_KNOWN_TYPE, FLD_JUNC_TYPE_UP, FLD_JUNC_TYPE_DOWN, FLD_SPLIT_FRAGS, FLD_REALIGN_FRAGS, FLD_DISCORD_FRAGS);
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
        String filename = determineFileName(sampleId, fileSources);
        List<RnaFusion> fusions = RnaFusionFile.read(filename);
        if(fusions == null)
        {
            CMP_LOGGER.warn("sample({}) failed to load Isofox Gene data", sampleId);
            return null;
        }

        List<ComparableItem> comparableItems = Lists.newArrayList();

        for(RnaFusion fusion : fusions)
        {
            BasePosition comparisonPositionUp = CommonUtils.determineComparisonGenomePosition(
                    fusion.chromosomeUp(), fusion.positionUp(), fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

            BasePosition comparisonPositionDown = CommonUtils.determineComparisonGenomePosition(
                    fusion.chromosomeDown(), fusion.positionDown(), fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

            comparableItems.add(new RnaFusionData(fusion, comparisonPositionUp, comparisonPositionDown));
        }

        return comparableItems;
    }

    private static String determineFileName(final String sampleId, final FileSources fileSources)
    {
        String current_file_name = RnaFusionFile.generateFilename(fileSources.Isofox, sampleId);
        String old_file_name = current_file_name.replace(".tsv", ".csv");

        if(!Files.exists(Paths.get(current_file_name)) && Files.exists(Paths.get(old_file_name)))
        {
            return old_file_name;
        }
        else
        {
            return current_file_name;
        }
    }
}
