package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.isofox.RnaFusionData.FLD_JUNC_TYPE_DOWN;
import static com.hartwig.hmftools.compar.isofox.RnaFusionData.FLD_JUNC_TYPE_UP;
import static com.hartwig.hmftools.compar.isofox.RnaFusionData.FLD_KNOWN_TYPE;
import static com.hartwig.hmftools.compar.isofox.RnaFusionData.FLD_SPLIT_FRAGS;

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
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.IntField;
import com.hartwig.hmftools.compar.common.field.StringField;
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
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return List.of(
                new StringField(FLD_KNOWN_TYPE, i -> ((RnaFusionData) i).RnaFusion().knownType().toString(), true),
                new StringField(FLD_JUNC_TYPE_UP, i -> ((RnaFusionData) i).RnaFusion().junctionTypeUp(), true),
                new StringField(FLD_JUNC_TYPE_DOWN, i -> ((RnaFusionData) i).RnaFusion().junctionTypeDown(), true),
                new IntField(FLD_SPLIT_FRAGS, i -> ((RnaFusionData) i).RnaFusion().splitFragments(), true, 5., 0.05)
        );
    }

    @Override
    public List<String> displayFieldNames()
    {
        return List.of(FLD_KNOWN_TYPE, FLD_JUNC_TYPE_UP, FLD_JUNC_TYPE_DOWN, FLD_SPLIT_FRAGS);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
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
        String filename = RnaFusionFile.generateFilename(fileSources.Isofox, sampleId);
        String oldFilename = filename.replace(".tsv", ".csv");

        if(!Files.exists(Paths.get(filename)) && Files.exists(Paths.get(oldFilename)))
        {
            return oldFilename;
        }
        else
        {
            return filename;
        }
    }
}
