package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_ALT_SJ_TYPE;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_FRAG_COUNT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_START;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.nio.file.Files;
import java.nio.file.Paths;
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
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.IntField;
import com.hartwig.hmftools.compar.common.field.StringField;
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
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches, final FieldConfig fieldConfig)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches, fieldConfig);
    }

    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return List.of(
                new StringField(FLD_ALT_SJ_TYPE, i -> ((NovelSpliceJunctionData) i).NovelSpliceJunction().type().toString(), true),
                new IntField(FLD_FRAG_COUNT, i -> ((NovelSpliceJunctionData) i).NovelSpliceJunction().fragmentCount(), true, 5., 0.05),
                new StringField(FLD_REGION_START, i -> ((NovelSpliceJunctionData) i).NovelSpliceJunction().regionStart().toString(), true),
                new StringField(FLD_REGION_END, i -> ((NovelSpliceJunctionData) i).NovelSpliceJunction().regionEnd().toString(), true)
        );
    }

    @Override
    public List<String> displayFieldNames()
    {
        return List.of(FLD_ALT_SJ_TYPE, FLD_FRAG_COUNT, FLD_REGION_START, FLD_REGION_END);
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
        List<NovelSpliceJunction> junctions = NovelSpliceJunctionFile.read(filename);
        if(junctions == null)
        {
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

    private static String determineFileName(final String sampleId, final FileSources fileSources)
    {
        String filename = NovelSpliceJunctionFile.generateFilename(fileSources.Isofox, sampleId);
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
