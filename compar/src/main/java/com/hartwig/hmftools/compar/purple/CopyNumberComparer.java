package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.CategoryType.COPY_NUMBER;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;
import static com.hartwig.hmftools.compar.purple.CopyNumberData.FLD_COPY_NUMBER;
import static com.hartwig.hmftools.compar.purple.CopyNumberData.FLD_MAJOR_ALLELE_CN;
import static com.hartwig.hmftools.compar.purple.CopyNumberData.FLD_METHOD;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.region.BasePosition;
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
import com.hartwig.hmftools.compar.common.field.StringField;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class CopyNumberComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public CopyNumberComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType category() { return COPY_NUMBER; }

    @Override
    public List<Field> fields(MatchLevel matchLevel)
    {
        return List.of(
                new DoubleField(FLD_COPY_NUMBER, i -> ((CopyNumberData) i).copyNumber(), true, 0.5, 0.15, "%.2f"),
                new DoubleField(FLD_MAJOR_ALLELE_CN, i -> ((CopyNumberData) i).majorAlleleCopyNumber(), true, 0.5, 0.15, "%.2f"),
                new StringField(FLD_METHOD, i -> ((CopyNumberData) i).method().toString(), true)
        );
    }

    @Override
    public boolean hasReportable() { return false; }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Lists.newArrayList(FLD_COPY_NUMBER, FLD_MAJOR_ALLELE_CN, FLD_METHOD);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        final List<PurpleCopyNumber> copyNumbers = dbAccess.readCopynumbers(sampleId);
        List<ComparableItem> items = Lists.newArrayList();
        copyNumbers.forEach(x -> items.add(createCopyNumberData(x, sourceType)));
        return items;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilenameForReading(
                    fileSources.Purple, sampleId));

            copyNumbers.forEach(x -> comparableItems.add(createCopyNumberData(x, fileSources.Source)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to read copy number data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }

    private CopyNumberData createCopyNumberData(final PurpleCopyNumber copyNumber, final SourceType sourceType)
    {
        BasePosition comparisonPositionStart = determineComparisonGenomePosition(
                copyNumber.chromosome(), copyNumber.start(), sourceType, mConfig.RequiresLiftover, mConfig.LiftoverCache);

        BasePosition comparisonPositionEnd = determineComparisonGenomePosition(
                copyNumber.chromosome(), copyNumber.end(), sourceType, mConfig.RequiresLiftover, mConfig.LiftoverCache);

        return new CopyNumberData(
                copyNumber.chromosome(), copyNumber.start(), copyNumber.end(),
                copyNumber.averageTumorCopyNumber(), copyNumber.majorAlleleCopyNumber(),
                copyNumber.method(), comparisonPositionStart, comparisonPositionEnd
        );
    }
}
