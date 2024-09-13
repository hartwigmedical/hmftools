package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.common.Category.GERMLINE_SV;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;
import static com.hartwig.hmftools.compar.linx.GermlineSvData.FLD_GERMLINE_FRAGS;
import static com.hartwig.hmftools.compar.linx.LinxCommon.FLD_JUNCTION_COPY_NUMBER;
import static com.hartwig.hmftools.compar.linx.LinxCommon.FLD_UNDISRUPTED_COPY_NUMBER;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxGermlineSv;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class GermlineSvComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public GermlineSvComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return GERMLINE_SV; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_QUAL, 20, 0.2);
        thresholds.addFieldThreshold(FLD_GERMLINE_FRAGS, 5, 0.1);
        thresholds.addFieldThreshold(FLD_JUNCTION_COPY_NUMBER, 0.5, 0.2);
        thresholds.addFieldThreshold(FLD_UNDISRUPTED_COPY_NUMBER, 0.5, 0.2);
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        List<String> fieldNames = LinxCommon.comparedFieldNamesBreakends();
        fieldNames.add(FLD_GERMLINE_FRAGS);
        fieldNames.add(FLD_QUAL);
        return fieldNames;
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        // final List<StructuralVariantData> svDataList = dbAccess.readStructuralVariantGe(sampleId);
        // currently unsupported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        List<ComparableItem> items = Lists.newArrayList();

        try
        {
            String germlineSvFile = LinxGermlineSv.generateFilename(fileSources.LinxGermline, sampleId);
            List<LinxGermlineSv> germlineSvs = LinxGermlineSv.read(germlineSvFile);

            String germlineBreakendFile = LinxBreakend.generateFilename(fileSources.LinxGermline, sampleId, true);

            List<LinxBreakend> germlineBreakends = LinxBreakend.read(germlineBreakendFile);

            CMP_LOGGER.debug("sample({}) loaded {} germline breakends", sampleId, germlineBreakends.size());

            for(LinxGermlineSv germlineSv : germlineSvs)
            {
                List<LinxBreakend> matchingBreakends =
                        germlineBreakends.stream().filter(x -> x.svId() == germlineSv.SvId).collect(Collectors.toList());
                for(LinxBreakend breakend : matchingBreakends)
                {
                    germlineBreakends.remove(breakend);

                    BasePosition comparisonStartPosition = determineComparisonGenomePosition(
                            germlineSv.ChromosomeStart, germlineSv.PositionStart, fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);
                    BasePosition comparisonEndPosition = determineComparisonGenomePosition(
                            germlineSv.ChromosomeEnd, germlineSv.PositionEnd, fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                    items.add(new GermlineSvData(germlineSv, breakend, comparisonStartPosition, comparisonEndPosition));
                }
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Linx germline SV data: {}", sampleId, e.toString());
            return null;
        }

        return items;
    }
}
