package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.compar.common.Category.GERMLINE_SV;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;
import static com.hartwig.hmftools.compar.linx.GermlineSvData.FLD_GERMLINE_FRAGS;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
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
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_REPORTED, FLD_GERMLINE_FRAGS);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        // final List<StructuralVariantData> svDataList = dbAccess.readStructuralVariantGe(sampleId);
        // currently unsupported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        List<ComparableItem> items = Lists.newArrayList();

        try
        {
            String germlineSvFile = LinxGermlineSv.generateFilename(fileSources.LinxGermline, sampleId);
            List<LinxGermlineSv> germlineSvs = LinxGermlineSv.read(germlineSvFile);

            String germlineBreakendFile = LinxBreakend.generateFilename(fileSources.LinxGermline, sampleId, true);

            // germline breakend file was introduced in v5.32, for old versions extract reported from the SV file
            if(Files.exists(Paths.get(germlineBreakendFile)))
            {
                List<LinxBreakend> germlineBreakends = LinxBreakend.read(germlineBreakendFile).stream()
                        .filter(x -> x.reportedDisruption()).collect(Collectors.toList());

                CMP_LOGGER.debug("sample({}) loaded {} germline SVs", sampleId, germlineSvs.size());


                for(LinxGermlineSv germlineSv : germlineSvs)
                {
                    boolean isReported = germlineBreakends.stream().anyMatch(x -> x.svId() == germlineSv.SvId);

                    BasePosition comparisonStartPosition = determineComparisonGenomePosition(
                            germlineSv.ChromosomeStart, germlineSv.PositionStart, fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                    BasePosition comparisonEndPosition = determineComparisonGenomePosition(
                            germlineSv.ChromosomeEnd, germlineSv.PositionEnd, fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                    items.add(new GermlineSvData(germlineSv, isReported, comparisonStartPosition, comparisonEndPosition));
                }
            }
            else
            {
                List<String> rawGermlineSvs = Files.readAllLines(Paths.get(germlineSvFile));
                Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(rawGermlineSvs.get(0), TSV_DELIM);
                Integer reportedIndex = fieldsIndexMap.get("reported");

                if(reportedIndex == null)
                    return null;

                rawGermlineSvs.remove(0);

                for(int i = 0; i < germlineSvs.size(); ++i)
                {
                    LinxGermlineSv germlineSv = germlineSvs.get(i);
                    String[] values = rawGermlineSvs.get(i).split(TSV_DELIM, -1);
                    boolean isReported = Boolean.parseBoolean(values[reportedIndex]);

                    BasePosition comparisonPositionStart = determineComparisonGenomePosition(
                            germlineSv.ChromosomeStart, germlineSv.PositionStart, fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                    BasePosition comparisonPositionEnd = determineComparisonGenomePosition(
                            germlineSv.ChromosomeEnd, germlineSv.PositionEnd, fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                    items.add(new GermlineSvData(germlineSv, isReported, comparisonPositionStart, comparisonPositionEnd));
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
