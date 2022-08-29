package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.GERMLINE_SV;
import static com.hartwig.hmftools.compar.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.linx.GermlineSvData.FLD_GERMLINE_FRAGS;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.LinxGermlineSv;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.Mismatch;
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
    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_REPORTED, FLD_GERMLINE_FRAGS);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        // final List<StructuralVariantData> svDataList = dbAccess.readStructuralVariantGe(sampleId);
        // currently unsupported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        try
        {
            final List<LinxGermlineSv> germlineSvs = LinxGermlineSv.read(LinxGermlineSv.generateFilename(fileSources.LinxGermline, sampleId));

            CMP_LOGGER.debug("sample({}) loaded {} germline SVs", sampleId, germlineSvs.size());

            return germlineSvs.stream().map(x -> new GermlineSvData(x)).collect(Collectors.toList());
        }
        catch(IOException e)
        {
            CMP_LOGGER.info("sample({}) failed to load Linx germline SV data: {}", sampleId, e.toString());
            return Lists.newArrayList();
        }
    }
}
