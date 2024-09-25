package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.common.Category.GERMLINE_SV;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;
import static com.hartwig.hmftools.compar.linx.DisruptionData.FLD_BREAKEND_INFO;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxGermlineDisruption;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.MatchLevel;
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
    public void registerThresholds(final DiffThresholds thresholds) {}

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_REPORTED, FLD_BREAKEND_INFO);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        // currently unsupported
        return Lists.newArrayList();
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        List<ComparableItem> items = Lists.newArrayList();

        Map<String,List<BreakendData>> geneBreakendMap = Maps.newHashMap();

        MatchLevel matchLevel = mConfig.Categories.getOrDefault(GERMLINE_SV, MatchLevel.REPORTABLE);

        try
        {
            String germlineSvFile = LinxGermlineDisruption.generateFilename(fileSources.LinxGermline, sampleId);
            List<LinxGermlineDisruption> germlineSvs = null; // loads on demand

            String germlineBreakendFile = LinxBreakend.generateFilename(fileSources.LinxGermline, sampleId, true);

            List<LinxBreakend> germlineBreakends = LinxBreakend.read(germlineBreakendFile).stream()
                    .filter(x -> x.reportedDisruption()).collect(Collectors.toList());

            boolean reportedOnly = mConfig.Categories.get(GERMLINE_SV) == MatchLevel.REPORTABLE;

            for(LinxBreakend breakend : germlineBreakends)
            {
                if(reportedOnly && !breakend.reportedDisruption())
                    continue;

                if(germlineSvs == null)
                {
                    germlineSvs = LinxGermlineDisruption.read(germlineSvFile);
                    CMP_LOGGER.debug("sample({}) loaded {} germline SVs", sampleId, germlineSvs.size());
                }

                LinxGermlineDisruption var = germlineSvs.stream().filter(x -> x.SvId == breakend.svId()).findFirst().orElse(null);

                if(var == null)
                    continue; // implies an inconsistency

                List<BreakendData> geneBreakends = geneBreakendMap.get(breakend.gene());

                if(geneBreakends == null)
                {
                    geneBreakends = Lists.newArrayList();
                    geneBreakendMap.put(breakend.gene(), geneBreakends);
                }

                boolean usesStart = breakend.isStart();

                int[] homologyOffsets = {0, 0}; // not available at the moment unless VCF is read

                String chromosome = usesStart ? var.ChromosomeStart : var.ChromosomeEnd;
                int position = usesStart ? var.PositionStart : var.PositionEnd;

                BasePosition comparisonPosition = determineComparisonGenomePosition(
                        chromosome, position, fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                BreakendData breakendData = new BreakendData(
                        breakend, var.VcfId, var.Type,
                        comparisonPosition.Chromosome, comparisonPosition.Position,
                        usesStart ? var.OrientStart : var.OrientEnd, homologyOffsets);

                geneBreakends.add(breakendData);
            }

            for(Map.Entry<String,List<BreakendData>> entry : geneBreakendMap.entrySet())
            {
                String geneName = entry.getKey();
                List<BreakendData> geneBreakends = entry.getValue();

                DisruptionData disruptionData = new DisruptionData(GERMLINE_SV, geneName, geneBreakends);
                items.add(disruptionData);
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
