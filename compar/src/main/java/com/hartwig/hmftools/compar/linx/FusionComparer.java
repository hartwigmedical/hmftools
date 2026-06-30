package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.common.CategoryType.FUSION;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.linx.DisruptionComparer.buildBreakendData;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_CHAIN_LINKS;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_CHAIN_TERM;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_DOMAINS_KEPT;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_DOMAINS_LOST;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_EXON_DOWN;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_EXON_UP;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_JUNCTION_COPY_NUMBER;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_LIKELIHOOD;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_PHASED;
import static com.hartwig.hmftools.compar.linx.FusionData.FLD_REPORTED_TYPE;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class FusionComparer implements ItemComparer
{
    private final ComparConfig mConfig;
    private DisruptionComparer mDisruptionComparer;

    public FusionComparer(final ComparConfig config)
    {
        mConfig = config;
        mDisruptionComparer = null;
    }

    public void setDisruptionComparer(final DisruptionComparer disruptionComparer) { mDisruptionComparer = disruptionComparer; }

    @Override
    public CategoryType category() { return FUSION; }

    @Override
    public void registerThresholds(final FieldConfig fieldConfig)
    {
        fieldConfig.addFieldThreshold(category(), FLD_JUNCTION_COPY_NUMBER, 0.5, 0.2);
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(
                FLD_REPORTED, FLD_REPORTED_TYPE, FLD_PHASED, FLD_LIKELIHOOD, FLD_EXON_UP,
                FLD_EXON_DOWN, FLD_CHAIN_LINKS, FLD_CHAIN_TERM, FLD_DOMAINS_KEPT, FLD_DOMAINS_LOST);

        // excluded unless matching breakends can be loaded: FLD_TRANSCRIPT_UP, FLD_TRANSCRIPT_DOWN, FLD_JUNCTION_COPY_NUMBER
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        List<LinxFusion> fusions = dbAccess.readFusions(sampleId);
        return processFusions(fusions, sourceType);
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        try
        {
            List<LinxFusion> fusions = LinxFusion.read(LinxFusion.generateFilename(fileSources.Linx, sampleId));
            return processFusions(fusions, fileSources.Source);
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Linx fusion data: {}", sampleId, e.toString());
            return null;
        }
    }

    private List<ComparableItem> processFusions(final List<LinxFusion> fusions, final SourceType sourceType)
    {
        List<ComparableItem> comparableItems = Lists.newArrayList();

        List<LinxBreakend> breakends = mDisruptionComparer != null ?
                mDisruptionComparer.breakends().get(sourceType) : Collections.emptyList();


        List<StructuralVariantData> svDataList = mDisruptionComparer != null ?
                mDisruptionComparer.svDataList().get(sourceType) : Collections.emptyList();

        for(LinxFusion fusion : fusions)
        {
            if(fusion.reportedType().equals(KnownFusionType.NONE.toString()))
                continue;

            String[] geneNames = LinxFusion.geneNames(fusion.name());

            BreakendData breakendStart = buildBreakend(fusion.fivePrimeBreakendId(), geneNames[0], breakends, svDataList);
            BreakendData breakendEnd = buildBreakend(fusion.threePrimeBreakendId(), geneNames[1], breakends, svDataList);

            comparableItems.add(new FusionData(fusion, fusion.name(), breakendStart, breakendEnd));
        }

        return comparableItems;
    }

    private BreakendData buildBreakend(
            final int breakendId, final String geneName, final List<LinxBreakend> breakends, final List<StructuralVariantData> svDataList)
    {
        LinxBreakend breakend = breakends.stream()
                .filter(x -> x.id() == breakendId && x.gene().equals(geneName)).findFirst().orElse(null);

        if(breakend == null)
        {
            breakend = breakends.stream().filter(x -> x.id() == breakendId).findFirst().orElse(null);

            if(breakend == null)
                return null;
        }

        int svId = breakend.svId();

        StructuralVariantData svData = svDataList.stream().filter(x -> x.id() == svId).findFirst().orElse(null);

        if(svData == null)
            return null;

        return buildBreakendData(breakend, svData, null);
    }
}
