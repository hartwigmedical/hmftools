package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SV_VCF_SUFFIX;
import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.compar.Category.LINX_DATA;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.sv.linx.LinxCluster;
import com.hartwig.hmftools.common.sv.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.FileSources;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class LinxSvComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public LinxSvComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return LINX_DATA; }

    @Override
    public void processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        final List<LinxCluster> clusters = dbAccess.readClusters(sampleId);
        final List<StructuralVariantData> svDataList = dbAccess.readStructuralVariantData(sampleId);
        final List<LinxSvAnnotation> annotations = dbAccess.readSvAnnotations(sampleId);

        return buildComparableItems(svDataList, annotations, clusters);
    }

    public List<ComparableItem> buildComparableItems(
            final List<StructuralVariantData> svDataList, final List<LinxSvAnnotation> annotations, final List<LinxCluster> clusters)
    {
        final List<ComparableItem> svItems = Lists.newArrayList();

        for(StructuralVariantData svData : svDataList)
        {
            if(!svData.filter().isEmpty() && !svData.filter().equals(PASS))
                continue;

            LinxSvAnnotation annotation = annotations.stream().filter(x -> x.svId() == svData.id()).findFirst().orElse(null);

            if(annotation == null)
                continue;

            LinxCluster cluster = clusters.stream().filter(x -> x.clusterId() == annotation.clusterId()).findFirst().orElse(null);

            if(cluster == null)
                continue;

            svItems.add(new LinxSvData(svData, annotation, cluster));
        }

        return svItems;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
    {
        try
        {
            List<LinxSvAnnotation> annotations = LinxSvAnnotation.read(LinxSvAnnotation.generateFilename(fileSources.Linx, sampleId));

            List<LinxCluster> clusters = LinxCluster.read(LinxCluster.generateFilename(fileSources.Linx, sampleId));

            final List<StructuralVariantData> svDataList = Lists.newArrayList();

            String vcfFile = fileSources.Purple + sampleId + PURPLE_SV_VCF_SUFFIX;
            List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(vcfFile, new AlwaysPassFilter());
            List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(variants);

            int svId = 0;

            for(EnrichedStructuralVariant variant : enrichedVariants)
            {
                svDataList.add(convertSvData(variant, svId++));
            }

            return buildComparableItems(svDataList, annotations, clusters);
        }
        catch(IOException e)
        {
            CMP_LOGGER.info("sample({}) failed to load Linx SV and cluster data: {}", sampleId, e.toString());
            return Lists.newArrayList();
        }
    }
}
