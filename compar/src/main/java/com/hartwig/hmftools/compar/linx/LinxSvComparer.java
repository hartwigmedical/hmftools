package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SV_VCF_SUFFIX;
import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.INFERRED;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.compar.Category.LINX_DATA;
import static com.hartwig.hmftools.compar.Category.SV;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.MismatchType.PRESENCE;

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
        final List<Mismatch> svMismatches = Lists.newArrayList();

        CommonUtils.processSample(this, mConfig, sampleId, svMismatches);

        if(svMismatches.isEmpty())
            return;

        for(Mismatch mismatch : svMismatches)
        {
            if(mismatch.MismatchType == PRESENCE)
            {
                if(mConfig.Categories.containsKey(SV))
                {
                    mismatches.add(new Mismatch(
                            mismatch.RefItem, mismatch.OtherItem, PRESENCE, Lists.newArrayList()));
                }
            }
            else if(mConfig.Categories.containsKey(LINX_DATA))
            {
                mismatches.add(mismatch);
            }
        }
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess)
    {
        final List<LinxCluster> clusters = dbAccess.readClusters(sampleId);
        final List<StructuralVariantData> svDataList = dbAccess.readStructuralVariantData(sampleId);
        final List<LinxSvAnnotation> annotations = dbAccess.readSvAnnotations(sampleId);

        return buildComparableItems(svDataList, annotations, clusters);
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

            CMP_LOGGER.debug("sample({}) loaded {} SVs, {} annotations and {} clusters",
                    sampleId, svDataList.size(), annotations.size(), clusters.size());

            List<ComparableItem> comparableItems = buildComparableItems(svDataList, annotations, clusters);

            CMP_LOGGER.debug("sample({}) loaded {} annotated SVs", sampleId, comparableItems.size());

            return comparableItems;
        }
        catch(IOException e)
        {
            CMP_LOGGER.info("sample({}) failed to load Linx SV and cluster data: {}", sampleId, e.toString());
            return Lists.newArrayList();
        }
    }

    private List<ComparableItem> buildComparableItems(
            final List<StructuralVariantData> svDataList, final List<LinxSvAnnotation> annotations, final List<LinxCluster> clusters)
    {
        final List<ComparableItem> svItems = Lists.newArrayList();

        for(StructuralVariantData svData : svDataList)
        {
            if(!svData.filter().isEmpty() && !svData.filter().equals(PASS) && !svData.filter().equals(INFERRED))
                continue;

            LinxSvAnnotation annotation = annotations.stream().filter(x -> x.vcfId().equals(svData.vcfId())).findFirst().orElse(null);

            if(annotation == null)
            {
                CMP_LOGGER.debug("sv({}:{}) missing Linx annotation", svData.vcfId(), svData.id());
                continue;
            }

            LinxCluster cluster = clusters.stream().filter(x -> x.clusterId() == annotation.clusterId()).findFirst().orElse(null);

            if(cluster == null)
            {
                CMP_LOGGER.debug("sv({}:{}) missing Linx cluster", svData.vcfId(), annotation.svId());
                continue;
            }

            svItems.add(new LinxSvData(svData, annotation, cluster));

            annotations.remove(annotation);
        }

        if(CMP_LOGGER.isDebugEnabled())
        {
            for(LinxSvAnnotation annotation : annotations)
            {
                CMP_LOGGER.debug("sv({}:{}) unmatched Linx annotation", annotation.vcfId(), annotation.svId());
            }
        }

        return svItems;
    }
}
