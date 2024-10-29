package com.hartwig.hmftools.esvee.depth;

import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.HighDepthRegion;
import com.hartwig.hmftools.common.region.UnmappedRegions;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class UnmappedRegionDepth
{
    private final Map<String,List<HighDepthRegion>> mChrLocationsMap;

    public UnmappedRegionDepth(final String unmappedRegionsFile)
    {
        if(unmappedRegionsFile != null && Files.exists(Paths.get(unmappedRegionsFile)))
            mChrLocationsMap = UnmappedRegions.loadUnmapRegions(unmappedRegionsFile);
        else
            mChrLocationsMap = Collections.emptyMap();
    }

    public void setUnmappedRegionsDepth(final int sampleCount, final List<DepthTask> depthTasks)
    {
        // work out median ref depth for each sample
        List<List<Integer>> sampleVariantDepths = Lists.newArrayList();

        for(int sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex)
        {
            sampleVariantDepths.add(Lists.newArrayList());
        }

        List<VariantContext> unmappedLowRefDepthVarints = Lists.newArrayList();

        for(DepthTask depthTask : depthTasks)
        {
            for(VariantContext variantContext : depthTask.variants())
            {
                for(int sampleIndex = 0; sampleIndex < variantContext.getNSamples(); ++sampleIndex)
                {
                    Genotype genotype = variantContext.getGenotype(sampleIndex);
                    int refDepth = getGenotypeAttributeAsInt(genotype, REF_DEPTH, 0);
                    int refPairDepth = getGenotypeAttributeAsInt(genotype, REF_DEPTH_PAIR, 0);
                    int tumorFrags = getGenotypeAttributeAsInt(genotype, TOTAL_FRAGS, 0);
                    double af = tumorFrags / (double)(refDepth + refPairDepth + tumorFrags);

                    if(af < 0.9)
                    {
                        sampleVariantDepths.get(sampleIndex).add(refDepth);
                        continue;
                    }

                    if(inHighDepthUnmappedRegion(variantContext.getContig(), variantContext.getStart()))
                    {
                        unmappedLowRefDepthVarints.add(variantContext);
                    }
                }
            }
        }

        if(unmappedLowRefDepthVarints.isEmpty())
            return;

        List<Integer> sampleMedianDepths = Lists.newArrayList();

        for(int sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex)
        {
            List<Integer> depthValues = sampleVariantDepths.get(sampleIndex);

            if(depthValues.size() < 2)
            {
                sampleMedianDepths.add(0);
                continue;
            }

            Collections.sort(depthValues);
            int medianIndex = depthValues.size() / 2;
            int medianDepth = depthValues.get(medianIndex);

            SV_LOGGER.info("sample median depth({}) from {} variants", medianDepth, depthValues.size());

            sampleMedianDepths.add(medianDepth);
        }

        SV_LOGGER.info("setting ref depth to median for {} variants", unmappedLowRefDepthVarints.size());

        for(VariantContext variantContext : unmappedLowRefDepthVarints)
        {
            int totalRefDepth = 0;

            for(int sampleIndex = 0; sampleIndex < variantContext.getNSamples(); ++sampleIndex)
            {
                Genotype genotype = variantContext.getGenotype(sampleIndex);
                int refDepth = getGenotypeAttributeAsInt(genotype, REF_DEPTH, 0);
                int medianDepth = sampleMedianDepths.get(sampleIndex);

                if(refDepth < medianDepth)
                {
                    totalRefDepth += medianDepth;
                    genotype.getExtendedAttributes().put(REF_DEPTH, medianDepth);
                }
                else
                {
                    totalRefDepth += refDepth;
                }
            }

            variantContext.getCommonInfo().removeAttribute(REF_DEPTH);
            variantContext.getCommonInfo().putAttribute(REF_DEPTH, totalRefDepth);
        }
    }

    private boolean inHighDepthUnmappedRegion(final String chromosome, final int position)
    {
        List<HighDepthRegion> regions = mChrLocationsMap.get(chromosome);

        if(regions == null)
            return false;

        return regions.stream().filter(x -> x.maxDepth() > 0).anyMatch(x -> x.containsPosition(position));

    }
}
