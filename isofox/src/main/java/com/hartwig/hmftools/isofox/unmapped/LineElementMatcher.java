package com.hartwig.hmftools.isofox.unmapped;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.INFERRED;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.sv.linx.LinxCluster;
import com.hartwig.hmftools.common.sv.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;

public class LineElementMatcher
{
    private final String mLinxDir;
    private final String mSvVcfFile;

    private final Map<String,List<StructuralVariant>> mUmrMatchedVariants; // UMR position key to matched SVs

    private static final int MAX_DISTANCE = 10000;

    public LineElementMatcher(final String linxDir, final String svVcfFile)
    {
        mLinxDir = linxDir;
        mSvVcfFile = svVcfFile;
        mUmrMatchedVariants = Maps.newHashMap();
    }

    public String formUmrMatchString(final String posKey)
    {
        List<StructuralVariant> matchedSVs = getUmrMatch(posKey);

        if(matchedSVs == null)
            return "";

        StringJoiner sj = new StringJoiner(ITEM_DELIM);

        // VcfId:Type:Position:OtherPos
        matchedSVs.forEach(x -> sj.add(String.format("%s:%s:%d:%d",
                x.id(), x.type(), x.position(true), x.end() != null ? x.position(false) : 0)));

        return sj.toString();
    }

    public List<StructuralVariant> getUmrMatch(final String posKey)
    {
        return mUmrMatchedVariants.get(posKey);
    }

    public void findMatches(
            final String sampleId, final Map<String,Map<String, Map<String,List<UnmappedRead>>>> unmappedReads)
    {
        List<StructuralVariant> variants = loadSampleLinxData(sampleId);
        mUmrMatchedVariants.clear();

        for(Map<String, Map<String, List<UnmappedRead>>> positionMaps : unmappedReads.values())
        {
            for(Map.Entry<String, Map<String, List<UnmappedRead>>> umrSampleMap : positionMaps.entrySet())
            {
                List<UnmappedRead> umReads = umrSampleMap.getValue().get(sampleId);

                findMatch(umReads.get(0), variants);
            }
        }
    }

    private void findMatch(final UnmappedRead umRead, final List<StructuralVariant> variants)
    {
        for(StructuralVariant variant : variants)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                if((variant.type() == SGL || variant.type() == INF) && se == SE_END)
                    continue;

                StructuralVariantLeg leg = se == SE_START ? variant.start() : variant.end();

                if(!umRead.ReadRegion.Chromosome.equals(leg.chromosome()))
                    continue;

                if(abs(leg.position() - umRead.ExonBoundary) <= MAX_DISTANCE)
                {
                    List<StructuralVariant> matchedSVs = mUmrMatchedVariants.get(umRead);
                    if(matchedSVs == null)
                    {
                        matchedSVs = Lists.newArrayList();
                        mUmrMatchedVariants.put(umRead.positionKey(), matchedSVs);
                    }

                    matchedSVs.add(variant);
                }
            }
        }
    }

    private List<StructuralVariant> loadSampleLinxData(final String sampleId)
    {
        List<StructuralVariant> variants = Lists.newArrayList();

        String sampleLinxDir = mLinxDir.contains("*") ? mLinxDir.replaceAll("\\*", sampleId) : mLinxDir;
        String svVcfFile = mSvVcfFile.contains("*") ? mSvVcfFile.replaceAll("\\*", sampleId) : mSvVcfFile;

        try
        {
            List<Integer> lineClusterIds = LinxCluster.read(LinxCluster.generateFilename(sampleLinxDir, sampleId)).stream()
                    .filter(x -> x.resolvedType().equals("LINE"))
                    .map(x -> x.clusterId())
                    .collect(Collectors.toList());

            List<String> lineVcfIds = LinxSvAnnotation.read(LinxSvAnnotation.generateFilename(sampleLinxDir, sampleId)).stream()
                    .filter(x -> lineClusterIds.contains(x.clusterId()))
                    .map(x -> x.vcfId())
                    .collect(Collectors.toList());

            StructuralVariantFileLoader.fromFile(svVcfFile, new AlwaysPassFilter()).stream()
                    .filter(x -> lineVcfIds.contains(x.id()))
                    .filter(x -> x.filter().isEmpty() || x.filter().equals(PASS) || x.filter().equals(INFERRED))
                    .forEach(x -> variants.add(x));

            ISF_LOGGER.info("sample({}) loaded {} line variants", sampleId, variants.size());
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to read Linx files for sample({}): {}", sampleId, e.toString());
        }

        return variants;
    }
}
