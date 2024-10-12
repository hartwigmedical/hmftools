package com.hartwig.hmftools.chord.sv;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;
import static com.hartwig.hmftools.chord.sv.SvContext.SV_TYPES;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.chord.ChordConfig;
import com.hartwig.hmftools.chord.common.LoggingOptions;
import com.hartwig.hmftools.chord.common.MutContextCount;
import com.hartwig.hmftools.chord.common.VariantTypePrep;
import com.hartwig.hmftools.chord.common.VcfFile;
import com.hartwig.hmftools.chord.indel.IndelPrep;
import com.hartwig.hmftools.common.sv.SvVcfTags;

import htsjdk.variant.variantcontext.VariantContext;

public class SvPrep implements VariantTypePrep<StructuralVariant>, LoggingOptions
{
    private final ChordConfig mConfig;
    private String mLogPrefix = "";

    private static final String SV_DETAILS_FILE_SUFFIX = ".chord.sv.details.tsv";

    public SvPrep(ChordConfig config)
    {
        mConfig = config;
    }

    @Override
    public SvPrep logPrefix(String logPrefix)
    {
        mLogPrefix = logPrefix;
        return this;
    }

    @Override
    public List<StructuralVariant> loadVariants(String sampleId) throws NoSuchFileException
    {
        VcfFile vcfFile = new VcfFile(mConfig.svVcfFile(sampleId), mConfig.IncludeNonPass).logPrefix(mLogPrefix);
        List<VariantContext> variantContexts = vcfFile.loadVariants();

        // Only keep the first mate of breakend pairs
        Set<String> idSet = new HashSet<>();
        List<StructuralVariant> variants = new ArrayList<>();
        for(VariantContext variantContext : variantContexts)
        {
            String id = variantContext.getID();
            String mateId = variantContext.getAttributeAsString(SvVcfTags.MATE_ID, id); // For SGLs, use id as mate id

            if(!idSet.contains(id))
            {
                variants.add(new StructuralVariant(variantContext));
                CHORD_LOGGER.trace("{}Selected SV with id({}) mateId({}): {} ", mLogPrefix, id, mateId, variantContext);
            }
            else
            {
                CHORD_LOGGER.trace("{}Skipped SV with id({}) mateId({}): {} ", mLogPrefix, id, mateId, variantContext);
            }

            idSet.add(id);
            idSet.add(mateId);
        }

        return variants
                .stream()
                .filter(v -> SV_TYPES.contains(v.Type))
                .collect(Collectors.toList());
    }

    @Override
    public List<MutContextCount> countMutationContexts(String sampleId)
    {
        try
        {
            CHORD_LOGGER.debug("{}Running {} - counting SVs by type and length", mLogPrefix, this.getClass().getSimpleName());

            List<StructuralVariant> svList = loadVariants(sampleId);
            CHORD_LOGGER.debug("{}Found {} SVs", mLogPrefix, svList.size());

            Map<String, Integer> contextCountsMap = SvContext.initializeCounts();

            for(StructuralVariant sv : svList)
            {
                SvContext svContext = SvContext.from(sv);
                String svContextName = svContext.getContextName();

                contextCountsMap.compute(svContextName, (k,v) -> v + 1);
            }

            List<MutContextCount> counts = new ArrayList<>();

            if(mConfig.WriteDetailedFiles)
            {
                String svDetailsPath = mConfig.OutputDir + "/" + sampleId + SV_DETAILS_FILE_SUFFIX;
                CHORD_LOGGER.debug("{}Writing indel details to: {}", mLogPrefix, svDetailsPath);
                writeDetails(svDetailsPath, svList);
            }

            for (String indelContextName: contextCountsMap.keySet()) {
                int count = contextCountsMap.get(indelContextName);

                MutContextCount mutTypeCount = new MutContextCount(indelContextName, count);

                CHORD_LOGGER.trace(mutTypeCount);

                counts.add(mutTypeCount);
            }

            CHORD_LOGGER.debug("{}Completed {}", mLogPrefix, this.getClass().getSimpleName());

            return counts;
        }
        catch(Exception e)
        {
            CHORD_LOGGER.error("{}{} failed: {}", mLogPrefix, this.getClass().getSimpleName(), e.toString());
            e.printStackTrace();
            return null;
        }
    }

    private static void writeDetails(String path, List<StructuralVariant> svList) throws IOException
    {
        BufferedWriter writer = SvDetails.initializeWriter(path);

        for(StructuralVariant sv : svList)
            SvDetails.writeLine(writer, sv);

        writer.close();
    }
}
