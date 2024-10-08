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
import com.hartwig.hmftools.chord.common.MutTypeCount;
import com.hartwig.hmftools.chord.common.VcfFile;
import com.hartwig.hmftools.common.sv.SvVcfTags;

import htsjdk.variant.variantcontext.VariantContext;

public class SvPrep
{
    private final ChordConfig mConfig;

    private static final String SV_DETAILS_FILE_SUFFIX = ".chord.sv.details.tsv";

    public SvPrep(ChordConfig config)
    {
        mConfig = config;
    }

    public List<StructuralVariant> loadVariants(String sampleId) throws NoSuchFileException
    {
        VcfFile vcfFile = new VcfFile(mConfig.purpleSvVcfFile(sampleId), mConfig.IncludeNonPass);
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
                CHORD_LOGGER.trace("  Added SV with id({}) mateId({}): {} ", id, mateId, variantContext);
            }
            else
            {
                CHORD_LOGGER.trace("Skipped SV with id({}) mateId({}): {} ", id, mateId, variantContext);
            }

            idSet.add(id);
            idSet.add(mateId);
        }

        return variants
                .stream()
                .filter(v -> SV_TYPES.contains(v.Type))
                .collect(Collectors.toList());
    }

    public List<MutTypeCount> extractSampleData(String sampleId)
    {
        try
        {
            CHORD_LOGGER.info("Extract SV type/length contexts");

            List<StructuralVariant> svList = loadVariants(sampleId);
            CHORD_LOGGER.debug("Loaded {} SVs", svList.size());

            CHORD_LOGGER.debug("Initializing counts");
            Map<String, Integer> contextCountsMap = SvContext.initializeCounts();

            CHORD_LOGGER.debug("Populating counts");
            for(StructuralVariant sv : svList)
            {
                SvContext svContext = SvContext.from(sv);
                String svContextName = svContext.getContextName();

                contextCountsMap.compute(svContextName, (k,v) -> v + 1);
            }

            List<MutTypeCount> counts = new ArrayList<>();

            if(mConfig.WriteDetailedFiles)
            {
                String svDetailsPath = mConfig.OutputDir + "/" + sampleId + SV_DETAILS_FILE_SUFFIX;
                CHORD_LOGGER.info("Writing SV details to: {}", svDetailsPath);
                writeDetails(svDetailsPath, svList);
            }

            for (String indelContextName: contextCountsMap.keySet()) {
                int count = contextCountsMap.get(indelContextName);

                MutTypeCount mutTypeCount = new MutTypeCount(indelContextName, count);

                CHORD_LOGGER.trace(mutTypeCount);

                counts.add(mutTypeCount);
            }

            return counts;
        }
        catch(Exception e)
        {
            CHORD_LOGGER.error("sample({}) failed to count SVs by type and length:", sampleId);
            e.printStackTrace();
            System.exit(1);
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
