package com.hartwig.hmftools.wisp.purity.variant;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.common.CommonUtils.FLD_CATEGORY;
import static com.hartwig.hmftools.wisp.common.CommonUtils.FLD_TUMOR_ID;
import static com.hartwig.hmftools.wisp.common.CommonUtils.FLD_VARIANT;
import static com.hartwig.hmftools.wisp.probe.CategoryType.GERMLINE_MUTATION;
import static com.hartwig.hmftools.wisp.probe.CategoryType.OTHER_CLONAL_MUTATION;
import static com.hartwig.hmftools.wisp.probe.CategoryType.OTHER_CODING_MUTATION;
import static com.hartwig.hmftools.wisp.probe.CategoryType.OTHER_MUTATION;
import static com.hartwig.hmftools.wisp.probe.CategoryType.REPORTABLE_MUTATION;
import static com.hartwig.hmftools.wisp.probe.CategoryType.SUBCLONAL_MUTATION;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.wisp.probe.CategoryType;

public class ProbeVariantCache
{
    private final Map<String, List<ProbeVariant>> mSampleVariants;

    public ProbeVariantCache(final String filename)
    {
        mSampleVariants = Maps.newHashMap();

        if(filename != null)
            loadVariants(filename);
    }

    public List<ProbeVariant> getSampleVariants(final String sampleId) { return mSampleVariants.get(sampleId); }

    private void loadVariants(final String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String header = lines.get(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, CSV_DELIM);

            int sampleIndex = fieldsIndexMap.get(FLD_TUMOR_ID);
            int categoryIndex = fieldsIndexMap.get(FLD_CATEGORY);
            int variantIndex = fieldsIndexMap.get(FLD_VARIANT);

            lines.remove(0);

            List<ProbeVariant> sampleVariants = null;
            String currentSampleId = "";

            for(String line : lines)
            {
                String[] values = line.split(CSV_DELIM, -1);

                String sampleId = values[sampleIndex];

                if(!sampleId.equals(currentSampleId))
                {
                    currentSampleId = sampleId;
                    sampleVariants = Lists.newArrayList();
                    mSampleVariants.put(sampleId, sampleVariants);
                }

                CategoryType category = CategoryType.valueOf(values[categoryIndex]);

                if(category == REPORTABLE_MUTATION || category == GERMLINE_MUTATION || category == OTHER_CODING_MUTATION
                || category == OTHER_CLONAL_MUTATION || category == OTHER_MUTATION || category == SUBCLONAL_MUTATION)
                {
                    String variantStr = values[variantIndex];

                    // parse into details: chr3:6289637 A>C SNP
                    String[] parts = variantStr.split(" ", 3);

                    String[] coords = parts[0].split(":");
                    String[] mutation = parts[1].split(">", 2);
                    VariantType type = VariantType.valueOf(parts[2]);

                    sampleVariants.add(new ProbeVariant(
                            coords[0], Integer.parseInt(coords[1]), mutation[0], mutation[1], type));
                }
            }

            CT_LOGGER.info("loaded {} sample probe variants from file({})",
                    mSampleVariants.values().stream().mapToInt(x -> x.size()).sum(), filename);
        }
        catch(Exception e)
        {
            CT_LOGGER.error("failed to load sample probe variants file({}): {}", filename, e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }
}
