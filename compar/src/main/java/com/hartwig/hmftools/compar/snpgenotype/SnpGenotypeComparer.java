package com.hartwig.hmftools.compar.snpgenotype;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CategoryType.SNP_GENOTYPE;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;
import static com.hartwig.hmftools.compar.common.FieldCheckCache.getOrMakeFieldCheck;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.variant.variantcontext.VariantContext;

public class SnpGenotypeComparer extends ItemComparer
{
    private static final String FILE_NAME = "snp_genotype_output.vcf";

    protected static final String FLD_GENOTYPE = "Genotype";
    protected static final String FLD_VCF_SAMPLE_ID = "VcfSampleId";

    protected enum Fields
    {
        Alt,
        Genotype,
        VcfSampleId;
    }

    public SnpGenotypeComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.Alt.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Alt.toString()), null));

        mFields.add(new FieldInfo(
                Fields.Genotype.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Genotype.toString()), null));

        mFields.add(new FieldInfo(
                Fields.VcfSampleId.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.VcfSampleId.toString()), null));
    }

    @Override
    public CategoryType category() { return SNP_GENOTYPE; }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        String vcfFile = checkAddDirSeparator(fileSources.SnpGenotype) + FILE_NAME;

        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

        if(!vcfFileReader.fileValid())
        {
            CMP_LOGGER.warn("failed to read SNP genotype VCF file({})", vcfFile);
            return null;
        }

        final List<ComparableItem> items = Lists.newArrayList();
        try (CloseableTribbleIterator<VariantContext> variantReader = vcfFileReader.iterator())
        {
            for(VariantContext variantContext : variantReader)
            {
                String chromosome = variantContext.getContig();
                int position = variantContext.getStart();
                String ref = variantContext.getReference().getBaseString();
                String alt = !variantContext.getAlternateAlleles().isEmpty() ? variantContext.getAlternateAlleles().get(0).toString() : ".";

                List<String> vcfSampleIds = variantContext.getSampleNamesOrderedByName();
                if(vcfSampleIds.size() != 1)
                {
                    throw new RuntimeException("sample(" + sampleId + ") SNPcheck VCF has more than one sample ID: " + vcfSampleIds);
                }
                String vcfSampleId = vcfSampleIds.get(0);
                String genotype = variantContext.getGenotype(vcfSampleId).getType().name();

                BasePosition comparisonPosition = determineComparisonGenomePosition(
                        chromosome, position, fileSources.Source, mConfig.RequiresLiftover, mConfig.LiftoverCache);

                items.add(new SnpGenotypeData(chromosome, position, ref, alt, genotype, vcfSampleId, comparisonPosition, mFields));
            }
        }
        catch(Exception e)
        {
            CMP_LOGGER.warn("sample({}) failed to load SNP genotype data: {}", sampleId, e.toString());
            return null;
        }

        return items;
    }

}
