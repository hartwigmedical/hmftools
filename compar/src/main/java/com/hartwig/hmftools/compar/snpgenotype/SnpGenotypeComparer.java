package com.hartwig.hmftools.compar.snpgenotype;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.SNP_GENOTYPE;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;
import static com.hartwig.hmftools.compar.snpgenotype.SnpGenotypeData.FLD_GENOTYPE;
import static com.hartwig.hmftools.compar.snpgenotype.SnpGenotypeData.FLD_VCF_SAMPLE_ID;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.variant.variantcontext.VariantContext;

public class SnpGenotypeComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    private static final String FILE_NAME = "snp_genotype_output.vcf";

    public SnpGenotypeComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category()
    {
        return SNP_GENOTYPE;
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_ALT, FLD_GENOTYPE, FLD_VCF_SAMPLE_ID);
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
        String vcfFile = checkAddDirSeparator(fileSources.SnpGenotype) + FILE_NAME;

        VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

        if(!vcfFileReader.fileValid())
        {
            CMP_LOGGER.error("failed to read SNP genotype VCF file({})", vcfFile);
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

                items.add(new SnpGenotypeData(chromosome, position, ref, alt, genotype, vcfSampleId, comparisonPosition));
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
