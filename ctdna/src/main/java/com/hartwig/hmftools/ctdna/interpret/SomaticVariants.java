package com.hartwig.hmftools.ctdna.interpret;

import static com.hartwig.hmftools.common.variant.SageVcfTags.LIST_SEPARATOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariants
{
    private final InterpretConfig mConfig;
    private final ResultsWriter mResultsWriter;

    public SomaticVariants(final InterpretConfig config, final ResultsWriter resultsWriter)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
    }

    public void processPatientVcf(final String patientId, final String somaticVcf)
    {
        CT_LOGGER.info("patient({}) reading somatic VCF: {}", patientId, somaticVcf);

        AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                somaticVcf, new VCFCodec(), false);

        VCFHeader vcfHeader = (VCFHeader)reader.getHeader();
        List<GenotypeInfo> genotypeInfos = Lists.newArrayList();

        for(int i = 0; i < vcfHeader.getGenotypeSamples().size(); ++i)
        {
            genotypeInfos.add(new GenotypeInfo(i, vcfHeader.getGenotypeSamples().get(i)));
        }

        int variantCount = 0;

        try
        {
            for(VariantContext variantContext : reader.iterator())
            {
                if(variantContext.isFiltered())
                    continue;

                try
                {
                    processVariant(patientId, genotypeInfos, variantContext);
                }
                catch(Exception e)
                {
                    e.printStackTrace();
                    CT_LOGGER.error("patient({}) error processing variant", patientId, e.toString());
                    return;
                }
                ++variantCount;

                if(variantCount > 0 && (variantCount % 100000) == 0)
                {
                    CT_LOGGER.info("processed {} variants", variantCount);
                }
            }

        }
        catch(IOException e)
        {
            CT_LOGGER.error("error reading vcf files: {}", e.toString());
        }

        CT_LOGGER.info("patient({}) processed {} somatic variants", patientId, variantCount);

    }

    private void processVariant(final String patientId, final List<GenotypeInfo> genotypeInfos, final VariantContext variantContext)
    {
        VariantContextDecorator variant = new VariantContextDecorator(variantContext);

        for(GenotypeInfo genotypeInfo : genotypeInfos)
        {
            Genotype genotype = variantContext.getGenotype(genotypeInfo.Index);

            if(genotype == null || genotype.getExtendedAttributes().isEmpty())
            {
                //CT_LOGGER.warn("patientId({}) genotypeInfo({}) missing for variant({}:{}:{}>{})",
                //        patientId, genotypeInfo, variant.chromosome(), variant.position(), variant.ref(), variant.alt());
                continue;
            }

            int alleleCount = genotype.getAD()[1];
            int depth = genotype.getDP();
            double qualPerAlleleCount = 0;

            if(alleleCount > 0)
            {
                final String[] qualCounts = genotype.getExtendedAttribute(READ_CONTEXT_QUALITY, 0).toString()
                        .split(LIST_SEPARATOR, -1);

                int qualTotal = 0;
                for(int i = 0; i < 4; ++i)
                {
                    qualTotal += Integer.parseInt(qualCounts[i]);
                }

                qualPerAlleleCount = qualTotal / (double) alleleCount;
            }

            mResultsWriter.writeSampleVariant(patientId, genotype.getSampleName(), variant, alleleCount, depth, qualPerAlleleCount);
        }
    }
}
