package com.hartwig.hmftools.cup.sigs;

import static com.hartwig.hmftools.common.sigs.SigUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.sigs.RefSomatics.populateRefPercentileData;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.vcf.VCFCodec;

public class SomaticDataLoader
{
    public static SigMatrix loadSampleCountsFromFile(final String filename, final Map<String,Integer> sampleCountsIndex)
    {
        if(filename.isEmpty())
            return null;

        SigMatrix sampleCounts = loadMatrixDataFile(filename, sampleCountsIndex, Lists.newArrayList("BucketName"));
        sampleCounts.cacheTranspose();

        return sampleCounts;
    }

    public static SigMatrix loadSamplePosFreqFromFile(final String filename, final Map<String,Integer> sampleCountsIndex)
    {
        if(filename.isEmpty())
            return null;

        SigMatrix sampleCounts = loadMatrixDataFile(filename, sampleCountsIndex, null);
        sampleCounts.cacheTranspose();

        return sampleCounts;
    }

    public static boolean loadSigContribsFromDatabase(
            final DatabaseAccess dbAccess, final List<String> sampleIds, final Map<String,Map<String,Double>> sampleSigContributions)
    {
        if(dbAccess == null)
            return false;

        for(final String sampleId : sampleIds)
        {
            Map<String,Double> sigContribs = sampleSigContributions.get(sampleId);

            if(sigContribs == null)
            {
                sigContribs = Maps.newHashMap();
                sampleSigContributions.put(sampleId, sigContribs);
            }

            final List<SignatureAllocation> sigAllocations = dbAccess.readSignatureAllocations(sampleId);

            for(final SignatureAllocation sigAllocation : sigAllocations)
            {
                sigContribs.put(sigAllocation.signature(), sigAllocation.allocation());
            }
        }

        return true;
    }

    public static boolean loadSigContribsFromCohortFile(final String filename, final Map<String,Map<String,Double>> sampleSigContributions)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            for(final String line : fileData)
            {
                // SampleId,SigName,SigContrib,SigPercent
                final String[] items = line.split(DATA_DELIM, -1);
                String sampleId = items[0];
                String sigName = items[1];
                double sigContrib = Double.parseDouble(items[2]);

                Map<String,Double> sigContribs = sampleSigContributions.get(sampleId);

                if(sigContribs == null)
                {
                    sigContribs = Maps.newHashMap();
                    sampleSigContributions.put(sampleId, sigContribs);
                }

                sigContribs.put(sigName, sigContrib);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sig contribution data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    public static boolean loadRefSignaturePercentileData(
            final String filename, final Map<String,Map<String,double[]>> refCancerSigContribs, final Map<String,double[]> refCancerSnvCounts)
    {
        if(filename.isEmpty())
            return true;

        return populateRefPercentileData(filename, refCancerSigContribs, refCancerSnvCounts);
    }

    public static SigMatrix loadRefSampleCounts(final String filename, final List<String> refSampleNames)
    {
        if(filename.isEmpty())
            return null;

        SigMatrix refSampleCounts = loadMatrixDataFile(filename, refSampleNames);
        refSampleCounts.cacheTranspose();
        return refSampleCounts;
    }

    public static SigMatrix loadRefSnvPosFrequences(final String filename, final List<String> refCancerTypes)
    {
        if(filename.isEmpty())
            return null;

        SigMatrix refSnvPosFrequences = loadMatrixDataFile(filename, refCancerTypes);
        refSnvPosFrequences.cacheTranspose();
        return refSnvPosFrequences;
    }

    public static List<SomaticVariant> loadSomaticVariants(final String sampleId, final DatabaseAccess dbAccess)
    {
        final List<SomaticVariant> somaticVariants = dbAccess.readSomaticVariants(sampleId, VariantType.SNP);
        somaticVariants.addAll(dbAccess.readSomaticVariants(sampleId, VariantType.INDEL));
        return somaticVariants;
    }

    public static List<SomaticVariant> loadSomaticVariants(final String sampleId, final String vcfFile, final List<VariantType> types)
    {
        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());

        SomaticVariantFactory variantFactory = new SomaticVariantFactory(filter);
        final List<SomaticVariant> variantList = Lists.newArrayList();

        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false);

            for (VariantContext variant : reader.iterator())
            {
                if (filter.test(variant))
                {
                    final SomaticVariant somaticVariant = variantFactory.createVariant(sampleId, variant).orElse(null);

                    if(somaticVariant == null)
                        continue;

                    if(!types.isEmpty() && !types.contains(variant.getType()))
                        continue;

                    variantList.add(somaticVariant);
                }
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error(" failed to read somatic VCF file({}): {}", vcfFile, e.toString());
        }

        return variantList;
    }
}
