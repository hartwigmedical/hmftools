package com.hartwig.hmftools.cup.somatics;

import static com.hartwig.hmftools.common.stats.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.utils.MatrixFile.loadMatrixDataFile;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.somatics.RefSomatics.REF_SIG_TYPE_SNV_COUNT;
import static com.hartwig.hmftools.cup.somatics.SomaticSigs.convertSignatureName;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Somaticvariant.SOMATICVARIANT;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.jooq.Record;
import org.jooq.Record8;
import org.jooq.Result;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;
import htsjdk.variant.vcf.VCFCodec;

public class SomaticDataLoader
{
    public static Matrix loadSampleCountsFromFile(final String filename, final Map<String,Integer> sampleCountsIndex)
    {
        if(filename.isEmpty() || !Files.exists(Paths.get(filename)))
            return null;

        return loadMatrixDataFile(filename, sampleCountsIndex, Lists.newArrayList("BucketName"), true);
    }

    public static Matrix loadSampleMatrixData(final String filename, final Map<String,Integer> sampleCountsIndex)
    {
        if(filename.isEmpty() || !Files.exists(Paths.get(filename)))
            return null;

        return loadMatrixDataFile(filename, sampleCountsIndex, null, true);
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
                final String sigName = convertSignatureName(sigAllocation.signature());
                sigContribs.put(sigName, sigAllocation.allocation());
            }
        }

        return true;
    }

    public static boolean loadRefSignaturePercentileData(
            final String filename, final Map<String,Map<String,double[]>> refCancerSigContribs, final Map<String,double[]> refCancerSnvCounts)
    {
        if(filename.isEmpty())
            return false;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            for(final String line : fileData)
            {
                // SampleId,DataType,Pct_0.00 etc
                final String[] items = line.split(DATA_DELIM, -1);
                String cancerType = items[0];

                String dataType = items[1];

                double[] percentileData = new double[PERCENTILE_COUNT];

                int startIndex = 2;

                for(int i = startIndex; i < items.length; ++i)
                {
                    double value = Double.parseDouble(items[i]);
                    percentileData[i - startIndex] = value;
                }

                if(dataType.equals(REF_SIG_TYPE_SNV_COUNT))
                {
                    refCancerSnvCounts.put(cancerType, percentileData);
                }
                else
                {
                    String sigName = dataType;

                    Map<String, double[]> sigContribsMap = refCancerSigContribs.get(cancerType);

                    if(sigContribsMap == null)
                    {
                        sigContribsMap = Maps.newHashMap();
                        refCancerSigContribs.put(cancerType, sigContribsMap);
                    }

                    sigContribsMap.put(sigName, percentileData);
                }
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sig contrib percentile data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    public static Matrix loadRefSampleCounts(
            final String filename, final List<String> refSampleNames, final List<String> ignoreCols)
    {
        if(filename.isEmpty() || !Files.exists(Paths.get(filename)))
            return null;

        return loadMatrixDataFile(filename, refSampleNames, ignoreCols, true);
    }

    public static List<SomaticVariant> loadSomaticVariants(final String sampleId, final DatabaseAccess dbAccess)
    {
        List<SomaticVariant> variants = Lists.newArrayList();

        Result<Record8<String, String, Integer, String, String, String, String, String>>
                result = dbAccess.context()
                .select(SOMATICVARIANT.GENE, SOMATICVARIANT.CHROMOSOME, SOMATICVARIANT.POSITION,
                        SOMATICVARIANT.REF, SOMATICVARIANT.ALT, SOMATICVARIANT.TYPE, SOMATICVARIANT.GENE,
                        SOMATICVARIANT.TRINUCLEOTIDECONTEXT)
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(SOMATICVARIANT.SAMPLEID.eq(sampleId))
                .and(SOMATICVARIANT.TYPE.eq(VariantType.SNP.toString()))
                .orderBy(SOMATICVARIANT.CHROMOSOME, SOMATICVARIANT.POSITION)
                .fetch();

        for(Record record : result)
        {
            variants.add(SomaticVariant.fromRecord(record));
        }

        return variants;
    }

    public static List<SomaticVariant> loadSomaticVariants(final String vcfFile, final List<VariantType> types)
    {
        CompoundFilter filter = new CompoundFilter(true);
        filter.add(new PassingVariantFilter());

        List<SomaticVariant> variants = Lists.newArrayList();

        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false);

            for(VariantContext variantContext : reader.iterator())
            {
                if(!filter.test(variantContext))
                    continue;

                if(types.isEmpty() || types.contains(VariantType.type(variantContext)))
                {
                    variants.add(SomaticVariant.fromContext(variantContext));
                }
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error(" failed to read somatic VCF file({}): {}", vcfFile, e.toString());
        }

        return variants;
    }

    public static boolean loadSigContribsFromCohortFile(final String filename, final Map<String,Map<String,Double>> sampleSigContributions)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            fileData.remove(0);

            for(final String line : fileData)
            {
                // SampleId,SigName,SigContrib,SigPercent
                final String[] items = line.split(DATA_DELIM, -1);
                String sampleId = items[0];
                String sigName = convertSignatureName(items[1]);
                double sigContrib = Double.parseDouble(items[2]);

                Map<String,Double> sigContribs = sampleSigContributions.get(sampleId);

                if(sigContribs == null)
                {
                    sigContribs = Maps.newHashMap();
                    sampleSigContributions.put(sampleId, sigContribs);
                }

                sigContribs.put(sigName, sigContrib);
            }

            CUP_LOGGER.info("loaded {} sample sig-contributions records from file({})", sampleSigContributions.size(), filename);
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sig contribution data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }
}
