package com.hartwig.hmftools.lilac.variant;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.coverage.AlleleCoverage;
import com.hartwig.hmftools.lilac.coverage.ComplexCoverage;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaGene;

public class CopyNumberAssignment
{
    private final Map<String, List<CopyNumberData>> mSampleCopyNumberData;

    public CopyNumberAssignment()
    {
        mSampleCopyNumberData = Maps.newHashMap();
    }

    public void loadCopyNumberData(final LilacConfig config)
    {
        if(config.CopyNumberFile.isEmpty())
            return;

        try
        {

            Predicate<GeneCopyNumber> filterPredicate = x ->
            {
                HlaGene gene = HlaGene.fromString(x.geneName());
                return gene != null && GENE_CACHE.GeneNames.contains(gene);
            };
            List<GeneCopyNumber> hlaGeneCopyNumbers = GeneCopyNumberFile.read(config.CopyNumberFile).stream()
                    .filter(filterPredicate).toList();

            List<CopyNumberData> cnDataList = hlaGeneCopyNumbers.stream()
                    .map(x -> new CopyNumberData(HlaGene.fromString(x.geneName()), x.minCopyNumber(), x.MinMinorAlleleCopyNumber))
                    .collect(Collectors.toList());

            mSampleCopyNumberData.put(config.Sample, cnDataList);
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to read gene copy number file({}): {}", config.CopyNumberFile, e.toString());
        }
    }

    public static List<Double> formEmptyAlleleCopyNumber(final Iterable<HlaAllele> winners)
    {
        List<Double> alleleCN = Lists.newArrayList();
        winners.forEach(x -> alleleCN.add(0.0));
        return alleleCN;
    }

    public void assign(final String sampleId, final ComplexCoverage refCoverage, final ComplexCoverage tumorCoverage,
            final Collection<Double> copyNumbers)
    {
        LL_LOGGER.info("calculating tumor copy number of winning alleles");

        if(refCoverage.getAlleleCoverage().size() != GENE_CACHE.ExpectAlleleCount
                || tumorCoverage.getAlleleCoverage().size() != GENE_CACHE.ExpectAlleleCount)
        {
            return;
        }

        List<CopyNumberData> cnDataList = mSampleCopyNumberData.get(sampleId);

        if(cnDataList == null || cnDataList.isEmpty())
            return;

        copyNumbers.clear();

        for(int index = 0; index < GENE_CACHE.ExpectAlleleCount; index = index + 2)
        {
            AlleleCoverage refCoverage1 = refCoverage.getAlleleCoverage().get(index);
            AlleleCoverage refCoverage2 = refCoverage.getAlleleCoverage().get(index + 1);

            AlleleCoverage tumorCoverage1 = tumorCoverage.getAlleleCoverage().get(index);
            AlleleCoverage tumorCoverage2 = tumorCoverage.getAlleleCoverage().get(index + 1);

            HlaGene gene = refCoverage1.Allele.Gene;
            CopyNumberData cnData = cnDataList.stream().filter(x -> x.Gene == gene).findFirst().orElse(null);

            if(cnData == null)
            {
                LL_LOGGER.warn("missing gene({}) copy number data", gene.toString());
                copyNumbers.add(0.0);
                copyNumbers.add(0.0);
                continue;
            }

            alleleCopyNumber(cnData, copyNumbers, refCoverage1, refCoverage2, tumorCoverage1, tumorCoverage2);
        }
    }

    private static void alleleCopyNumber(
            final CopyNumberData copyNumberData, final Collection<Double> alleleCopyNumbers,
            final AlleleCoverage refCoverage1, final AlleleCoverage refCoverage2,
            final AlleleCoverage tumorCoverage1, final AlleleCoverage tumorCoverage2)
    {
        double minor = copyNumberData.MinMinorAlleleCopyNumber;
        double major = copyNumberData.MinCopyNumber - minor;

        if(refCoverage1.Allele.equals(refCoverage2.Allele))
        {
            alleleCopyNumbers.add(major);
            alleleCopyNumbers.add(minor);
            return;
        }

        double ratio1 = tumorCoverage1.TotalCoverage / refCoverage1.TotalCoverage;
        double ratio2 = tumorCoverage2.TotalCoverage / refCoverage2.TotalCoverage;

        if(ratio1 >= ratio2)
        {
            alleleCopyNumbers.add(major);
            alleleCopyNumbers.add(minor);
        }
        else
        {
            alleleCopyNumbers.add(minor);
            alleleCopyNumbers.add(major);
        }
    }

}
