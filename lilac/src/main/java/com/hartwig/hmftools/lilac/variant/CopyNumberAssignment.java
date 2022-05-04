package com.hartwig.hmftools.lilac.variant;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.DELIM;
import static com.hartwig.hmftools.lilac.LilacConstants.EXPECTED_ALLELE_COUNT;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_GENES;
import static com.hartwig.hmftools.lilac.LilacConstants.longGeneName;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.coverage.AlleleCoverage;
import com.hartwig.hmftools.lilac.coverage.ComplexCoverage;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class CopyNumberAssignment
{
    private final Map<String,List<CopyNumberData>> mSampleCopyNumberData;

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
            if(config.CopyNumberFile.contains(config.Sample))
            {
                List<GeneCopyNumber> hlaGeneCopyNumbers = GeneCopyNumberFile.read(config.CopyNumberFile).stream()
                        .filter(x -> HLA_GENES.contains(x.geneName())).collect(Collectors.toList());

                List<CopyNumberData> cnDataList = hlaGeneCopyNumbers.stream()
                        .map(x -> new CopyNumberData(x.geneName(), x.minCopyNumber(), x.minMinorAlleleCopyNumber()))
                        .collect(Collectors.toList());

                mSampleCopyNumberData.put(config.Sample, cnDataList);
            }
            else
            {
                // load a cohort file - for now only retain the required sample's data
                final List<String> fileData = Files.readAllLines(new File(config.CopyNumberFile).toPath());
                String header = fileData.get(0);
                Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIM);
                fileData.remove(0); // remove header
                int sampleIndex = fieldsIndexMap.get("SampleId");
                int geneIndex = fieldsIndexMap.get("Gene");
                int minorIndex = fieldsIndexMap.get("MinCopyNumber");
                int minorMinIndex = fieldsIndexMap.get("MinMinorAlleleCopyNumber");

                boolean loadAll = false;

                if(loadAll)
                {
                    List<CopyNumberData> cnDataList = null;
                    String currentSample = "";

                    for(final String line : fileData)
                    {
                        String[] items = line.split(DELIM, -1);

                        String sampleId = items[sampleIndex];

                        if(!currentSample.equals(sampleId))
                        {
                            cnDataList = Lists.newArrayList();
                            currentSample = sampleId;
                            mSampleCopyNumberData.put(sampleId, cnDataList);
                        }

                        cnDataList.add(new CopyNumberData(
                                items[geneIndex], Double.parseDouble(items[minorIndex]), Double.parseDouble(items[minorMinIndex])));
                    }
                }
                else
                {
                    // only load the configured sample
                    List<CopyNumberData> cnDataList = Lists.newArrayList();

                    for(final String line : fileData)
                    {
                        String[] items = line.split(DELIM, -1);

                        String sampleId = items[sampleIndex];

                        if(!sampleId.equals(config.Sample))
                        {
                            if(!cnDataList.isEmpty())
                                break;
                            else
                                continue;
                        }

                        cnDataList.add(new CopyNumberData(
                                items[geneIndex], Double.parseDouble(items[minorIndex]), Double.parseDouble(items[minorMinIndex])));
                    }

                    mSampleCopyNumberData.put(config.Sample, cnDataList);
                }
            }
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to read gene copy number file({}): {}", config.CopyNumberFile, e.toString());
        }
    }

    public static List<Double> formEmptyAlleleCopyNumber(final List<HlaAllele> winners)
    {
        List<Double> alleleCN = Lists.newArrayList();
        winners.forEach(x -> alleleCN.add(0.0));
        return alleleCN;
    }

    public void assign(
            final String sampleId, final List<HlaAllele> winners,
            final ComplexCoverage refCoverage, final ComplexCoverage tumorCoverage, final List<Double> copyNumbers)
    {
        LL_LOGGER.info("calculating tumor copy number of winning alleles");

        if(refCoverage.getAlleleCoverage().size() != EXPECTED_ALLELE_COUNT || tumorCoverage.getAlleleCoverage().size() != EXPECTED_ALLELE_COUNT)
            return;

        List<CopyNumberData> cnDataList = mSampleCopyNumberData.get(sampleId);

        if(cnDataList == null || cnDataList.isEmpty())
            return;

        copyNumbers.clear();

        for(int index = 0; index < EXPECTED_ALLELE_COUNT; index = index + 2)
        {
            AlleleCoverage refCoverage1 = refCoverage.getAlleleCoverage().get(index);
            AlleleCoverage refCoverage2 = refCoverage.getAlleleCoverage().get(index + 1);

            AlleleCoverage tumorCoverage1 = tumorCoverage.getAlleleCoverage().get(index);
            AlleleCoverage tumorCoverage2 = tumorCoverage.getAlleleCoverage().get(index + 1);

            String gene = longGeneName(refCoverage1.Allele.Gene);
            CopyNumberData cnData = cnDataList.stream().filter(x -> x.Gene.equals(gene)).findFirst().orElse(null);

            if(cnData == null)
            {
                LL_LOGGER.warn("missing gene({}) copy number data", gene);
                copyNumbers.add(0.0);
                copyNumbers.add(0.0);
                continue;
            }

            alleleCopyNumber(cnData, copyNumbers, refCoverage1, refCoverage2, tumorCoverage1, tumorCoverage2);
        }
    }

    private void alleleCopyNumber(
            final CopyNumberData copyNumberData, final List<Double> alleleCopyNumbers,
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
